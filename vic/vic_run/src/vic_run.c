/******************************************************************************
* @section DESCRIPTION
*
* This subroutine controls the model core, it solves both the energy and water
* balance models, as well as frozen soils. #include <vic_run.h>
******************************************************************************/

#include "vic_driver_shared_all.h"

/******************************************************************************
* @brief        This subroutine controls the model core, it solves both the
*               energy and water balance models, as well as frozen soils.
******************************************************************************/
int
vic_run(force_data_struct   *force,
        all_vars_struct     *all_vars,
        global_param_struct *gp,
        soil_con_struct     *soil_con,
        veg_con_struct      *veg_con,
        veg_lib_struct      *veg_lib)
{
    extern option_struct     options;
    extern param_set_struct  param_set;
    extern parameters_struct param;

    size_t     iveg;
    size_t     Nveg;
    size_t     veg_class;
    size_t     band;
    int        ErrorFlag;
    double     step_dt;
    double     Tfoliage;
    double     Tair; 
    double     step_prec;
    double     wind;
    double     Cv;
    double     gauge_correction[2];

    cell_data_struct  *cell;
    veg_var_struct    *veg_var;
    energy_bal_struct *energy;
    snow_data_struct  *snow;

    /* Set number of vegetation tiles */
    Nveg = veg_con[0].vegetat_type_num;

    /** Set global step times **/
    step_dt = gp->step_dt;
    double rainfall = 0.0;
    double snowfall = 0.0;
    double Canopy_Lower = 0.0;
    double Canopy_Upper = 0.0;
    /* Compute gauge undercatch correction factors
       - this assumes that the gauge is free of vegetation effects, so gauge
       correction is constant for the entire grid cell */
    if (options.CORRPREC && force->prec[NR] > 0) {
        correct_precip(gauge_correction, force->wind[NR], gp->wind_h,
                       param.SOIL_ROUGH, param.SNOW_ROUGH);
    }
    else {
        gauge_correction[0] = 1;
        gauge_correction[1] = 1;
    }

    /**************************************************
       Solve Energy and Water Balance for Each
       Vegetation Tile
    **************************************************/
    for (iveg = 0; iveg <= Nveg; iveg++) {
        
        /** Solve Veg Tile only if Coverage Greater than 0% **/
        if (veg_con[iveg].Cv > 0.0) {

            /** Define vegetation Coverage **/
            Cv = veg_con[iveg].Cv;

            /** Define vegetation class number **/
            veg_class = veg_con[iveg].veg_class;

            /**************************************************
               Initialize Model Parameters
            **************************************************/
            /* Set local pointers */
            veg_var = &(all_vars->veg_var[iveg]);
            cell    = &(all_vars->cell[iveg]);
            snow    = &(all_vars->snow[iveg]);
            energy  = &(all_vars->energy[iveg]);

            /**************************************************
               Define vegetation elevation bands
            **************************************************/
            band = veg_con[iveg].BandIndex;
            Tair = force->air_temp[NR] + soil_con->Tfactor[band];
            int has_prec = param_set.TYPE[PREC].SUPPLIED;
            int has_rain = param_set.TYPE[RAINF].SUPPLIED;
            int has_snow = param_set.TYPE[SNOWF].SUPPLIED;
            // 根据输入的气象数据计算降雨和降雪
            if (has_prec && !has_rain && !has_snow) {
                /* set air temperature and precipitation for this snow band */
                step_prec = force->prec[NR] * soil_con->Pfactor[band];

                /** Calculate Fraction of Precipitation that falls as Rain **/
                calc_rainonly(Tair, step_prec, force->vp[NR],
                                    force->pressure[NR],
                                    force->rel_humid[NR],
                                    soil_con->elevation,
                                    &snowfall, &rainfall);
                force->rainf[NR] = rainfall;
                force->snowf[NR] = snowfall;
            }
            else if (has_rain && has_snow) {
                force->rainf[NR] *= soil_con->Pfactor[band];
                force->snowf[NR] *= soil_con->Pfactor[band];
                snowfall = force->snowf[NR];
                rainfall = force->rainf[NR];
            }
            // 校正降雨和降雪
            snowfall *= gauge_correction[SNOW];
            rainfall *= gauge_correction[RAIN];

            // 计算新雪密度
            if (snowfall > 0.0) {
                snow->new_snow_density = new_snow_density(Tair);
            }
            else {
                snow->new_snow_density = 0.0;
            }
           
            /* Calculate the snow and rain interception */
            if (veg_con[iveg].IS_GLAC == false) {

                // initialize canopy terms
                Tfoliage = energy->Tfoliage;
            
                /* Initialize wind speeds */
                wind = force->wind[NR];

                /** Assign Canopy **/
                Canopy_Upper = veg_lib[veg_class].Canopy_Upper;
                Canopy_Lower = veg_lib[veg_class].Canopy_Lower;

                // 分配子网格根区分数
                cell->Nroot = veg_con[iveg].Nroot;

                /* Calculate net LAI and SAI */
                calc_net_veg(Canopy_Upper, 
                             Canopy_Lower,
                             snow->snow_depth,
                             veg_var);

                ErrorFlag = snow_intercept(step_dt, Tfoliage,
                                           &snowfall, &rainfall,
                                           wind, snow, veg_var);
                if (ErrorFlag == ERROR) {
                    return (ERROR);
                }
            }
            // 设置土地类型指示器
            if (veg_var->fcanopy > 0.0 && veg_var->NetLAI + veg_var->NetSAI > 0.0) {
                cell->IS_VEG = true;
            }
            else {
                cell->IS_VEG = false;
            }
            cell->IS_GLAC = veg_con[iveg].IS_GLAC;
            
            /* Initialize snow coverage */
            calc_snow_coverage(Cv, veg_con[iveg].IS_GLAC,
                               snowfall * step_dt,  // (mm H2O)
                               snow, soil_con);

            // 初始化粗糙度
            initialize_roughness(veg_con[iveg].IS_GLAC, 
                                 Canopy_Upper,
                                 snow->coverage,
                                 veg_con[iveg].root,
                                 cell, veg_var);

            /******************************
              Solve ground surface fluxes
            ******************************/
            if (veg_con[iveg].IS_GLAC == false) {   // 非冰川HRU
                ErrorFlag = surface_fluxes(step_dt, Tair, 
                                           snowfall, rainfall,
                                           force, energy,
                                           cell, snow,
                                           soil_con, veg_var,
                                           &veg_lib[veg_class]);

            }
            else {  // 冰川HRU
                ErrorFlag = surface_fluxes_glac(Tair, snowfall,
                                                rainfall,
                                                force, energy, 
                                                gp, cell,
                                                snow, soil_con);
            }

            if (ErrorFlag == ERROR) {
                return (ERROR);
            }
        } /** end non-zero area veg tile **/
    } /** end of vegetation loop **/

    return (0);
}
