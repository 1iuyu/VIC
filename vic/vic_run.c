/******************************************************************************
* @section DESCRIPTION
*
* This subroutine controls the model core, it solves both the energy and water
* balance models, as well as frozen soils. #include <vic_run.h>
******************************************************************************/

#include <vic_driver_shared_all.h>

veg_lib_struct *vic_run_veg_lib;

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

    unsigned short     iveg;
    size_t             Nveg;
    unsigned short     veg_class;
    unsigned short     band;
    int                ErrorFlag;
    double             step_dt;
    double             Tcanopy;
    double             Tair; 
    double             step_prec;
    double             Canopy_Lower;
    double             Canopy_Upper;
    double             wind;
    double             displacement[2];
    double             roughness[2];
    double             ref_height;
    double             Cv;
    double             gauge_correction[2];

    cell_data_struct  *cell;
    veg_var_struct    *veg_var;
    energy_bal_struct *energy;
    snow_data_struct  *snow;

    // assign vic_run_veg_lib to veg_lib, so that the veg_lib for the correct
    // grid cell is used within vic_run. For simplicity sake, use vic_run_veg_lib
    // everywhere within vic_run
    vic_run_veg_lib = veg_lib;

    /* Set number of vegetation tiles */
    Nveg = veg_con[0].vegetat_type_num;

    /** Set global step times **/
    step_dt = gp->dt;
    double rainfall = 0.;
    double snowfall = 0.;
    /* Compute gauge undercatch correction factors
       - this assumes that the gauge is free of vegetation effects, so gauge
       correction is constant for the entire grid cell */
    if (options.CORRPREC && force->prec[NR] > 0) {
        correct_precip(gauge_correction, force->wind[NR], gp->wind_h,
                       soil_con->rough, soil_con->snow_rough);
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
            else {
                if (has_rain && has_snow) {
                    force->rainf[NR] *= soil_con->Pfactor[band];
                    force->snowf[NR] *= soil_con->Pfactor[band];
                    snowfall = force->snowf[NR];
                    rainfall = force->rainf[NR];
                }
            }
            // 校正降雨和降雪
            snowfall *= gauge_correction[SNOW];
            rainfall *= gauge_correction[RAIN];

            // 计算新雪密度
            snow->new_snow_density = new_snow_density(Tair);

            // initialize canopy terms
            Tcanopy = energy->Tcanopy;
            cell->VPcanopy = force->vp[NR];
            
            /* Initialize wind speeds */
            wind = force->wind[NR];

            /** Assign Canopy **/
            Canopy_Upper = vic_run_veg_lib[veg_class].Canopy_Upper;
            Canopy_Lower = vic_run_veg_lib[veg_class].Canopy_Lower;

            /* Calculate the snow and rain interception */
            if (veg_class != options.GLACIER_ID) {
                /* Calculate net LAI and SAI */
                calc_net_veg(Canopy_Upper, Canopy_Lower,
                             snow->snow_depth, veg_var);

                ErrorFlag = snow_intercept(step_dt, Tcanopy,
                                           &snowfall, &rainfall,
                                           wind, snow, veg_var);
                if (ErrorFlag == ERROR) {
                    return (ERROR);
                }
            }
            
            /* Initialize snow coverage */
            calc_snow_coverage(Cv, veg_class, snow, soil_con);

            /* Set surface descriptive variables [0]grnd, [1] veg */
            if (veg_class != options.GLACIER_ID) {
                // Ground surface parameters
                roughness[0] = soil_con->rough * (1 - snow->coverage) +
                                    soil_con->snow_rough * snow->coverage;
                displacement[0] = snow->snow_depth;
                // Vegetation surface parameters
                roughness[1] = veg_var->roughness;
                displacement[1] = calc_veg_displacement(snow->snow_depth,
                                                        Canopy_Upper);
                /* Estimate reference height */
                if (veg_var->NetLAI + veg_var->NetSAI > 0.) {
                    ref_height = displacement[1] + param.REF_HEIGHT;
                }
                else {
                    ref_height = displacement[0] + param.REF_HEIGHT;
                }
                if (ref_height <= displacement[0]) {
                    ref_height = displacement[0] + param.REF_HEIGHT;
                }
            }
            else {
                /* Glacier surface descriptive variables */
                roughness[0] = soil_con->snow_rough;
                displacement[0] = snow->snow_depth;
                // For glaciers, both surfaces use the same parameters
                roughness[1] = roughness[0];
                displacement[1] = displacement[0];
                ref_height = displacement[1] + param.REF_HEIGHT;
            }

            // 分配子网格根区分数
            options.Nroot = veg_con->Nroot;

            /* Soil and ice thermal properties for the layer */
            prepare_full_energy(veg_class, step_dt, cell,
                                energy, snow, soil_con);

            /******************************
              Solve ground surface fluxes
            ******************************/
            if (veg_class != options.GLACIER_ID) {
                ErrorFlag = surface_fluxes(Tair, snowfall, 
                                           rainfall,
                                           ref_height, roughness,
                                           displacement,
                                           force, energy, gp, 
                                           cell, snow,
                                           soil_con, veg_var,
                                           &vic_run_veg_lib[veg_class]);

            }
            else {
                ErrorFlag = surface_fluxes_glac(Tair, snowfall,
                                                rainfall, 
                                                ref_height,
                                                roughness,
                                                displacement,
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
