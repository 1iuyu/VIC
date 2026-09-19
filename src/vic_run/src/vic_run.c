/******************************************************************************
* @section DESCRIPTION
*
* This subroutine controls the model core, it solves both the energy and water
* balance models, as well as frozen soils. #include <vic_run.h>
******************************************************************************/

#include "vic_run.h"

/******************************************************************************
* @brief        This subroutine controls the model core, it solves both the
*               energy and water balance models, as well as frozen soils.
******************************************************************************/
int
vic_run(all_vars_struct     *all_vars,
        global_param_struct *gp,
        soil_con_struct     *soil_con,
        veg_con_struct      *veg_con,
        veg_lib_struct      *veg_lib)
{
    extern option_struct     options;
    extern parameters_struct param;

    size_t     iveg;
    int        ErrorFlag;
    double     Tfoliage;
    double     gauge_correction[2];

    cell_data_struct  *cell;
    veg_var_struct    *veg_var;
    energy_bal_struct *energy;
    snow_data_struct  *snow;
    force_data_struct *force;

    /* Set number of vegetation tiles */
    size_t Nveg = veg_con[0].vegetat_type_num;

    /** Set global step times **/
    double step_dt = gp->step_dt;
    double Canopy_Lower = 0.0;
    double Canopy_Upper = 0.0;

    /**************************************************
       Solve Energy and Water Balance for Each
       Vegetation Tile
    **************************************************/
    for (iveg = 0; iveg <= Nveg; iveg++) {
        
        /** Solve Veg Tile only if Coverage Greater than 0% **/
        if (veg_con[iveg].Cv > 0.0) {

            /** Define vegetation class number **/
            size_t veg_class = veg_con[iveg].veg_class;

            /**************************************************
               Initialize Model Parameters
            **************************************************/
            /* Set local pointers */
            veg_var = &(all_vars->veg_var[iveg]);
            cell    = &(all_vars->cell[iveg]);
            snow    = &(all_vars->snow[iveg]);
            energy  = &(all_vars->energy[iveg]);
            force   = &(all_vars->force[iveg]);

            /* Compute gauge undercatch correction factors
            - this assumes that the gauge is free of vegetation effects, so gauge
            correction is constant for the entire grid cell */
            if (options.CORRPREC && force->prec > 0) {
                correct_precip(gauge_correction, force->wind,
                               param.REF_HEIGHT_WIND,
                               param.SOIL_ROUGH, param.SNOW_ROUGH);
            }
            else {
                gauge_correction[0] = 1;
                gauge_correction[1] = 1;
            }           

            // 校正降雨和降雪
            if (options.CORRPREC) {
                force->rainf *= gauge_correction[0];
                force->snowf *= gauge_correction[1];
            }
            // 计算新雪密度
            snow->new_snow_density = new_snow_density(force->air_temp);
           
            /* Calculate the snow and rain interception */
            if (cell->IS_VEG) {

                // initialize canopy terms
                Tfoliage = energy->Tfoliage;

                /** Assign Canopy **/
                Canopy_Upper = veg_lib[veg_class].Canopy_Upper;
                Canopy_Lower = veg_lib[veg_class].Canopy_Lower;

                /* Calculate net LAI and SAI */
                calc_net_veg(Canopy_Upper, 
                             Canopy_Lower,
                             snow, veg_var);

                ErrorFlag = snow_intercept(step_dt, Tfoliage,
                                           force, energy, 
                                           cell, snow, 
                                           soil_con, veg_var);

                if (ErrorFlag == ERROR) {
                    return (ERROR);
                }
            }
            // 设置土地类型指示器
            if (veg_var->fcanopy > 0.0 && veg_var->NetLAI + 
                                          veg_var->NetSAI > 0.0) {
                cell->IS_VEG = true;
            }
            else {
                cell->IS_VEG = false;
            }
            
            /* Initialize snow coverage */
            calc_snow_coverage(step_dt, force->snowf, cell,
                               snow, soil_con, SNOW_ACCU);

            // 初始化粗糙度
            initialize_roughness(Canopy_Upper,
                                 snow->coverage,
                                 cell, veg_var);

            /******************************
              Solve ground surface fluxes
            ******************************/
            ErrorFlag = surface_fluxes(step_dt, force, 
                                       energy, cell, 
                                       snow, soil_con, 
                                       veg_var,
                                       &veg_lib[veg_class]);

            if (ErrorFlag == ERROR) {
                return (ERROR);
            }
        } /** end non-zero area veg tile **/
    } /** end of vegetation loop **/

    return (0);
}