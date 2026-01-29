/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes the model state (energy balance, water balance, and
 * snow components) to default values.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    Initialize the model state (energy balance, water balance, and
 *           snow components) to default values.
 *****************************************************************************/
void
generate_default_state(all_vars_struct *all_vars,
                       soil_con_struct *soil_con,
                       veg_con_struct  *veg_con)
{
    extern option_struct     options;
    extern global_param_struct global_param;

    size_t                   Nveg;
    size_t                   veg;
    size_t                   band;
    size_t                   lidx;
    size_t                   k, layer;
    double                   Cv, Bexp;
    cell_data_struct        *cell;
    snow_data_struct        *snow;
    energy_bal_struct       *energy;
    veg_var_struct          *veg_var;

    cell = all_vars->cell;
    energy = all_vars->energy;
    snow = all_vars->snow;
    veg_var = all_vars->veg_var;
    Nveg = veg_con[0].vegetat_type_num;

    /*********************************
       Initialize snowpack temperatures
    *********************************/
    for (veg = 0; veg <= Nveg; veg++) {
        Cv = veg_con[veg].Cv;
        if (Cv > 0) {
            snow[veg].swq = 10.0;   // [mm]
            snow[veg].old_swq = 10.0;   // [mm]
            snow[veg].density = 250.;   // [kg/m3]
            if (snow[veg].density > 0.) {
                snow[veg].snow_depth = snow[veg].swq /
                                    snow[veg].density;
            }
            // set snow layer properties
            distribute_snow_state(veg_con[veg].IS_GLAC, 
                                  global_param.dt, &snow[veg]);
        }
    }

    /******************************************************************************
       Initialize soil moistures
       currently setting moist to porosity as default and eliminate init_moist 
       (require user to use a state file if they want control over initial moist)
    ******************************************************************************/
    for (veg = 0; veg <= Nveg; veg++) {
        Cv = veg_con[veg].Cv;
        if (Cv > 0) {
            /* Initialize soil moistures */
            for (lidx = 0; lidx < options.Nnode; lidx++) {
                cell[veg].moist[lidx] =
                        soil_con->Wfc_node[lidx];
                if (cell[veg].moist[lidx] >
                    soil_con->Wsat_node[lidx]) {
                    cell[veg].moist[lidx] =
                            soil_con->Wsat_node[lidx];
                }
            }
        }
    }

    /*********************************
       Initialize soil temperatures
    *********************************/
    for (veg = 0; veg <= Nveg; veg++) {
        Cv = veg_con[veg].Cv;
        if (Cv > 0) {
            band = veg_con[veg].BandIndex;
            /* Initialize soil node temperatures */
            for (k = 0; k < options.Nnode; k++) {
                if (options.FROZEN_SOIL) {
                    cell[veg].soil_T[k] = 250.0;
                }
                else {
                    cell[veg].soil_T[k] = CONST_TKFRZ;
                }
            }
            /* Initial estimate of temperatures */
            energy[veg].Tgrnd = 250.0;
            energy[veg].Tair = 250.0 + soil_con->Tfactor[band];
            energy[veg].Tcanopy = 250.0 + soil_con->Tfactor[band];
        }
    }

    return;
}
