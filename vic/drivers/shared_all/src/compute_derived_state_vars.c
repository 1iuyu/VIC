/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine computes the state variables (energy balance, water balance,
 * and snow components) that are derived from the variables that are stored in
 * state files.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    Compute the state variables (energy balance, water balance,
 *           and snow components) that are derived from the variables that
 *           are stored in state files.
 *****************************************************************************/
void
compute_derived_state_vars(all_vars_struct *all_vars,
                           soil_con_struct *soil_con,
                           veg_con_struct  *veg_con)
{
    extern global_param_struct global_param;
    extern option_struct       options;

    size_t                     Nveg;
    size_t                     veg;
    int                        ErrorFlag;
    double                     Cv;
    double                     dt_thresh;

    cell_data_struct          *cell;
    energy_bal_struct         *energy;
    snow_data_struct          *snow;

    cell = all_vars->cell;
    energy = all_vars->energy;
    snow = all_vars->snow;
    Nveg = veg_con[0].vegetat_type_num;

    /******************************************
       Compute soil thermal node properties
    ******************************************/
    for (veg = 0; veg <= Nveg; veg++) {
        // Initialize soil for existing vegetation types
        Cv = veg_con[veg].Cv;
        bool IS_GLAC = veg_con[veg].IS_GLAC;
        if (Cv > 0) {
            // set soil moisture properties for all soil thermal nodes
            if (options.FROZEN_SOIL) {
                ErrorFlag =
                        distribute_node_moisture_properties(
                            IS_GLAC,
                            cell[veg].moist,
                            cell[veg].ice,
                            cell[veg].liq,
                            cell[veg].soil_T,
                            soil_con->Wsat_node,
                            soil_con->expt_node,
                            soil_con->bubble_node,
                            soil_con->FS_ACTIVE);
                if (ErrorFlag == ERROR) {
                    log_err("Error setting physical properties for "
                            "soil thermal nodes");
                }
            }

            // Check node spacing v time step
            // (note this is only approximate since heat capacity and
            // conductivity can change considerably during the
            // simulation depending on soil moisture and ice content)
            if (options.FROZEN_SOIL) {
                    // in seconds
                    dt_thresh = 0.5 * energy[veg].Cs_node[1] /
                                energy[veg].kappa_node[1] *
                                pow((soil_con->dz_soil[1]),
                                    2);
                    if (global_param.dt > dt_thresh) {
                        log_err("You are currently running FROZEN SOIL "
                                "with an explicit method (IMPLICIT is "
                                "set to FALSE).  For the explicit method "
                                "to be stable, time step %f seconds is too "
                                "large for the given thermal node spacing "
                                "%f m, soil heat capacity %f J/m3/K, and "
                                "soil thermal conductivity %f J/m/s/K.  "
                                "Either set IMPLICIT to TRUE in your "
                                "global parameter file (this is the "
                                "recommended action), or decrease time "
                                "step length to <= %f seconds, or decrease "
                                "the number of soil thermal nodes.",
                                global_param.dt,
                                soil_con->dz_soil[1],
                                energy[veg].Cs_node[1],
                                energy[veg].kappa_node[1], dt_thresh);
                    }
            }

            /* Find freezing and thawing front depths */
            if (soil_con->FS_ACTIVE) {
                    find_0_degree_fronts(&energy[veg],
                                         soil_con->Zsum_soil,
                                         cell[veg].soil_T,
                                         options.Nnode);
            }
        }
    }
}
