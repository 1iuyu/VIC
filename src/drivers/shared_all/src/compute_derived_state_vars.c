/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine computes the state variables (energy balance, water balance,
 * and snow components) that are derived from the variables that are stored in
 * state files.
 *****************************************************************************/

#include "vic_driver_shared_all.h"

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

    size_t veg, Nveg;
    double Cv;
    double dt_thresh;

    energy_bal_struct *energy;
    energy = all_vars->energy;
    Nveg = veg_con[0].vegetat_type_num;

    /******************************************
       Compute soil thermal node properties
    ******************************************/
    for (veg = 0; veg <= Nveg; veg++) {
        // Initialize soil for existing vegetation types
        Cv = veg_con[veg].Cv;
        if (Cv > 0) {
            // Check node spacing v time step
            // (note this is only approximate since heat capacity and
            // conductivity can change considerably during the
            // simulation depending on soil moisture and ice content)
            // in seconds
            dt_thresh = 0.5 * energy[veg].Cs_node[0] /
                        energy[veg].kappa_node[0] *
                        pow((soil_con->dz_soil[0]), 2.0);

            if (global_param.step_dt > dt_thresh) {
                log_warn("Crank-Nicolson stability condition is violated. "
                            "The current dimensionless parameter "
                            "exceeds the critical threshold of 0.5. "
                            "This may cause spurious oscillations or unphysical "
                            "decay in the numerical solution. "
                            "Now switch to the implicit method "
                            "Current settings: time step = %f s, "
                            "first layer node spacing = %f m, "
                            "soil heat capacity = %f J/m3/K, "
                            "soil thermal conductivity = %f J/m/s/K. "
                            "Recommended maximum stable time step: %f s. "
                            "To resolve this: Reduce the time step to <= %f s,",
                            global_param.step_dt,
                            soil_con->dz_soil[0],
                            energy[veg].Cs_node[0],
                            energy[veg].kappa_node[0],
                            dt_thresh, dt_thresh);
            }
        }
    }
}
