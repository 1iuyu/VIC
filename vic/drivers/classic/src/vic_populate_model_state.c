/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes the model state (energy balance, water balance, and
 * snow components).
 *
 * If a state file is provided to the model then its
 * contents are checked to see if it agrees with the current simulation set-up,
 * if so it is used to initialize the model state.  If no state file is
 * provided the model initializes all variables with defaults and the user
 * should expect to throw out the beginning of the simulation period as model
 * start-up.
 *****************************************************************************/

#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    Initialize the model state (energy balance, water balance, and
 *           snow components).
 *****************************************************************************/
void
vic_populate_model_state(all_vars_struct *all_vars,
                         filep_struct     filep,
                         size_t           cellnum,
                         soil_con_struct *soil_con,
                         veg_con_struct  *veg_con)
{
    extern option_struct options;

    size_t               Nveg;

    cell_data_struct   *cell;
    energy_bal_struct  *energy;
    snow_data_struct   *snow;
    veg_var_struct     *veg_var;

    cell = all_vars->cell;
    energy = all_vars->energy;
    snow = all_vars->snow;
    veg_var = all_vars->veg_var;

    Nveg = veg_con[0].vegetat_type_num;

    // Initialize all data structures to 0
    initialize_soil(cell, Nveg);
    initialize_snow(snow, Nveg);
    initialize_veg(veg_var, Nveg);
    initialize_energy(energy, Nveg);

    // Read initial state from a file if provided
    if (options.INIT_STATE) {
        read_initial_model_state(filep.init_state, all_vars, Nveg,
                                 options.SNOW_BAND, cellnum, soil_con);
    }
    else {
        // else generate a default state
        generate_default_state(all_vars, soil_con, veg_con);
    }

    // compute those state variables that are derived from the others
    compute_derived_state_vars(all_vars, soil_con, veg_con);
}
