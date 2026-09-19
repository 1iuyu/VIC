/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine creates an array of structures that contain information about a
 * cell's states and fluxes.
 *****************************************************************************/

#include "vic_driver_shared_all.h"

/******************************************************************************
 * @brief    Creates an array of structures that contain information about a
 *           cell's states and fluxes.
 *****************************************************************************/
all_vars_struct
make_all_vars(size_t    Nvegs,
              size_t    Nbands)
{
    all_vars_struct all_vars;

    all_vars.snow = calloc(Nvegs, sizeof(*all_vars.snow));
    check_alloc_status(all_vars.snow, "Memory allocation error.");

    all_vars.energy = calloc(Nvegs, sizeof(*all_vars.energy));
    check_alloc_status(all_vars.energy, "Memory allocation error.");

    all_vars.veg_var = calloc(Nvegs, sizeof(*all_vars.veg_var));
    check_alloc_status(all_vars.veg_var, "Memory allocation error.");

    all_vars.cell = calloc(Nvegs, sizeof(*all_vars.cell));
    check_alloc_status(all_vars.cell, "Memory allocation error.");

    all_vars.force = calloc(Nbands, sizeof(*all_vars.force));   
    check_alloc_status(all_vars.force, "Memory allocation error.");

    return (all_vars);
}
