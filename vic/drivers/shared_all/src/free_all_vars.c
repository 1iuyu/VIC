/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine frees all memory allocated down the all_vars data structure.
 *
 * This include all grid cell specific variables (soil, vegetation, energy,
 * snow).
 *****************************************************************************/

#include "vic_driver_shared_all.h"

/******************************************************************************
 * @brief    Free all variables.
 *****************************************************************************/
void
free_all_vars(all_vars_struct *all_vars)
{

    free(all_vars->cell);
    free(all_vars->veg_var);
    free(all_vars->energy);
    free(all_vars->snow);
}
