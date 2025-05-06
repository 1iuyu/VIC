/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine frees all memory allocated down the all_vars data structure.
 *
 * This include all grid cell specific variables (soil, vegetation, energy,
 * snow).
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    Free all variables.
 *****************************************************************************/
void
free_all_vars(all_vars_struct *all_vars,
              int              Nveg)
{
    extern option_struct options;

    int                  j, Nitems;

    Nitems = Nveg + 1;

    free((char *) all_vars->cell);
    for (j = 0; j < Nitems; j++) {
        if (options.CARBON) {
            free((char *) all_vars->veg_var[j].NscaleFactor);
            free((char *) all_vars->veg_var[j].aPARLayer);
            free((char *) all_vars->veg_var[j].CiLayer);
            free((char *) all_vars->veg_var[j].rsLayer);
        }
    }
    free((char *) all_vars->veg_var);

    free((char *) all_vars->energy);

    free((char *) all_vars->snow);

    free((char *) all_vars->glacier);
}
