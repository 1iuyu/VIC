/******************************************************************************
* @section DESCRIPTION
*
* This subroutine updates data structures with values for the current
* time step.
******************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
* @brief        This subroutine controls the model core, it solves both the
*               energy and water balance models, as well as frozen soils.
******************************************************************************/
int
update_step_vars(all_vars_struct *all_vars,
                 veg_con_struct  *veg_con,
                 veg_hist_struct *veg_hist)
{
    extern option_struct options;

    unsigned short       iveg;
    size_t               Nveg;
    veg_var_struct      *veg_var;

    /* set local pointers */
    veg_var = all_vars->veg_var;

    /* Set number of vegetation types */
    Nveg = veg_con[0].vegetat_type_num;

    /* Assign current veg characteristics */
    for (iveg = 0; iveg <= Nveg; iveg++) {
        veg_var[iveg].albedo = veg_hist[iveg].albedo[NR];
        veg_var[iveg].displacement = veg_hist[iveg].displacement[NR];
        veg_var[iveg].fcanopy = veg_hist[iveg].fcanopy[NR];
        veg_var[iveg].LAI = veg_hist[iveg].LAI[NR];
        veg_var[iveg].roughness = veg_hist[iveg].roughness[NR];
    }

    return (0);
}
