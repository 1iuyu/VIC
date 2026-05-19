/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initailizes the vegetation variable array.
 *****************************************************************************/

#include "vic_driver_shared_all.h"

/******************************************************************************
 * @brief    This routine initailizes the vegetation variable array.
 *****************************************************************************/
void
initialize_veg(veg_var_struct *veg_var,
               size_t          Nveg)
{
    for (size_t i = 0; i < Nveg; i++) {
        // Prognostic states
        veg_var[i].fcanopy = 0.0;
        veg_var[i].LAI = 0.0;
        veg_var[i].SAI = 0.0;
        veg_var[i].NetLAI = 0.0;
        veg_var[i].NetSAI = 0.0;
        veg_var[i].int_snow = 0.0;
        veg_var[i].int_rain = 0.0;
        veg_var[i].canopy_swq = 0.0;
        veg_var[i].wetFrac = 0.0;
        veg_var[i].dryFrac = 0.0;
        // Diagnostic states
        veg_var[i].MaxSnowInt = 0.0;
        veg_var[i].MaxRainInt = 0.0;
        veg_var[i].Wdew = 0.1;
        // Fluxes
        veg_var[i].RainThroughFall = 0.0;
        veg_var[i].SnowThroughFall = 0.0;
        veg_var[i].RainDrip = 0.0;
        veg_var[i].SnowDrip = 0.0;
        veg_var[i].SnowUnload = 0.0;
        veg_var[i].leaf_sun = 0.0;
        veg_var[i].leaf_shade = 0.0;
        veg_var[i].RS_sunlit = 0.0;
        veg_var[i].RS_shade = 0.0;

    }
}
