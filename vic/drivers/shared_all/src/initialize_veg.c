/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initailizes the vegetation variable array.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    This routine initailizes the vegetation variable array.
 *****************************************************************************/
void
initialize_veg(veg_var_struct *veg_var,
               size_t          Nveg)
{
    extern option_struct options;

    size_t               i, j;

    for (i = 0; i < Nveg; i++) {
        // Prognostic states
        veg_var[i].albedo = 0.0;
        veg_var[i].displacement = 0.0;
        veg_var[i].fcanopy = 0.0;
        veg_var[i].LAI = 0.0;
        veg_var[i].SAI = 0.0;
        veg_var[i].NetLAI = 0.0;
        veg_var[i].NetSAI = 0.0;
        veg_var[i].int_snow = 0.0;
        veg_var[i].int_rain = 0.0;
        veg_var[i].wetFrac = 0.0;
        veg_var[i].MaxSnowInt = 0.0;
        veg_var[i].MaxRainInt = 0.0;
        veg_var[i].roughness = 0.0;
        veg_var[i].Wdew = 0.1;
        // Fluxes
        veg_var[i].RainThroughFall = 0.0;
        veg_var[i].SnowThroughFall = 0.0;
        veg_var[i].RainDrip = 0.0;
        veg_var[i].SnowDrip = 0.0;
    }
    if (options.CARBON) {

        // Carbon states
        veg_var[i].AnnualNPP = 0.0;
        veg_var[i].AnnualNPPPrev = 0.0;
        veg_var[i].Ci = 0.0;
        veg_var[i].NPPfactor = 0.0;
        veg_var[i].rc = 0.0;
        for (j = 0; j < options.Ncanopy; j++) {
            veg_var[i].CiLayer[j] = 0.0;
            veg_var[i].NscaleFactor[j] = 0.0;
            veg_var[i].rsLayer[j] = 0.0;
        }
        // Carbon fluxes
        veg_var[i].aPAR = 0.0;
        for (j = 0; j < options.Ncanopy; j++) {
            veg_var[i].aPARLayer[j] = 0.0;
        }
        veg_var[i].GPP = 0.0;
        veg_var[i].Litterfall = 0.0;
        veg_var[i].NPP = 0.0;
        veg_var[i].Raut = 0.0;
        veg_var[i].Rdark = 0.0;
        veg_var[i].Rgrowth = 0.0;
        veg_var[i].Rmaint = 0.0;
        veg_var[i].Rphoto = 0.0;
    }
}
