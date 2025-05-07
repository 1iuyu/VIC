/******************************************************************************
 * @section DESCRIPTION
 *
 * Initialize energy structure.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    Initialize energy structure.
 *****************************************************************************/
void
initialize_energy(energy_bal_struct *energy,
                  size_t             Nveg)
{
    extern option_struct options;
    
    size_t               index;
    size_t               veg;

    // initialize miscellaneous energy balance terms
    for (veg = 0; veg <= Nveg; veg++) {

            // Prognostic states
            energy[veg].AlbedoLake = 0.0;
            energy[veg].AlbedoOver = 0.0;
            energy[veg].AlbedoUnder = 0.0;
            energy[veg].Cs[0] = 0.0;
            energy[veg].Cs[1] = 0.0;
            energy[veg].frozen = false;
            energy[veg].kappa[0] = 0.0;
            energy[veg].kappa[1] = 0.0;
            energy[veg].Nfrost = 0;
            energy[veg].Nthaw = 0;
            energy[veg].T1_index = 0;
            energy[veg].Tcanopy = 0;
            energy[veg].Tcanopy_fbflag = false;
            energy[veg].Tcanopy_fbcount = 0;
            energy[veg].Tfoliage = 0.0;
            energy[veg].Tfoliage_fbflag = false;
            energy[veg].Tfoliage_fbcount = 0;
            energy[veg].Tsurf = 0.0;
            energy[veg].Tsurf_fbflag = false;
            energy[veg].Tsurf_fbcount = 0;
            energy[veg].unfrozen = 0.0;
            for (index = 0; index < options.Nnode - 1; index++) {
                energy[veg].Cs_node[index] = 0.0;
                energy[veg].ice[index] = 0.0;
                energy[veg].kappa_node[index] = 0.0;
                energy[veg].moist[index] = 0.0;
                energy[veg].T[index] = 0.0;
                energy[veg].T_fbflag[index] = false;
                energy[veg].T_fbcount[index] = 0;
            }
            for (index = 0; index < MAX_FRONTS - 1; index++) {
                energy[veg].fdepth[index] = 0.0;
                energy[veg].tdepth[index] = 0.0;
            }
            // Fluxes
            energy[veg].advected_sensible = 0.0;
            energy[veg].advection = 0.0;
            energy[veg].AtmosError = 0.0;
            energy[veg].AtmosLatent = 0.0;
            energy[veg].AtmosLatentSub = 0.0;
            energy[veg].AtmosSensible = 0.0;
            energy[veg].canopy_advection = 0.0;
            energy[veg].canopy_latent = 0.0;
            energy[veg].canopy_latent_sub = 0.0;
            energy[veg].canopy_refreeze = 0.0;
            energy[veg].canopy_sensible = 0.0;
            energy[veg].deltaCC = 0.0;
            energy[veg].deltaH = 0.0;
            energy[veg].error = 0.0;
            energy[veg].fusion = 0.0;
            energy[veg].grnd_flux = 0.0;
            energy[veg].latent = 0.0;
            energy[veg].latent_sub = 0.0;
            energy[veg].longwave = 0.0;
            energy[veg].LongOverIn = 0.0;
            energy[veg].LongUnderIn = 0.0;
            energy[veg].LongUnderOut = 0.0;
            energy[veg].melt_energy = 0.0;
            energy[veg].NetLongAtmos = 0.0;
            energy[veg].NetLongOver = 0.0;
            energy[veg].NetLongUnder = 0.0;
            energy[veg].NetShortAtmos = 0.0;
            energy[veg].NetShortGrnd = 0.0;
            energy[veg].NetShortOver = 0.0;
            energy[veg].NetShortUnder = 0.0;
            energy[veg].out_long_canopy = 0.0;
            energy[veg].out_long_surface = 0.0;
            energy[veg].refreeze_energy = 0.0;
            energy[veg].sensible = 0.0;
            energy[veg].shortwave = 0.0;
            energy[veg].ShortOverIn = 0.0;
            energy[veg].ShortUnderIn = 0.0;
            energy[veg].snow_flux = 0.0;
            energy[veg].glacier_flux = 0.0;
            energy[veg].deltaCC_glac = 0.0;
            energy[veg].glacier_melt_energy = 0.0;
    }
}
