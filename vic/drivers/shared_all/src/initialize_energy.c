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
    
    size_t               i, index;
    size_t               veg;

    // initialize miscellaneous energy balance terms
    for (veg = 0; veg <= Nveg; veg++) {

        // Prognostic states
        for (i = 0; i < MAX_SWBANDS; i++) {
            energy[veg].AlbedoSnowDir[i] = 0.0;
            energy[veg].AlbedoSnowDfs[i] = 0.0;
            energy[veg].AlbedoSoilDir[i] = 0.0;
            energy[veg].AlbedoSoilDfs[i] = 0.0;
            energy[veg].AlbedoGrndDir[i] = 0.0;
            energy[veg].AlbedoGrndDfs[i] = 0.0;
            energy[veg].AlbedoSurfDir[i] = 0.0;
            energy[veg].AlbedoSurfDfs[i] = 0.0;
            energy[veg].ReflGrndDir[i] = 0.0;
            energy[veg].ReflGrndDfs[i] = 0.0;
            energy[veg].ReflSubDir[i] = 0.0;
            energy[veg].ReflSubDfs[i] = 0.0;
            energy[veg].ShortDir2Dir[i] = 0.0;
            energy[veg].ShortDfs2Dir[i] = 0.0;
            energy[veg].ShortDir2Dfs[i] = 0.0;
            energy[veg].ShortDfs2Dfs[i] = 0.0;
            energy[veg].AbsSubDir[i] = 0.0;
            energy[veg].AbsSubDfs[i] = 0.0;
            energy[veg].ReflectVeg[i] = 0.0;
            energy[veg].TransmitVeg[i] = 0.0;
        }
        energy[veg].FrozenOver = false;
        energy[veg].FrozenGrnd = false;
        energy[veg].Nthaw = 0.;
        energy[veg].Nfrost = 0.;
        energy[veg].Tcanopy = 0.0;
        energy[veg].Tsurf = 0.0;
        energy[veg].Tgrnd = 0.0;
        energy[veg].Tair = 0.0;
        energy[veg].Tupper = 0.0;
        energy[veg].Tlower = 0.0;
        energy[veg].fdepth = 0.0;
        energy[veg].tdepth = 0.0;

        for (index = 0; index < MAX_NODES + MAX_SNOWS; index++) {
            energy[veg].Cs_node[index] = 0.0;
            energy[veg].kappa_node[index] = 0.0;
            energy[veg].T[index] = 0.0;
            energy[veg].fusion_fact[index] = 0.0;
            energy[veg].fusion_flux[index] = 0.0;
            energy[veg].kappa_int[index] = 0.0;
            energy[veg].fact[index] = 0.0;
            energy[veg].fn[index] = 0.0;
        }
        // Fluxes
        energy[veg].advection = 0.0;
        energy[veg].AdvectOver = 0.0;
        energy[veg].AdvectGrnd = 0.0;
        energy[veg].AdvectSub = 0.0;
        energy[veg].grnd_flux = 0.0;
        energy[veg].GroundSub = 0.0;
        energy[veg].GroundGrnd = 0.0;
        energy[veg].latent = 0.0;
        energy[veg].LatentSub = 0.0;
        energy[veg].LatentGrnd = 0.0;
        energy[veg].sensible = 0.0;
        energy[veg].SensibleSub = 0.0;
        energy[veg].SensibleGrnd = 0.0;
        energy[veg].SensibleOver = 0.0;
        energy[veg].LatentCanopy = 0.0;
        energy[veg].LatentTransp = 0.0;
        energy[veg].LatentVapGrnd = 0.0;
        energy[veg].LatentVapOver = 0.0;

        energy[veg].EmissLongSub = 0.0;
        energy[veg].EmissLongGrnd = 0.0;
        energy[veg].EmissLongSurf = 0.0;
        energy[veg].ReflShortSurf = 0.0;
        energy[veg].ReflShortGrnd = 0.0;
        energy[veg].ReflShortSub = 0.0;
        
        energy[veg].NetLongSurf = 0.0;
        energy[veg].NetLongGrnd = 0.0;
        energy[veg].NetLongOver = 0.0;
        energy[veg].NetLongSub = 0.0;
        energy[veg].shortwave = 0.0;
        energy[veg].NetShortSurf = 0.0;
        energy[veg].NetShortGrnd = 0.0;
        energy[veg].NetShortSub = 0.0;
        energy[veg].NetShortSoil = 0.0;
        energy[veg].NetShortSnow = 0.0;
    }
}