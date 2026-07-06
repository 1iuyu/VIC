/******************************************************************************
 * @section DESCRIPTION
 *
 * Initialize energy structure.
 *****************************************************************************/

#include "vic_driver_shared_all.h"

/******************************************************************************
 * @brief    Initialize energy structure.
 *****************************************************************************/
void
initialize_energy(energy_bal_struct *energy,
                  size_t             Nveg)
{

    size_t               i, lidx;
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
        for (i = 0; i < MAX_CANOPYS; i++) {
            energy[veg].AbsDirSun[i] = 0.0;
            energy[veg].AbsDirSha[i] = 0.0;
            energy[veg].AbsDfsSun[i] = 0.0;
            energy[veg].AbsDfsSha[i] = 0.0;
        }
        energy[veg].FrozenOver = false;
        energy[veg].FrozenGrnd = false;
        energy[veg].Tcanopy = 0.0;
        energy[veg].Tgrnd = 0.0;
        energy[veg].Tsurf = 0.0;
        energy[veg].Tfoliage = 0.0;
        energy[veg].Tstem = 0.0;
        for (lidx = 0; lidx < MAX_NODES; lidx++) {
            energy[veg].Cs_node[lidx] = 0.0;
            energy[veg].last_Cs[lidx] = 0.0;
            energy[veg].kappa_node[lidx] = 0.0;
            energy[veg].T[lidx] = 0.0;
            energy[veg].last_T[lidx] = 0.0;
            energy[veg].kappa_int[lidx] = 0.0;
        }
        energy[veg].energy_flag = false;
        energy[veg].moist_flag = false;
        // Fluxes
        energy[veg].advection = 0.0;
        energy[veg].AdvectOver = 0.0;
        energy[veg].AdvectGrnd = 0.0;
        energy[veg].AdvectSub = 0.0;
        // 地表热通量
        energy[veg].grnd_flux = 0.0;
        // 潜热通量
        energy[veg].latent = 0.0;
        energy[veg].LatentLeaf = 0.0;
        energy[veg].LatentVapGrnd = 0.0;
        energy[veg].LatentVapOver = 0.0;
        // 感热通量
        energy[veg].sensible = 0.0;
        energy[veg].SensibleLeaf = 0.0;
        energy[veg].SensibleStem = 0.0;
        // 辐射项
        energy[veg].EmissLongSub = 0.0;
        energy[veg].EmissLongGrnd = 0.0;
        energy[veg].EmissLongSurf = 0.0;
        energy[veg].ReflShortSurf = 0.0;
        energy[veg].ReflShortGrnd = 0.0;
        energy[veg].ReflShortSub = 0.0;
        // 长波辐射项
        energy[veg].longwave = 0.0;
        energy[veg].NetLongOver = 0.0;
        // 短波辐射项
        energy[veg].shortwave = 0.0;
        energy[veg].NetShortSurf = 0.0;
        energy[veg].NetShortGrnd = 0.0;
        energy[veg].NetShortSub = 0.0;
        energy[veg].NetShortSoil = 0.0;
        energy[veg].NetShortSnow = 0.0;
    }
}
