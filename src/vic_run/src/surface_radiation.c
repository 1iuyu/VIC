/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine calculate the solar fluxes absorbed by vegetation 
 * and ground surface
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief  
 * Compute the solar fluxes absorbed by vegetation and ground surface.
 *****************************************************************************/
int
surface_radiation(double            *shortwave_dir,
                  double            *shortwave_dfs,
                  energy_bal_struct *energy,
                  cell_data_struct  *cell,
                  veg_var_struct    *veg_var)
{
    extern option_struct     options;
    size_t i;
    double NetShortSub = 0.0;
    double NetShortSurf = 0.0;
    double NetShortGrnd = 0.0;
    double NetShortSoil = 0.0;
    double NetShortSnow = 0.0;
    double transmit_dir = 0.0;
    double transmit_dfs = 0.0;
    double tmp_absorb_grnd = 0.0;
    double ShortOverDir[MAX_SWBANDS];
    double ShortOverDfs[MAX_SWBANDS];
    double *aPAR_sun = veg_var->aPAR_sun;
    double *aPAR_sha = veg_var->aPAR_sha;
    double *AbsSubDir = energy->AbsSubDir;
    double *AbsSubDfs = energy->AbsSubDfs;
    double *ShortDir2Dir = energy->ShortDir2Dir;
    double *ShortDfs2Dir = energy->ShortDfs2Dir;
    double *ShortDfs2Dfs = energy->ShortDfs2Dfs;
    double *AlbedoGrndDir = energy->AlbedoGrndDir;
    double *AlbedoGrndDfs = energy->AlbedoGrndDfs;
    double *AlbedoSoilDir = energy->AlbedoSoilDir;
    double *AlbedoSoilDfs = energy->AlbedoSoilDfs;
    double *AlbedoSnowDir = energy->AlbedoSnowDir;
    double *AlbedoSnowDfs = energy->AlbedoSnowDfs;
    double *AlbedoSurfDir = energy->AlbedoSurfDir;
    double *AlbedoSurfDfs = energy->AlbedoSurfDfs;

    for (i = 0; i < options.Nswband; i++) {
        if (cell->IS_VEG) {
            // absorbed by canopy
            ShortOverDir[i] = shortwave_dir[i] * AbsSubDir[i];
            ShortOverDfs[i] = shortwave_dfs[i] * AbsSubDfs[i];
            NetShortSub += ShortOverDir[i] + ShortOverDfs[i];
            // transmitted solar fluxes incident on grnd.
            transmit_dir = shortwave_dir[i] * ShortDir2Dir[i];
            transmit_dfs = shortwave_dir[i] * ShortDfs2Dir[i] +
                                    shortwave_dfs[i] * ShortDfs2Dfs[i];
            // solar radiation absorbed by ground surface.
            tmp_absorb_grnd = transmit_dir * (1.0 - AlbedoGrndDir[i]) + 
                                transmit_dfs * (1.0 - AlbedoGrndDfs[i]);
            NetShortGrnd += tmp_absorb_grnd;

            // calculate absorbed solar by soil/snow separately
            NetShortSoil += transmit_dir * (1.0 - AlbedoSoilDir[i]) + 
                                        transmit_dfs * (1.0 - AlbedoSoilDfs[i]);
            NetShortSnow += transmit_dir * (1.0 - AlbedoSnowDir[i]) + 
                                        transmit_dfs * (1.0 - AlbedoSnowDfs[i]);
        }
        else if (cell->IS_GLAC) {
            NetShortGrnd += shortwave_dir[i] * (1.0 - AlbedoGrndDir[i]) +
                            shortwave_dfs[i] * (1.0 - AlbedoGrndDfs[i]);
        }
    }
    NetShortSurf = NetShortGrnd + NetShortSub;
    energy->NetShortGrnd = NetShortGrnd;
    energy->NetShortSurf = NetShortSurf;
    energy->shortwave = NetShortGrnd;
    energy->NetShortSub = NetShortSub;
    energy->NetShortSoil = NetShortSoil;
    energy->NetShortSnow = NetShortSnow;
    // absorbed PAR for sunlit and shade leaves in canopy layer
    if (cell->IS_VEG) {
        double LAI_sun = 0.0;
        double LAI_sha = 0.0;
        for (i = 0; i < cell->Ncanopy; i++) {
            LAI_sun += veg_var->LAIsun_z[i];
            LAI_sha += veg_var->LAIsha_z[i];
        }
        veg_var->LAI_sun = LAI_sun;
        veg_var->LAI_sha = LAI_sha;
        // Absorbed PAR profile through canopy
        for (i = 0; i < cell->Ncanopy; i++) {
            aPAR_sun[i] = shortwave_dir[0] * energy->AbsDirSun[i] + 
                                        shortwave_dfs[0] * energy->AbsDfsSun[i];
            aPAR_sha[i] = shortwave_dir[0] * energy->AbsDirSha[i] + 
                                        shortwave_dfs[0] * energy->AbsDfsSha[i];
        }
    }

    /* reflected solar radiation */
    double ReflShortSurf = 0.0;
    double ReflShortGrnd= 0.0;
    double ReflShortSub = 0.0;
    for (i = 0; i < options.Nswband; i++) {
        if (cell->IS_VEG) {
            ReflShortSurf += AlbedoSurfDir[i] * shortwave_dir[i] + 
                                AlbedoSurfDfs[i] * shortwave_dfs[i];

            ReflShortGrnd += energy->ReflGrndDir[i] * shortwave_dir[i] + 
                            energy->ReflGrndDfs[i] * shortwave_dfs[i];
            ReflShortSub += energy->ReflSubDir[i] * shortwave_dir[i] +
                            energy->ReflSubDfs[i] * shortwave_dfs[i];
        }
        else if (cell->IS_GLAC) {
            ReflShortSurf += shortwave_dir[i] * AlbedoGrndDir[i] +
                                shortwave_dfs[i] * AlbedoGrndDfs[i];
            ReflShortGrnd = ReflShortSurf;
        }
    }
    energy->ReflShortGrnd = ReflShortGrnd;
    energy->ReflShortSurf = ReflShortSurf;
    energy->ReflShortSub = ReflShortSub;
    
    return (0);
}