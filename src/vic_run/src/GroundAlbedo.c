/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate ground albedo based on soil and snow albedo.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
* @brief    This subroutine computes the ground surface albedo.'
*
* @note     Computes albedo as a function of snow age and season, based on the
*           algorithm of the US Army Corps of Engineers.
******************************************************************************/
void
GroundAlbedo(double             moist,
             double             coverage,
             energy_bal_struct *energy,
             soil_con_struct   *soil_con)
{
    extern option_struct options;
    size_t      i; 
    double      swcf_albedo;
    double *AlbedoSat = soil_con->AlbedoSat;
    double *AlbedoDry = soil_con->AlbedoDry;
    double *AlbedoSoilDir = energy->AlbedoSoilDir;
    double *AlbedoSoilDfs = energy->AlbedoSoilDfs;
    double *AlbedoGrndDir = energy->AlbedoGrndDir;
    double *AlbedoGrndDfs = energy->AlbedoGrndDfs;
    double *AlbedoSnowDir = energy->AlbedoSnowDir;
    double *AlbedoSnowDfs = energy->AlbedoSnowDfs;

    for (i = 0; i < options.Nswband; i++) {

        swcf_albedo = max(0.11 - 0.40 * moist, 0.0);
        AlbedoSoilDir[i] = min(AlbedoSat[i] +
                                swcf_albedo, AlbedoDry[i]);
        
        AlbedoSoilDfs[i] = AlbedoSoilDir[i];
        AlbedoGrndDir[i] = AlbedoSoilDir[i] * (1.0 - coverage) + 
                                AlbedoSnowDir[i] * coverage;
        AlbedoGrndDfs[i] = AlbedoSoilDfs[i] * (1.0 - coverage) + 
                                AlbedoSnowDfs[i] * coverage;
    }

}