/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate ground albedo based on soil and snow albedo.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
* @brief    This subroutine computes the ground surface albedo.
*
* @note     Computes ground albedo as a function of soil and snow albedo. 
******************************************************************************/
void
GroundAlbedo(double             moist,
             double             coszen,
             double             coverage,
             energy_bal_struct *energy,
             cell_data_struct  *cell,
             soil_con_struct   *soil_con)
{
    extern parameters_struct param;
    extern option_struct     options;
    double  swcf_albedo;
    double *AlbedoSat = soil_con->AlbedoSat;
    double *AlbedoDry = soil_con->AlbedoDry;
    double *AlbedoSoilDir = energy->AlbedoSoilDir;
    double *AlbedoSoilDfs = energy->AlbedoSoilDfs;

    for (size_t i = 0; i < options.Nswband; i++) {
        /* 植被：基于湿度动态计算裸土反照率 */
        if (cell->IS_VEG) {
            swcf_albedo = max(0.11 - 0.40 * moist, 0.0);
            AlbedoSoilDir[i] = min(AlbedoSat[i] +
                                    swcf_albedo, AlbedoDry[i]);
            
            AlbedoSoilDfs[i] = AlbedoSoilDir[i];
        }
        else if (cell->IS_GLAC) {
            /* 冰川：使用固定冰川反照率 */
            AlbedoSoilDir[i] = param.GLAC_ALBEDO[i];
            AlbedoSoilDfs[i] = AlbedoSoilDir[i];           
        }
        else if (cell->IS_WET) {
            // 水体或湿地
            if (energy->Tgrnd > 0.0) {
                AlbedoSoilDir[i] = 0.05 / (max(0.001, coszen) + 0.15);
                AlbedoSoilDfs[i] = 0.1;
            }
            else {
                AlbedoSoilDir[i] = param.LAKE_ALBEDO[i];
                AlbedoSoilDfs[i] = AlbedoSoilDir[i];
            }
        }
    }
}