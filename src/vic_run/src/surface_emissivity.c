/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine calculate Compute ground, vegetation, and total surface 
 * longwave emissivity
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief  
 * Compute ground, vegetation, and total surface longwave emissivity.
 *****************************************************************************/
int
surface_emissivity(energy_bal_struct *energy,
                   cell_data_struct  *cell,
                   snow_data_struct  *snow,
                   veg_var_struct    *veg_var)
{
    extern parameters_struct param;
    /* Initialize variables */
    double coverage = snow->coverage;
    double NetLAI = veg_var->NetLAI;
    double NetSAI = veg_var->NetSAI;
    double fcanopy = veg_var->fcanopy;
    double EmissLongSub = 0.0;
    double EmissLongGrnd = 0.0;
    double EmissLongSurf = 0.0;

    if (cell->IS_VEG) {
        EmissLongSub = 1.0 - exp(-(NetLAI + NetSAI) / 1.0);
        EmissLongGrnd = param.EMISS_GRND * coverage +
                            param.EMISS_ICE * (1.0 - coverage);
    }
    else if (cell->IS_GLAC) {
        EmissLongGrnd = param.EMISS_ICE * (1.0 - coverage) + 
                                        param.EMISS_SNOW * coverage;      
    }
    else if (cell->IS_WET) {
        EmissLongGrnd = param.EMISS_H2O * (1.0 - coverage) + 
                                        param.EMISS_SNOW * coverage;           
    }
    // net surface emissivity
    EmissLongSurf = fcanopy * (EmissLongGrnd * (1.0 - EmissLongSub) + EmissLongSub + 
                        EmissLongSub * (1.0 - EmissLongSub) * 
                            (1.0 - EmissLongGrnd)) + (1.0 - fcanopy) * EmissLongGrnd;

    energy->EmissLongSub = EmissLongSub;
    energy->EmissLongGrnd = EmissLongGrnd;
    energy->EmissLongSurf = EmissLongSurf;

    return (0);
}