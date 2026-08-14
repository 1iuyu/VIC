/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate vegetation parameters.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    This subroutine estimates the displacement height of vegetation
 *           with a given average height based on equations on page 4.12 of the
 *           Handbook of Hydrology.
 *****************************************************************************/
double
calc_veg_displacement(double V,
                      double ratio,
                      double canopy_Upper)
{

    double                   value;

    value = ratio * canopy_Upper * V;

    return (value);
}


/******************************************************************************
 * @brief    This subroutine estimates the roughness height of vegetation with
 *           a given average height based on equations on page 4.12 of the
 *           Handbook of Hydrology.
 *****************************************************************************/
void
calc_net_veg(double            Canopy_Upper,
             double            Canopy_Lower,
             snow_data_struct *snow,
             veg_var_struct   *veg_var)
{
    extern parameters_struct param;

    double bury_depth = 0.0;
    double bury_frac = 0.0;
    double snow_depth = snow->snow_depth;
    double coverage = snow->coverage;
    double LAI = veg_var->LAI;
    double SAI = veg_var->SAI;
    // adjust lai and sai for burying by snow.
    // adjust for grasses, crops
    if (Canopy_Upper > 0.0 && Canopy_Upper <= 1.0) {
        bury_frac = 1.0 - (max(min(snow_depth, max(0.05, Canopy_Upper * 0.8)), 
                                            0.0) / (max(0.05, Canopy_Upper * 0.8)));
    }
    else {
        bury_depth = min(max(snow_depth - Canopy_Lower, 0.0),
                                          (Canopy_Upper - Canopy_Lower));
        bury_frac = 1.0 - bury_depth / max(param.TOL_A, (Canopy_Upper - Canopy_Lower));
    }
    // 
    double NetLAI = max(LAI * (1.0 - coverage) + LAI * bury_frac * coverage, 0.0);
    double NetSAI = max(SAI * (1.0 - coverage) + SAI * bury_frac * coverage, 0.0);
    // if exposed lai and sai
    // are less than 0.05, set equal to zero to prevent numerical
    // problems associated with very small lai and sai
    if (NetLAI < 0.05) {
        veg_var->NetLAI = 0.0;
    }
    else {
        veg_var->NetLAI = NetLAI;
    }
    if (NetSAI < 0.05) {
        veg_var->NetSAI = 0.0;
    }
    else {
        veg_var->NetSAI = NetSAI;
    }
}