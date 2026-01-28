/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate vegetation parameters.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    This subroutine estimates the displacement height of vegetation
 *           with a given average height based on equations on page 4.12 of the
 *           Handbook of Hydrology.
 *****************************************************************************/
double
calc_veg_displacement(double snow_depth,
                      double canopy_Upper)
{

    double                   value;

    value = max((0.65 * canopy_Upper), snow_depth);

    return (value);
}


/******************************************************************************
 * @brief    This subroutine estimates the roughness height of vegetation with
 *           a given average height based on equations on page 4.12 of the
 *           Handbook of Hydrology.
 *****************************************************************************/
void
calc_net_veg(double          Canopy_Upper,
             double          Canopy_Lower,
             double          snow_depth,
             veg_var_struct *veg_var)
{
    extern parameters_struct param;

    double      bury_depth;
    double      bury_frac;

    if (Canopy_Upper > 0.0 && Canopy_Upper <= 1.0) {
        bury_depth = Canopy_Upper * exp(-min(snow_depth, 10.0) / 0.2);
        bury_frac = min(snow_depth, bury_depth) / bury_depth;
    }
    else {
        bury_depth = min(max(snow_depth - Canopy_Lower, 0.0),
                                          (Canopy_Upper - Canopy_Lower));
        bury_frac = bury_depth / max(param.TOL_A, (Canopy_Upper - Canopy_Lower));
    }

    veg_var->NetLAI = veg_var->LAI * (1 - bury_frac);
    veg_var->NetSAI = veg_var->SAI * (1 - bury_frac);

}