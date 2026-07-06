/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine computes the slope of the liquid water-temperature curve
 * (dθdT) as well as the partial derivative of T*dθdT.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief  
 * Compute the slope of the liquid water-temperature curve.
 *****************************************************************************/
double
water_curve_deriv(size_t  		   nidx,
                  double           T,
                  double           liq,
                  double           matric,
                  soil_con_struct *soil_con)
{
    double tmp_deriv = 0.0;
    double liq_derivT = 0.0;

    // 计算dθ/dψ
    tmp_deriv = SoilWaterRetentionCurve(DERIV_FLAG, nidx, liq,
                                        matric, soil_con);
    if (tmp_deriv > 0.0) {
        liq_derivT = (CONST_LATICE * tmp_deriv) / (CONST_G * T);
    }
    else {
        liq_derivT = 0.0;
    }
    return (liq_derivT);
}