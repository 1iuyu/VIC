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
    double liq_deriv = 0.0;

    // 计算基质势对未冻水含量的导数 d(ψ_m)/dθ_l
    tmp_deriv = SoilWaterRetentionCurve(DERIV_FLAG, nidx, liq,
                                        matric, soil_con);
    if (tmp_deriv > 0.0) {
        liq_deriv = tmp_deriv =(CONST_LATICE / T) / (CONST_G / tmp_deriv);
    }
    else {
        liq_deriv = 0.0;
    }
    return (liq_deriv);
}