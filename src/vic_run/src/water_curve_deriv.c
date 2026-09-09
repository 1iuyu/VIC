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

/******************************************************************************
 * @brief  
 * Compute the slope of the snow layer liquid water-temperature curve dθ/dT.
 *****************************************************************************/
double 
snow_curve_deriv(double T,
                 double theta_ice,
                 double theta_liq,
                 double Cs_node)
{
    const double DT_PC = 5.0;
    const double halfDT = 0.5 * DT_PC;   // 半区间宽度

    if (theta_ice == 0.0 && theta_liq == 0.0) {
        return Cs_node;
    }

    double dtheta_liq_dT = 0.0;

    if (T > CONST_TKFRZ) {
        if (theta_ice > 0.0) {
            dtheta_liq_dT = (CONST_RHOICE / CONST_RHOFW) * theta_ice / halfDT;
        }
    }
    else if (T < CONST_TKFRZ) {
        if (theta_liq > 0.0) {
            dtheta_liq_dT = theta_liq / halfDT;
        }
    }
    else { // T == TKFRZ
        if (theta_ice > 0.0) {
            dtheta_liq_dT = (CONST_RHOICE / CONST_RHOFW) * theta_ice / halfDT;
        } else if (theta_liq > 0.0) {
            dtheta_liq_dT = theta_liq / halfDT;
        }
    }

    return Cs_node + CONST_RHOFW * CONST_LATICE * dtheta_liq_dT;
}