/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine defines the relationship between the temperature and the
 * fraction of liquid water in snow for the freezing curve.
 *  [1] Calculates the fraction of liquid water based on the temperature.
 *  [2] Calculates the temperature based on the fraction of liquid water
 *          (inverse of the freezing curve).
 *  [3] Calculates the derivative of the fraction of liquid water with
 *          respect to the temperature.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief  This subroutine defines the relationship between the temperature 
 *         and the fraction of liquid water for the freezing curve.
 *****************************************************************************/
double
SnowFreezingCurve(int    flag,
                  double Tk,
                  double theta_liq)
{
    extern parameters_struct param;
    double Tdep, Tdim;

    if (flag == SNOWLIQ_FLAG) {
        Tdep = CONST_TKFRZ - min(Tk, CONST_TKFRZ);
        return 1.0 / (1.0 + pow(param.SNOW_FRZ_SCALE * Tdep, 2.0));
    }
    else if (flag == TEMPLIQ_FLAG) {
        return CONST_TKFRZ -
               sqrt((1.0 / theta_liq - 1.0) / (param.SNOW_FRZ_SCALE * param.SNOW_FRZ_SCALE));
    }
    else if (flag == DERILIQ_FLAG) {
        Tdep = CONST_TKFRZ - min(Tk, CONST_TKFRZ);
        Tdim = param.SNOW_FRZ_SCALE * Tdep;
        return (param.SNOW_FRZ_SCALE * 2.0 * Tdim) /
               pow(1.0 + Tdim * Tdim, 2.0);
    }
    else {
        log_err("Unknown FLAG option");
    }

    return 0.0;
}