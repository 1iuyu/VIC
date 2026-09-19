/******************************************************************************
 * @section DESCRIPTION
 *
 * Determines from the air temperature what fraction of incoming precipitation
 * is frozen and unfrozen (snow and rain).
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    Determines from the air temperature what fraction of incoming
 *           precipitation is frozen and unfrozen (snow and rain).
 *****************************************************************************/
void
calc_rainonly(double  prec,
              double  Twet_bulb,
              double  rel_humid,
              double  elevation,
              double *snowfall,
              double *rainfall)
{   
    double      snowfrac;
    double      expr;
    double      T0;
    double      delta_T;
    double      delta_S;
    double      MIN_RAIN_TEMP;
    double      MAX_SNOW_TEMP;

    /* Initialize outputs */
    snowfrac = 0.;
    (*snowfall) = 0.;
    (*rainfall) = 0.;

    /* Ding et al. 2014, JH Eq.7.8.9 */
    delta_T = 0.215 - 0.099 * rel_humid + 1.018 * rel_humid * rel_humid;
    delta_S = 2.374 - 1.634 * rel_humid;
    T0 = 267.28 - 0.1042 * elevation + 0.0885 * elevation * elevation + 
                      16.06 * rel_humid - 9.614 * rel_humid * rel_humid;

    if (delta_S != 0.0 && delta_T / delta_S > log(2.0)) {
        expr = exp(delta_T / delta_S) - 2.0 * exp(-delta_T / delta_S);
        // avoid log(<=0)
        if (expr > 1e-12) {
            MIN_RAIN_TEMP = T0 - delta_S * log(expr);
            MAX_SNOW_TEMP = 2.0 * T0 - MIN_RAIN_TEMP;
        } else {
            MIN_RAIN_TEMP = T0;
            MAX_SNOW_TEMP = T0;
        }
    }
    else {
        MIN_RAIN_TEMP = T0;
        MAX_SNOW_TEMP = T0;     
    }

    /* Decide phase fraction */
    if (Twet_bulb < MAX_SNOW_TEMP && Twet_bulb > MIN_RAIN_TEMP) {
        snowfrac = 1.0 / (1.0 + exp((Twet_bulb - T0) / delta_S));
    }
    else if (Twet_bulb >= MAX_SNOW_TEMP) {
        snowfrac = 0.0;  // all rain
    }
    else if (Twet_bulb <= MIN_RAIN_TEMP) {
        snowfrac = 1.0;  // all snow
    }

    /* Final outputs */ 
    *snowfall = prec * snowfrac;
    *rainfall = prec * (1 - snowfrac);

}
