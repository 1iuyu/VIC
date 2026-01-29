/******************************************************************************
 * @section DESCRIPTION
 *
 * Determines from the air temperature what fraction of incoming precipitation
 * is frozen and unfrozen (snow and rain).
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    Determines from the air temperature what fraction of incoming
 *           precipitation is frozen and unfrozen (snow and rain).
 *****************************************************************************/
void
calc_rainonly(double  air_temp,
              double  prec,
              double  vp,
              double  pressure,
              double  rel_humid,
              double  elevation,
              double *snowfall,
              double *rainfall)
{   
    extern parameters_struct param;
    
    size_t      iter;
    double      snowfrac;
    double      svp;         
    double      Twet;        // wetbulb temperature ℃
    double      latent;      // latent heat of vapor or sublimation
    double      PsychConst;  // psychrometric constant
    double      delta;
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
    // 转换为摄氏温度 [K->℃]
    air_temp = K_TO_C(air_temp);
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

    if (air_temp > 0.) {
        latent = CONST_LATVAP;
    }
    else {
        latent = CONST_LATSUB;
    }
    /* Calculate psychrometric constant */
    PsychConst = CONST_CPDAIR * pressure / (0.622 * latent);

    /* initial guess */
    Twet = air_temp - 3.0;

    /* iterate wet-bulb temperature */
    iter = 0;
    do {
        // saturation vapor pressure
        svp = 610.8 * exp((17.27 * Twet) / (237.3 + Twet));
        delta = (svp - vp) / PsychConst;
        Twet -= delta;
        iter++;
    } 
    while (iter < param.MAX_ITER_WETBULB &&
        fabs(delta) > param.TOL_WETBULB);

    /* Decide phase fraction */
    if (Twet < MAX_SNOW_TEMP && Twet > MIN_RAIN_TEMP) {
        snowfrac = 1.0 / (1.0 + exp((Twet - T0) / delta_S));
    }
    else if (Twet >= MAX_SNOW_TEMP) {
        snowfrac = 0.;  // all rain
    }
    else if (Twet <= MIN_RAIN_TEMP) {
        snowfrac = 1.;  // all snow
    }

    /* Final outputs */ 
    *snowfall = prec * snowfrac;
    *rainfall = prec * (1 - snowfrac);

}
