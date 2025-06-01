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
double
calc_rainonly(double            air_temp,
              double            prec,
              double            vp,
              double            vpd,
              double            atmos_pressure,
              soil_con_struct  *soil_con)
{
    extern option_struct        options;

    double                      rainonly;
    double                      deltaT; 
    double                      deltaS;
    double                      T0;
    double                      Tw;
    double                      esat;
    double                      d_esat;
    double                      MIN_RAIN_TEMP; 
    double                      MAX_SNOW_TEMP;
    double                      rel_humid;
    double                      elevation;

    rainonly  = 0.;
    rel_humid = vp / (vp + vpd);
    elevation = soil_con->elevation;

    if (options.TEMP_TH_TYPE == CLASSIC) {

        MAX_SNOW_TEMP = soil_con->MAX_SNOW_TEMP;
        MIN_RAIN_TEMP = soil_con->MIN_RAIN_TEMP;

        if (MAX_SNOW_TEMP <= MIN_RAIN_TEMP) {
            log_err("MAX_SNOW_TEMP must be greater then MIN_RAIN_TEMP");
        }
        if (air_temp < MAX_SNOW_TEMP && air_temp > MIN_RAIN_TEMP) {
            rainonly = (air_temp - MIN_RAIN_TEMP) /
                       (MAX_SNOW_TEMP - MIN_RAIN_TEMP) * prec;
        }
        else if (air_temp >= MAX_SNOW_TEMP) {
            rainonly = prec;
        }
    }

    else if (options.TEMP_TH_TYPE == DING) {
        // 计算deltaT, deltaS和T0
        deltaT = 0.215 - 0.099 * rel_humid + 1.018 * rel_humid * rel_humid;
        deltaS = 2.374 - 1.634 * rel_humid;
        T0 = -5.87 - 0.1042 * elevation + 0.0885 * elevation * elevation + 16.06 * rel_humid - 9.614 * rel_humid * rel_humid;
        esat = svp(air_temp);
        d_esat = (17.67 * 243.5) / powf(air_temp + 243.5, 2) * esat;
        Tw = air_temp - esat * (1 - rel_humid) / (0.000643 * atmos_pressure + d_esat);

        // 计算MIN_RAIN_TEMP和MAX_SNOW_TEMP
        if (deltaT / deltaS > log(2)) {
            double log_argument = exp(deltaT / deltaS) - 2 * exp(-deltaT / deltaS);
            if (log_argument > 0) {
                MIN_RAIN_TEMP = T0 - deltaS * log(log_argument);
                MAX_SNOW_TEMP = 2 * T0 - MIN_RAIN_TEMP;
            } else {
                // 处理log_argument为负的情况
                log_err("log argument is non-positive");
            }
        } 
        else {
            MIN_RAIN_TEMP = T0;
            MAX_SNOW_TEMP = T0;
        }

        if (MAX_SNOW_TEMP <= MIN_RAIN_TEMP) {
            log_err("MAX_SNOW_TEMP must be greater then MIN_RAIN_TEMP");
        }
        if (Tw < MAX_SNOW_TEMP && Tw > MIN_RAIN_TEMP) {
            rainonly = (Tw - MIN_RAIN_TEMP) /
                       (MAX_SNOW_TEMP - MIN_RAIN_TEMP) * prec;
        }
        else if (Tw >= MAX_SNOW_TEMP) {
            rainonly = prec;
        }
    }
    return(rainonly);
}
