/******************************************************************************
* @section DESCRIPTION
*
* Snowpack compaction process
* Update snow depth via compaction due to destructive metamorphism, overburden, & melt.
******************************************************************************/

#include "vic_run.h"

/******************************************************************************
* @brief    Calculate snow depth via compaction due to destructive metamorphism, 
*           overburden, & melt.
******************************************************************************/
void 
snow_compaction(double            step_dt,
                double            air_temp,
                double            pressure,
                snow_data_struct *snow)
{
    extern parameters_struct   param;

    size_t i, Nsnow;                   
    double pack_press = 0.0;                 // Pressure of overlying snow [kg/m2]
    double TempDiff;                   // ConstFreezePoint - [K]
    double SnowVoid;                   // Void (1 - pack_ice - pack_water)
    double SnowMass = 0.0;                   // Water mass (ice + liquid) [kg/m2]
    double SnowAgeFac;                 // rate of compaction due to destructive metamorphism [1/s]
    double SnowAging = 0.0;       // rate of compaction due to destructive metamorphism [1/s]
    double SnowBurden = 0.0;      // rate of compaction of snowpack due to overburden [1/s]
    double SnowMelt = 0.0;        // rate of compaction of snowpack due to melt [1/s]
    double Compact_depth = 0.0;   // change in snow layers-thickness due to compaction [1/s]

    // 定义指针指向结构体中的数组
    double *dz_snow = snow->dz_snow;
    double *pack_T = snow->pack_T;
    double *pack_ice = snow->pack_ice;
    double *pack_liq = snow->pack_liq;
    double *snow_frac = snow->snow_frac;
    double *pack_melt = snow->pack_melt;
    double *theta_liq = snow->theta_liq;
    double *density = snow->density;  // Partial density of ice [kg/m3]
    double *last_snowfrac = snow->last_snowfrac;
    Nsnow = snow->Nsnow;
    double SNOW_COMPACT_P = -0.000695 * air_temp + 0.206067;
    if (pressure >= 85000.0) {
        SNOW_COMPACT_P = max(SNOW_COMPACT_P, 0.017);
    }
    else if (pressure >= 80000.0) {
        SNOW_COMPACT_P = max(SNOW_COMPACT_P, 0.018);
    }
    else {
        SNOW_COMPACT_P = max(SNOW_COMPACT_P, 0.019);
    }
    SNOW_COMPACT_P = min(SNOW_COMPACT_P, 0.0315);

    // Start snow compaction
    for (i = 0; i < Nsnow; i++) {
        SnowMass = pack_ice[i] + pack_liq[i];
        if (SnowMass < param.TOL_A) {
            snow_frac[i] = 0.0;
            continue;
        }
        else {
            snow_frac[i] = pack_ice[i] / SnowMass;
        }
        // 计算雪空隙率
        double ice_volume = pack_ice[i] / CONST_RHOICE;
        double liq_volume = pack_liq[i] / CONST_RHOFW;
        SnowVoid = 1.0 - (ice_volume + liq_volume) / dz_snow[i];

        // Allow compaction only for non-saturated node and higher ice lens node
        if (SnowVoid > 0.001 && pack_ice[i] > 0.1) {
            density[i] = pack_ice[i] / dz_snow[i];
            TempDiff = max(0.0, CONST_TKFRZ - pack_T[i]);

            // Settling/compaction as a result of destructive metamorphism
            SnowAgeFac = exp(-param.SNOW_COMPACT_B * TempDiff);
            SnowAging = -param.SNOW_COMPACT_A * SnowAgeFac;
            if (density[i] > param.SNOW_COMPACT_DM) {
                SnowAging *= exp(-46.0e-3 * (density[i] - 
                                            param.SNOW_COMPACT_DM));
            }
            if (pack_liq[i] > (0.01 * dz_snow[i])) {
                SnowAging *= param.SNOW_COMPACT_C; // Liquid water term
            }

            // Compaction due to overburden
            /* SnowBurden = -(pack_press + 0.5 * SnowMass) * 
                                    exp(-0.08 * TempDiff - SNOW_COMPACT_P * 
                                        density[i]) / param.SNOW_COMPACT_ETA; */
            double f1 = 1.0 / (1.0 + 60.0 * theta_liq[i]);
            double eta = 4.0 * f1 * (density[i] / 358.0) * exp(0.1 * TempDiff + 0.023 * density[i]) * 7622370.0;
            SnowBurden = -(pack_press + 0.5 * SnowMass) / eta;

            // Compaction occurring during melt
            if (pack_melt[i] > 0.0) {
                SnowMelt = max(0.0, (last_snowfrac[i] - snow_frac[i]) / 
                                                max(param.TOL_A, last_snowfrac[i]));
                SnowMelt = -SnowMelt / step_dt;
            } 
            else {
                SnowMelt = 0.0;
            }

            // Time rate of fractional change in snow thickness (units of s-1)
            Compact_depth = (SnowAging + SnowBurden + SnowMelt) * step_dt;
            Compact_depth = min(0.0, max(-0.5, Compact_depth));

            // Change in thickness due to compaction
            dz_snow[i] *= (1.0 + Compact_depth);

            // Constrain thickness to physical limits
            dz_snow[i] = max(dz_snow[i], pack_ice[i] / CONST_RHOICE + 
                                        pack_liq[i] / CONST_RHOFW);

            // Constrain snow density to a reasonable range (50~500 kg/m3)
            dz_snow[i] = min(max(dz_snow[i], (pack_ice[i] + pack_liq[i]) / 500.0), 
                                         (pack_ice[i] + pack_liq[i]) / 50.0);
            // Update snow density
            density[i] = pack_ice[i] / dz_snow[i];
        }
        // Update pressure of overlying snow
        pack_press += SnowMass;
    }
}