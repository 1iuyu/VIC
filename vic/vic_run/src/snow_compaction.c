/******************************************************************************
* @section DESCRIPTION
*
* Snowpack compaction process
* Update snow depth via compaction due to destructive metamorphism, overburden, & melt.
******************************************************************************/

#include <vic_run.h>

/******************************************************************************
* @brief    Calculate snow depth via compaction due to destructive metamorphism, 
*           overburden, & melt.
******************************************************************************/
void 
snow_compaction(double            step_dt,
                double            air_temp,
                double            pressure,
                cell_data_struct *cell,
                snow_data_struct *snow)
{

    extern parameters_struct   param;

    size_t      i;                          // Snow layer loop index
    double      pack_press;                 // Pressure of overlying snow [kg/m2]
    double      TempDiff;                   // ConstFreezePoint - [K]
    double      SnowVoid;                   // Void (1 - pack_ice - pack_water)
    double      old_SnowFrac[MAX_SNOWS];    // ice fraction in snow layers
    double      SnowMass;                   // Water mass (ice + liquid) [kg/m2]
    double      SnowAgeFac;                 // rate of compaction due to destructive metamorphism [1/s]
    double      SnowAging[MAX_SNOWS];       // rate of compaction due to destructive metamorphism [1/s]
    double      SnowBurden[MAX_SNOWS];      // rate of compaction of snowpack due to overburden [1/s]
    double      SnowMelt[MAX_SNOWS];        // rate of compaction of snowpack due to melt [1/s]
    double      Compact_depth[MAX_SNOWS];   // change in snow layers-thickness due to compaction [1/s]
    double      density;                    // Partial density of ice [kg/m3]

    // 定义指针指向结构体中的数组
    double *dz_snow = snow->dz_snow;
    double *pack_T = snow->pack_T;
    double *pack_ice = snow->pack_ice;
    double *pack_liq = snow->pack_liq;
    double *SnowFrac = cell->SnowFrac;
    int *PhaseChange = cell->PhaseChange;

    // Initialization for output variables
    for (i = 0; i < snow->Nsnow; i++) {
        SnowAging[i]  = 0.0;
        SnowBurden[i] = 0.0;
        SnowMelt[i]   = 0.0;
        old_SnowFrac[i] = SnowFrac[i];
    }
    pack_press = 0.0;
    SnowMass   = 0.0;
    double SNOW_COMPACT_P = -0.000695 * air_temp + 0.206067;
    if (pressure >= 85000.0) {
        SNOW_COMPACT_P = max(SNOW_COMPACT_P, 0.017);
    }
    if (pressure >= 80000.0 && pressure < 85000.0) {
        SNOW_COMPACT_P = max(SNOW_COMPACT_P, 0.018);
    }
    if (pressure < 80000.0) {
        SNOW_COMPACT_P = max(SNOW_COMPACT_P, 0.019);
    }
    param.SNOW_COMPACT_P = max(SNOW_COMPACT_P, 0.0315);

    // Start snow compaction
    for (i = 0; i < snow->Nsnow; i++) {
        SnowMass = pack_ice[i] + pack_liq[i];
        SnowFrac[i] = pack_ice[i] / SnowMass;
        cell->SnowFrac[i] = SnowFrac[i];
        // 计算雪空隙率
        double ice_volume = pack_ice[i] / CONST_RHOICE;
        double liq_volume = pack_liq[i] / CONST_RHOFW;
        SnowVoid = 1.0 - (ice_volume + liq_volume) / dz_snow[i];

        // Allow compaction only for non-saturated node and higher ice lens node
        if (SnowVoid > 0.001 && pack_ice[i] > 0.1) {
            density = pack_ice[i] / dz_snow[i];
            TempDiff = max(0.0, CONST_TKFRZ - pack_T[i]);

            // Settling/compaction as a result of destructive metamorphism
            SnowAgeFac = exp(-param.SNOW_COMPACT_B * TempDiff);
            SnowAging[i] = -param.SNOW_COMPACT_A * SnowAgeFac;
            if (density > param.SNOW_COMPACT_DM) {
                SnowAging[i] *= exp(-46.0e-3 * (density - 
                                            param.SNOW_COMPACT_DM));
            }
            if (pack_liq[i] > (0.01 * dz_snow[i])) {
                SnowAging[i] *= param.SNOW_COMPACT_C; // Liquid water term
            }

            // Compaction due to overburden
            SnowBurden[i] = -(pack_press + 0.5 * SnowMass) * 
                                    exp(-0.08 * TempDiff - param.SNOW_COMPACT_P * 
                                        density) / param.SNOW_COMPACT_ETA;

            // Compaction occurring during melt
            if (PhaseChange[i] == 1) {
                SnowMelt[i] = max(0.0, (old_SnowFrac[i] - 
                                                SnowFrac[i]) / 
                                                max(param.TOL_A, old_SnowFrac[i]));
                SnowMelt[i] = -SnowMelt[i] / step_dt;
            } else {
                SnowMelt[i] = 0.0;
            }

            // Time rate of fractional change in snow thickness (units of s^-1)
            Compact_depth[i] = (SnowAging[i] + SnowBurden[i] + SnowMelt[i]) * step_dt;
            Compact_depth[i] = max(-0.5, Compact_depth[i]);

            // Change in thickness due to compaction
            dz_snow[i] *= (1.0 + Compact_depth[i]);

            // Constrain thickness to physical limits
            dz_snow[i] = max(dz_snow[i], pack_ice[i] / CONST_RHOICE + 
                                        pack_liq[i] / CONST_RHOFW);

            // Constrain snow density to a reasonable range (50~500 kg/m3)
            dz_snow[i] = min(max(dz_snow[i], (pack_ice[i] + pack_liq[i]) / 500.0), 
                                         (pack_ice[i] + pack_liq[i]) / 50.0);
            snow->density = density;
        }
        // Update pressure of overlying snow
        pack_press += SnowMass;
    }
}