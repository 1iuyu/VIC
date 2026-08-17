/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate snowpack hydrology processes (sublimation/frost, evaporation/dew,
 * meltwater) and update snowpack ice and liquid water content.
 *****************************************************************************/

#include "vic_run.h"

/*****************************************************************************
  * @brief   Calculate snowpack hydrology processes (sublimation/frost,
  *          evaporation/dew, meltwater) and update snowpack ice and liquid
  *          water content.
 *****************************************************************************/
int
snow_hydrology(double             step_dt,
               double             air_temp,
               double             rainfall,
               double             pressure,
               double             wind,
               energy_bal_struct *energy,
               cell_data_struct  *cell,
               snow_data_struct  *snow,        
               soil_con_struct   *soil_con)
{
    extern parameters_struct param;
    size_t i, lidx;  
    double vapor_grnd = 0.0;
    double snow_sublim = 0.0;
    double snowfrost = 0.0;
    double snow_evap = 0.0;
    double snow_dew = 0.0;
    double soil_evap = 0.0;
    double soil_dew = 0.0;
    double soil_sublim = 0.0;
    double soilfrost = 0.0;
    // 指针赋值               
    double *liq = cell->liq;
    double *ice = cell->ice;
    double *dz_snow = snow->dz_snow;
    double *porosity = snow->porosity;
    double *theta_liq = snow->theta_liq;
    double *dz_soil = soil_con->dz_soil;
    double *pack_ice = snow->pack_ice;
    double *pack_liq = snow->pack_liq;
    double *pack_frze = snow->pack_frze;
    double *pack_melt = snow->pack_melt;
    double *last_thice = snow->last_thice;
    double *theta_ice = snow->theta_ice;
    double *pack_outflow = snow->pack_outflow;

    /* initialize */
    size_t Nsnow = snow->Nsnow;
    double snow_density = 0.0;
    double excess_flux = 0.0;
    double liquid_capacity = 0.0;
    double coverage = snow->coverage;
    double LatentVapGrnd = energy->LatentVapGrnd;
    // 重制雪层融化量和冻结量
    for (i = 0; i < MAX_SNOWS; i++) {
        pack_frze[i] = 0.0;
        pack_melt[i] = 0.0;
    }
    /** compute soil/snow surface evap,
        dew rate based on energy flux. **/
    // positive part of ground latent heat
    double snow_limit = 0.0, soil_limit = 0.0;
    if (Nsnow > 0) {
        vapor_grnd = energy->LatentSnow / LatentVapGrnd; // mm/s
        if (vapor_grnd >= 0.0) {
            if (pack_ice[0] + pack_liq[0] > 0.0) {
                snow_evap = max(vapor_grnd * (pack_liq[0] / 
                                (pack_ice[0] + pack_liq[0])), 0.0); // mm/s
            }
            else {
                snow_evap = 0.0;
            }
            snow_sublim = vapor_grnd - snow_evap;
        }
        else {
            if (energy->Tgrnd >= CONST_TKFRZ) {
                snowfrost = fabs(vapor_grnd); // mm/s
            }
            else {
                snow_dew = fabs(vapor_grnd); // mm/s
            }
        }
        snow_limit = (pack_liq[0] + pack_ice[0]) / (coverage * step_dt);
        if (vapor_grnd > snow_limit) {
            snow_evap = max(pack_liq[0] / (coverage * step_dt), 0.0);
            snow_sublim = max(pack_ice[0] / (coverage * step_dt), 0.0);
            energy->LatentSnow = snow_limit * LatentVapGrnd;
            energy->SensibleSnow += (vapor_grnd - snow_limit) * LatentVapGrnd;
        }
        // update soil energy fluxes
        vapor_grnd = energy->LatentSoil / LatentVapGrnd; // mm/s
        if (vapor_grnd >= 0.0) {
            if (liq[0] + ice[0] > 0.0) {
                soil_evap = max(vapor_grnd * (liq[0] / (liq[0] + ice[0])), 0.0); // mm/s
            }
            else {
                soil_evap = 0.0;
            }
            soil_sublim = vapor_grnd - soil_evap;
        }
        else {
            if (energy->Tgrnd >= CONST_TKFRZ) {
                soilfrost = fabs(vapor_grnd); // mm/s
            }
            else {
                soil_dew = fabs(vapor_grnd); // mm/s
            }
        }
        double soil_water = liq[0] * dz_soil[0] * CONST_RHOFW + ice[0] * dz_soil[0] * CONST_RHOICE;
        soil_limit = soil_water / ((1.0 - coverage) * step_dt);
        if (vapor_grnd > soil_limit) {
            soil_evap = min(soil_evap, liq[0] * dz_soil[0] * CONST_RHOFW / ((1.0 - coverage) * step_dt));
            soil_sublim = min(soil_sublim, ice[0] * dz_soil[0] * CONST_RHOICE / ((1.0 - coverage) * step_dt));
            energy->LatentSoil = soil_limit * LatentVapGrnd;
            energy->SensibleSoil += (vapor_grnd - soil_limit) * LatentVapGrnd;
        }
        energy->latent = energy->LatentSnow * coverage + energy->LatentSoil * (1.0 - coverage);
        energy->sensible = energy->SensibleSnow * coverage + energy->SensibleSoil * (1.0 - coverage);
    }
    else {
        vapor_grnd = energy->latent / LatentVapGrnd; // mm/s
        if (vapor_grnd >= 0.0) {
            snow_sublim = 0.0;
            snow_evap = 0.0;
            if (snow->swq > 0.0 && coverage > 0.0) {
                snow_sublim = vapor_grnd;
                snow_limit = snow->swq / (coverage * step_dt);
                if (snow_sublim > snow_limit) {
                    snow_sublim = snow_limit;
                }
            }
            // bare ground flux
            soil_evap = 0.0;
            soil_sublim = 0.0;
            if (1.0 - coverage > 0.0) {
                if (liq[0] + ice[0] > 0.0) {
                    soil_evap = max(vapor_grnd * (liq[0] / (liq[0] + ice[0])), 0.0);
                }
                soil_sublim = vapor_grnd - soil_evap;
                double soil_water = liq[0] * dz_soil[0] * CONST_RHOFW +
                                    ice[0] * dz_soil[0] * CONST_RHOICE;
                soil_limit = soil_water / ((1.0 - coverage) * step_dt);
                if (vapor_grnd > soil_limit) {
                    soil_evap = min(soil_evap, liq[0] * dz_soil[0] * CONST_RHOFW /
                                    ((1.0 - coverage) * step_dt));
                    soil_sublim = min(soil_sublim, ice[0] * dz_soil[0] * CONST_RHOICE /
                                    ((1.0 - coverage) * step_dt));
                }
            }
        }
        else {
            if (coverage > 0.0) {
                if (energy->Tgrnd >= CONST_TKFRZ) {
                    snowfrost = fabs(vapor_grnd);
                }
                else {
                    snow_dew = fabs(vapor_grnd);
                }
            }
            if (1.0 - coverage > 0.0) {
                if (energy->Tgrnd >= CONST_TKFRZ) {
                    soilfrost = fabs(vapor_grnd);
                }
                else {
                    soil_dew = fabs(vapor_grnd);
                }
            }
        }
        /* actual cell-average latent heat flux after snow/soil water limitation */
        double vapor_actual = coverage * (snow_evap + snow_sublim) +
                              (1.0 - coverage) * (soil_evap + soil_sublim);
        /* conserve total energy flux. The difference between evaporation demand
           and actual evaporation is returned to sensible heat. */
        if (vapor_actual < vapor_grnd) {
            energy->latent = vapor_actual * LatentVapGrnd;
            energy->sensible += (vapor_grnd - vapor_actual) * LatentVapGrnd;
        }
    }

    // 计算雪层冻结和融化量
    for (i = 0; i < Nsnow; i++) {
        double delta_ice = theta_ice[i] - last_thice[i];
        if (delta_ice > 0.0) {
            pack_frze[i] = delta_ice * dz_snow[i] * CONST_RHOICE * coverage;
            pack_melt[i] = 0.0;
        }
        else if (delta_ice < 0.0) {
            pack_frze[i] = 0.0;
            pack_melt[i] = -delta_ice * dz_snow[i] * CONST_RHOICE * coverage;
        }
        else {
            pack_frze[i] = 0.0;
            pack_melt[i] = 0.0;
        }
    }
    /*******************************
      Snowpack hydrology processes
    *******************************/
    if (Nsnow == 0 && snow->swq > 0.0) {
        double old_swq = snow->swq;
        snow->swq += (snowfrost - snow_sublim) * step_dt * coverage;
        snow->swq = max(snow->swq, 0.0);
        double ratio = snow->swq / old_swq;
        snow->snow_depth = max(0.0, ratio * snow->snow_depth);
        snow->snow_depth = min(max(snow->snow_depth, snow->swq / 500.0), snow->swq / 50.0);
        if (snow->swq < 0.0) {
            ice[0] += snow->swq / (dz_soil[0] * MM_PER_M);
            snow->swq = 0.0;
            snow->snow_depth = 0.0;
        }
        // soil layer evaporation/deposition
        double snow_out = (snow_dew - snow_evap) * step_dt * coverage;
        double soil_inflow = ((soil_dew - soil_evap) * 
                             (1.0 - coverage) + rainfall) * step_dt + snow_out;
        if (soil_inflow < 0.0) {
            liq[0] += soil_inflow / (dz_soil[0] * MM_PER_M);
            cell->soil_inflow = 0.0;
        }
        else {
            cell->soil_inflow = soil_inflow;
        }
        ice[0] += (soilfrost - soil_sublim) * step_dt * (1.0 - coverage) / (dz_soil[0] * MM_PER_M);
        if (ice[0] < 0.0) {
            liq[0] += ice[0];
            ice[0] = 0.0;
        }
        cell->moist[0] = liq[0] + ice[0] * CONST_RHOICE / CONST_RHOFW;
    }
    if (snow->snow_depth < param.TOL_B || snow->swq < param.TOL_A) {
        snow->snow_depth = 0.0;
        snow->swq = 0.0;
        snow->coverage = 0.0;
    }

    /* for multi-layer (>= 1) snow */
    if (Nsnow > 0) {
        // top layer total snow water before sublimation
        pack_ice[0] += (snowfrost - snow_sublim) * step_dt * coverage;
        pack_ice[0] = max(pack_ice[0], 0.0);
        pack_liq[0] += (snow_dew - snow_evap) * step_dt * coverage;
        pack_liq[0] = max(0.0, pack_liq[0]);
        /* soil layer evaporation/deposition */
        double soil_inflow = (soil_dew - soil_evap + rainfall) * step_dt * (1.0 - coverage);
        if (soil_inflow < 0.0) {
            liq[0] += soil_inflow / (dz_soil[0] * MM_PER_M);
            cell->soil_inflow = 0.0;
        }
        else {
            cell->soil_inflow = soil_inflow;
        }
        ice[0] += (soilfrost - soil_sublim) * step_dt * (1.0 - coverage) / (dz_soil[0] * MM_PER_M);
        if (ice[0] < 0.0) {
            liq[0] += ice[0];
            ice[0] = 0.0;
        }
        cell->moist[0] = liq[0] + ice[0] * CONST_RHOICE / CONST_RHOFW;
    }

    /* compute inter-layer snow water flow */
    if (Nsnow > 0) {
        double snow_inflow = 0.0;
        for (i = 0; i < Nsnow; i++) {
            snow_density = pack_ice[i] / (dz_snow[i] * coverage);
            if (snow_density >= 200.0) {
                liquid_capacity = param.SNOW_LIQUID_WATER_CAPACITY;
            }
            else {
                liquid_capacity = param.SNOW_LIQUID_WATER_CAPACITY + 0.07 * (200.0 - snow_density) / 200.0;
            }
            // 更新雪层水分含量和冰分数
            pack_liq[i] += snow_inflow;
            theta_ice[i] = min(1.0, pack_ice[i] / (dz_snow[i] * coverage * CONST_RHOICE));
            porosity[i] = 1.0 - theta_ice[i];
            theta_liq[i] = min(porosity[i], pack_liq[i] / (dz_snow[i] * coverage * CONST_RHOFW));
            /* excess liquid water snow-area basis */
            double outflow = max(0.0, (theta_liq[i] - liquid_capacity * 
                                porosity[i]) * dz_snow[i] * CONST_RHOFW);
            pack_outflow[i] = outflow * coverage;
            if (i == 0) {
                pack_outflow[i] = max((theta_liq[i] - porosity[i]) * dz_snow[i] * CONST_RHOFW * 
                                    coverage, param.SNOW_RELEASE_FAC * step_dt * pack_outflow[i]);
            }
            pack_liq[i] -= pack_outflow[i];
            if (pack_liq[i] / (pack_ice[i] + pack_liq[i]) > param.SNOW_MAX_LIQUID_FRAC) {
                double excess_liq = param.SNOW_MAX_LIQUID_FRAC / (1.0 - param.SNOW_MAX_LIQUID_FRAC) * pack_ice[i];
                pack_outflow[i] += pack_liq[i] - excess_liq;
                pack_liq[i] = excess_liq;
            }
            snow_inflow = pack_outflow[i];
            theta_liq[i] = pack_liq[i] / (dz_snow[i] * coverage * CONST_RHOFW);
        }
    }
    /* Liquid water from snow bottom to soil [mm/s] */
    if (Nsnow > 0) {
        lidx = Nsnow - 1;
        snow->snow_outflow = pack_outflow[lidx] / step_dt;
        for (i = 0; i < Nsnow; i++) {
            pack_outflow[i] /= step_dt;
        }
    }
    else {
        snow->snow_outflow = 0.0;
    }
    // obtain equilibrium state of snow in glacier region
    if (snow->swq > param.SNOW_MAX_SURFACE_SWE && Nsnow > 0) {
        snow_density = pack_ice[0] / (dz_snow[0] * coverage);
        excess_flux = snow->swq - param.SNOW_MAX_SURFACE_SWE;
        pack_ice[0] = pack_ice[0] - excess_flux;
        dz_snow[0] -= excess_flux / snow_density;
        excess_flux /= step_dt;
    }

    // snow layer compaction, combination, and division
    if (Nsnow > 0) {
        /* snow layer combination */
        snow_compaction(step_dt, 
                        air_temp, wind,
                        pressure, snow);
        /* snow layer combination */
        snow_combination(dz_soil[0], cell, snow);
        /* snow layer division */
        snow_division(snow);
    }

    // 计算累计厚度和中心位置
    double cum_depth = 0.0;
    for (i = 0; i < snow->Nsnow; i++) {
        cum_depth += snow->dz_snow[i];
        snow->Zsum_snow[i] = cum_depth;  // 累计厚度到该层底部
        snow->zc_snow[i] = cum_depth - snow->dz_snow[i] / 2.0;  // 中心位置
    }

    /* Update SnowDepth and snow mass for multi-layer snow */
    if (snow->Nsnow > 0) {
        snow->swq = 0.0;
        snow->snow_depth = 0.0;
        for (i = 0; i < snow->Nsnow; i++) {
            snow->swq += pack_ice[i] + pack_liq[i];
            snow->snow_depth += dz_snow[i];
        }
    }

    // update snow quantity
    if (snow->snow_depth < param.TOL_A || snow->swq < param.TOL_A) {
        snow->swq = 0.0;
        snow->snow_depth = 0.0;
    }

    /* accumulate glacier excessive flow [mm] */
    snow->glac_excess = excess_flux * step_dt;

    // Update node properties
    update_nodes(pressure, 
                 energy,  cell, 
                 snow, soil_con);

    // === 损失项（正值） ===
    cell->soil_evap = soil_evap * (1.0 - coverage);       // 裸土蒸发
    cell->soil_sublim = soil_sublim * (1.0 - coverage);   // 裸土升华
    cell->snow_evap = snow_evap * coverage;               // 雪面蒸发
    cell->snow_sublim = snow_sublim * coverage;           // 雪面升华

    // === 收入项（负值） ===
    cell->soil_dew = soil_dew * (1.0 - coverage);         // 裸土露水
    cell->soil_frost = soilfrost * (1.0 - coverage);      // 裸土霜凝华
    cell->snow_dew = snow_dew * coverage;                 // 雪面露水
    cell->snow_frost = snowfrost * coverage;              // 雪面霜凝华

    cell->evap = (cell->soil_evap + cell->soil_sublim + cell->snow_evap + cell->snow_sublim)
               - (cell->soil_dew + cell->soil_frost + cell->snow_dew + cell->snow_frost);

    return (0);
}