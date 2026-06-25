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
               double             snowfall,
               double             rainfall,
               double             pressure,
               energy_bal_struct *energy,
               cell_data_struct  *cell,
               snow_data_struct  *snow,        
               soil_con_struct   *soil_con)
{
    extern parameters_struct param;
    size_t      i;  
    double      vapor_grnd;
    double      conden_grnd;
    double      snow_sublim;
    double      snowfrost;
    // 指针赋值               
    double *liq = cell->liq;
    double *ice = cell->ice;
    double *dz_snow = snow->dz_snow;
    double *porosity = snow->porosity;
    double *theta_liq = snow->theta_liq;
    double *theta_ice = snow->theta_ice;
    double *dz_soil = soil_con->dz_soil;
    double *pack_ice = snow->pack_ice;
    double *pack_liq = snow->pack_liq;
    double *pack_outflow = snow->pack_outflow;

    /* initialize */
    double snow_density = 0.0;
    double snow_inflow = 0.0;
    double excess_flux = 0.0;
    double latent = energy->latent;
    double liquid_capacity = 0.0;
    double LatentVapGrnd = energy->LatentVapGrnd;
    /** compute soil/snow surface evap,
        dew rate based on energy flux. **/
    // positive part of ground latent heat
    vapor_grnd = max(latent / LatentVapGrnd, 0.);
    // negative part of ground latent heat
    conden_grnd = fabs(min(latent / LatentVapGrnd, 0.));
    snow->last_swq = snow->swq;
    for (i = 0; i < snow->Nsnow; i++) {
        snow->last_snowfrac[i] = snow->snow_frac[i];
    }
    cell->evap = vapor_grnd - conden_grnd;

    /* ground sublimation and evaporation */
    snow_sublim = 0.0;
    if (snow->swq > 0.0) {
        snow_sublim = min(vapor_grnd, snow->swq / step_dt);
    }

    double esoil = vapor_grnd - snow_sublim;
    /* ground frost and dew */
    snowfrost = 0.0;
    if (snow->swq > 0.0) {
        snowfrost = conden_grnd;
    }
    double dewsoil = conden_grnd - snowfrost;

    /* snowpack water processs */
    update_snow(air_temp, step_dt,
                snowfall, snow);
    
    if (snow->Nsnow > 0) {
        /* snow layer combination */
        snow_compaction(step_dt, 
                        air_temp,
                        pressure, snow);
        /* snow layer combination */
        snow_combination(dz_soil[0], cell, snow);
        /* snow layer division */
        snow_division(snow);
    }

    /*******************************
      Snowpack hydrology processes
    *******************************/
    size_t Nsnow = snow->Nsnow;
    if (snow->swq == 0.0) {
        ice[0] += (snowfrost - snow_sublim) * step_dt / (dz_soil[0] * MM_PER_M);
        if (ice[0] < 0.0) {
            liq[0] += ice[0];
            ice[0] = 0.0;
        }
    }
    if (Nsnow == 0 && snow->swq > 0.0) {
        double tmp_swq = snow->swq;
        snow->swq = snow->swq - snow_sublim * step_dt + snowfrost * step_dt;   
        double ratio = snow->swq / tmp_swq;
        snow->snow_depth = max(0.0, ratio * snow->snow_depth);
        snow->snow_depth = min(max(snow->snow_depth, snow->swq / 500.0), snow->swq / 50.0);
        if (snow->swq < 0.0) {
            ice[0] += snow->swq / (dz_soil[0] * MM_PER_M);
            snow->swq = 0.0;
            snow->snow_depth = 0.0;
        }
        if (ice[0] < 0.0) {
            liq[0] += ice[0];
            ice[0] = 0.0;
        }
    }
    if (snow->snow_depth < param.TOL_B || snow->swq < param.TOL_A) {
        snow->snow_depth = 0.0;
        snow->swq = 0.0;
    }

    /* for multi-layer (>= 1) snow */
    if (Nsnow > 0) {
        double tmp_liq = pack_ice[0] + pack_liq[0]; // top layer total snow water before sublimation
        double tmp_ice = pack_ice[0] - snow_sublim * step_dt + snowfrost * step_dt;
        pack_ice[0] = tmp_ice;
        if (tmp_ice < param.TOL_A && Nsnow > 0) {
            snow_combination(dz_soil[0], cell, snow);
        }
        if (tmp_ice > param.TOL_A && Nsnow > 0) {
            dz_snow[0] *= (pack_ice[0] + pack_liq[0]) / tmp_liq;
        }
        if (Nsnow > 0) {
            pack_liq[0] += rainfall * step_dt * snow->coverage;
            pack_liq[0] = max(0.0, pack_liq[0]);
        }
    }

    /* Porosity and partial volume */
    for (i = 0; i < Nsnow; i++) {
        theta_ice[i] = min(1.0, pack_ice[i] / (dz_snow[i] * CONST_RHOICE));
        porosity[i] = 1.0 - theta_ice[i];
    }
    /* compute inter-layer snow water flow */
    if (Nsnow > 0) {
        for (i = 0; i < Nsnow; i++) {
            snow_density = pack_ice[i] / dz_snow[i];
            if (snow_density >= 200.0) {
                liquid_capacity = param.SNOW_LIQUID_WATER_CAPACITY;
            }
            else {
                liquid_capacity = param.SNOW_LIQUID_WATER_CAPACITY + 0.07 * (200.0 - snow_density) / 200.0;
            }
            pack_liq[i] += snow_inflow;
            theta_liq[i] = pack_liq[i] / (dz_snow[i] * CONST_RHOFW);
            pack_outflow[i] = max(0.0, (theta_liq[i] - liquid_capacity *
                                                    porosity[i]) * dz_snow[i]);

            if (i == 0) {
                pack_outflow[i] = max((theta_liq[i] - porosity[i]) * dz_snow[i], 
                                                    param.SNOW_RELEASE_FAC * step_dt * pack_outflow[i]);
            }
            pack_outflow[i] *= CONST_RHOFW;
            pack_liq[i] -= pack_outflow[i];
            if (pack_liq[i] / (pack_ice[i] + pack_liq[i]) > param.SNOW_MAX_LIQUID_FRAC) {
                pack_outflow[i] += pack_liq[i] - param.SNOW_MAX_LIQUID_FRAC / 
                                                    (1.0 - param.SNOW_MAX_LIQUID_FRAC) * pack_ice[i];
                pack_liq[i] = param.SNOW_MAX_LIQUID_FRAC / 
                                                    (1.0 - param.SNOW_MAX_LIQUID_FRAC) * pack_ice[i];
            }
            snow_inflow = pack_outflow[i];
            theta_liq[i] = pack_liq[i] / (dz_snow[i] * CONST_RHOFW);
        }
    }
    else {
        for (i = 0; i < MAX_SNOWS; i++) {
            pack_outflow[i] = 0.0;  // no snow, no outflow
        }
    }
    for (i = 0; i < Nsnow; i++) {
        dz_snow[i] = max(dz_snow[i], pack_liq[i] / CONST_RHOFW + pack_ice[i] / CONST_RHOICE);
    }
    /* Liquid water from snow bottom to soil [mm/s] */
    if (Nsnow > 0) {
        snow->snow_outflow = pack_outflow[Nsnow-1] / step_dt;
        for (i = 0; i < Nsnow; i++) {
            pack_outflow[i] /= step_dt;
        }
    }
    else {
        snow->snow_outflow = 0.0;
    }
    // obtain equilibrium state of snow in glacier region
    if (snow->swq > param.SNOW_MAX_SURFACE_SWE && Nsnow > 0) {
        snow_density = pack_ice[0] / dz_snow[0];
        excess_flux = snow->swq - param.SNOW_MAX_SURFACE_SWE;
        pack_ice[0] = pack_ice[0] - excess_flux;
        dz_snow[0] -= excess_flux / snow_density;
        excess_flux /= step_dt;
    }
    // 计算累计厚度和中心位置
    double cum_depth = 0.0;
    for (i = 0; i < Nsnow; i++) {
        cum_depth += snow->dz_snow[i];
        snow->Zsum_snow[i] = -cum_depth;  // 累计厚度到该层底部
        snow->zc_snow[i] = -(cum_depth - snow->dz_snow[i]/2.0);  // 中心位置（负值）
    }

    /* sum up snow mass for layered snow */
    if (Nsnow > 0) {
        snow->swq = 0.0;
        for (i = 0; i < Nsnow; i++) {
            snow->swq += pack_ice[i] + pack_liq[i];
        }
    }

    /* Update SnowDepth for multi-layer snow */
    if (Nsnow > 0) {
        snow->snow_depth = 0.0;
        for (i = 0; i < Nsnow; i++) {
            snow->snow_depth += dz_snow[i];
        }
    }
    // update snow quantity
    if (snow->snow_depth < param.TOL_A || snow->swq < param.TOL_A) {
        snow->swq = 0.0;
        snow->snow_depth = 0.0;
    }
    // update snow fraction for multi-layer snow
    if (Nsnow > 0) {
        for (i = 0; i < Nsnow; i++) {
            snow->snow_frac[i] = pack_ice[i] / (pack_ice[i] + pack_liq[i]);
        }
    }
    else {
        for (i = 0; i < MAX_SNOWS; i++) {
            snow->snow_frac[i] = 0.0;
        }
    }
    /* accumulate glacier excessive flow [mm] */
    snow->glac_excess += excess_flux * step_dt;

    /* frozen ground and soil */
    if (energy->FrozenGrnd) {
        ice[0] += (dewsoil - esoil) * step_dt / (dz_soil[0] * MM_PER_M);
        dewsoil = 0.0;
        esoil = 0.0;
        if (ice[0] < 0.0) {
            liq[0] += ice[0];
            ice[0] = 0.0;
        }
        cell->moist[0] = liq[0] + ice[0];
    }

    /* set esoil and transpiration mm/s -> m/s */
    cell->esoil = esoil / MM_PER_M;
    cell->dewsoil = dewsoil;
    cell->snowfrost = snowfrost;
    cell->snow_sublim = snow_sublim;

    return (0);
}