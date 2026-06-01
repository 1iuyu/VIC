/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate glacier accumulation and melt
 *****************************************************************************/

#include "vic_run.h"

/*****************************************************************************
  * @brief   Calculate glacier accumulation and melt using an energy balance
  *          approach for a two layer model
 *****************************************************************************/
int
calc_water_bal_glac(double             step_dt,
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
    double      tmp_ice[MAX_SOILS];
    double      tmp_liq[MAX_SOILS];          
    double      snow_density;

    // 指针赋值            
    double *liq = cell->liq;
    double *ice = cell->ice;
    double *moist = cell->moist;
    double *dz_snow = snow->dz_snow;
    double *porosity = snow->porosity;
    double *pack_ice = snow->pack_ice;
    double *pack_liq = snow->pack_liq;
    double *theta_liq = snow->theta_liq;
    double *theta_ice = snow->theta_ice;
    double *dz_soil = soil_con->dz_soil;
    double latent = energy->latent;
    double LatentVapGrnd = energy->LatentVapGrnd;
    double *pack_outflow = snow->pack_outflow;  

    /* initialize */
    for (i = 0; i < cell->Nsoil; i++) {
        tmp_ice[i] = 0.;
        tmp_liq[i] = 0.;
    }

    /* prepare for water process */
    for (i = 0; i < cell->Nsoil; i++) {
        ice[i] = moist[i] - liq[i];
        if (ice[i] < 0.) {
            ice[i] = 0.;
        }
        tmp_liq[i] = liq[i];
        tmp_ice[i] = ice[i];
    }
    // record snow water content and snow fraction
    snow->last_swq = snow->swq;
    for (i = 0; i < snow->Nsnow; i++) {
        snow->last_snowfrac[i] = snow->snow_frac[i];
    }
    double snow_inflow = 0.;
    double soil_inflow = 0.;
    double liquid_capacity = 0.0;
    /* compute soil/snow surface evap/dew rate based on energy flux */
    double vapor_grnd = max(latent / LatentVapGrnd, 0.); // positive part of ground latent heat
    double conden_grnd = fabs(min(latent / LatentVapGrnd, 0.)); // negative part of ground latent heat
    cell->evap = vapor_grnd - conden_grnd;

    /* ground sublimation and evaporation */
    double snow_sublim = vapor_grnd;
    /* ground frost and dew */
    double snow_frost = conden_grnd;
    cell->snowfrost = snow_frost;
    cell->snow_sublim = snow_sublim;
    if (snowfall > 0.0) {
        snow->delta_depth = snowfall / snow->new_snow_density;
    }
    else {
        snow->delta_depth = 0.0;
    }
    
    /* snowpack water processs */
    update_snow(air_temp, step_dt,
                snowfall, snow);

    if (snow->Nsnow > 0) {
        /* snow layer combination */
        snow_compaction(step_dt,
                        air_temp,
                        pressure, energy,
                        snow);
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
        ice[0] += (snow_frost - snow_sublim) * step_dt / (dz_soil[0] * MM_PER_M);
        if (ice[0] < 0.0) {
            liq[0] += ice[0];
            ice[0] = 0.0;
        }
    }
    if (Nsnow == 0 && snow->swq > 0.0) {
        double tmp_swq = snow->swq;
        snow->swq -= snow_sublim * step_dt + snow_frost * step_dt;
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
        double tmp_packliq = pack_ice[0] + pack_liq[0];
        double tmp_packice = pack_ice[0] - snow_sublim *
                                step_dt + snow_frost * step_dt;
        pack_ice[0] = tmp_packice;
        if (tmp_packice < param.TOL_A && Nsnow > 0) {
            snow_combination(dz_soil[0], cell, snow);
        }
        if (tmp_packice > param.TOL_A && Nsnow > 0) {
            dz_snow[0] *= (pack_ice[0] + pack_liq[0]) / tmp_packliq;
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
                pack_outflow[i] = max((theta_liq[i] - porosity[i]) * 
                                        dz_snow[i], param.SNOW_RELEASE_FAC * 
                                                step_dt * pack_outflow[i]);
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
            pack_outflow[i] = 0.0;
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
    for (i = 0; i < Nsnow; i++) {
        pack_outflow[i] /= step_dt;
    }

    // to obtain equilibrium state of snow in glacier
    double excess_flux = 0.0;
    if (snow->swq > param.SNOW_MAX_SURFACE_SWE && Nsnow > 0) {
        snow_density = pack_ice[0] / dz_snow[0];
        excess_flux = snow->swq - param.SNOW_MAX_SURFACE_SWE;
        pack_ice[0] = pack_ice[0] - excess_flux;
        dz_snow[0] -= excess_flux / snow_density;
        excess_flux /= step_dt / MM_PER_M; // [m/s]
    }
    // 计算累计厚度和中心位置
    double cum_depth = 0.0;
    for (i = 0; i < Nsnow; i++) {
        cum_depth += snow->dz_snow[i];
        snow->Zsum_snow[i] = -cum_depth;  // 累计厚度到该层底部
        snow->zc_snow[i] = -(cum_depth - snow->dz_snow[i]/2.0);  // 中心位置（负值）
    }
    // 更新节点数据
    update_node(cell, snow, soil_con);

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

    // 限制雪层变量
    if (snow->snow_depth < param.TOL_A || snow->swq < param.TOL_A) {
        snow->swq = 0.0;
        snow->snow_depth = 0.0;
    }

    /* accumulate glacier excessive flow [m] */
    snow->glac_excess += excess_flux * step_dt;

    /* 计算进入冰层表面的水量 */
    // pack_melt: 雪层融化的水量[mm]
    // pack_comb: 薄雪层合并中液体所导致的表面积水[mm]
    // pack_transp: 从多层积雪转变为无积雪时由薄雪造成的地表积水深度[mm]
    soil_inflow = (snow->pack_melt + snow->pack_comb + 
                        snow->pack_transp) / step_dt / MM_PER_M;
    // 添加雪层多余水分,露水,降水
    soil_inflow += (snow->snow_outflow + rainfall * 
                    (1.0 - snow->coverage)) / MM_PER_M; // convert to [m/s]

    cell->soil_inflow = soil_inflow;
    cell->runoff = soil_inflow; // convert to [m/s]
    double liq_repl_sub = 0.0;
    for (i = 0; i < cell->Nsoil; i++) {
        liq_repl_sub += dz_soil[i] * 
                (ice[i] - tmp_ice[i] + liq[i] - tmp_liq[i]);
    }
    liq_repl_sub /= step_dt; // convert to [m/s]
    for (i = 0; i < cell->Nsoil; i++) {
        ice[i] = min(1.0, tmp_ice[i]);
        liq[i] = 1.0 - ice[i];
    }
    cell->baseflow = excess_flux + liq_repl_sub;
    
    if (cell->baseflow < 0.0) {
        cell->baseflow = 0.0;
    }

    return (0);

}
