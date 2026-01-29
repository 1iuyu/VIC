/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate glacier accumulation and melt
 *****************************************************************************/

#include <vic_run.h>

/*****************************************************************************
  * @brief   Calculate glacier accumulation and melt using an energy balance
  *          approach for a two layer model
 *****************************************************************************/
int
calc_water_bal(size_t             hidx,
               double             step_dt,
               double             air_temp,
               double             snowfall,
               double             rainfall,
               force_data_struct *force,
               energy_bal_struct *energy,
               cell_data_struct  *cell,
               snow_data_struct  *snow,
               veg_var_struct    *veg_var,               
               soil_con_struct   *soil_con)
{
    extern option_struct     options;
    extern parameters_struct param;

    size_t      lindex;
    size_t      i;
    int         ErrorFlag;
    double      LatentVapGrnd;      
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
    double *soil_transp = cell->soil_transp;
    double *transp_fact = cell->transp_fact;
    double *pack_outflow = snow->pack_outflow;

    /* initialize */
    double ppt = 0.;
    double snow_density = 0.0;
    double snow_inflow = 0.0;
    double soil_inflow = 0.0;
    double snow_outflow = 0.0;
    double snow_depth = 0.0;
    double excess_flux = 0.0;
    double latent = energy->latent;
    double pressure = force->pressure[hidx];
    LatentVapGrnd = energy->LatentVapGrnd;
    /** compute soil/snow surface evap,
        dew rate based on energy flux. **/
    // positive part of ground latent heat
    vapor_grnd = max(latent / LatentVapGrnd, 0.);
    // negative part of ground latent heat
    conden_grnd = fabs(min(latent / LatentVapGrnd, 0.)); 
    cell->evap = vapor_grnd - conden_grnd;
    cell->vapor_grnd = vapor_grnd;
    cell->conden_grnd = conden_grnd;
    /* Calculation of evaporation from the canopy */
    cell->canopyevap = canopy_evap(step_dt, energy,
                                   cell, snow, veg_var);

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
                        pressure,
                        cell, snow);
        /* snow layer combination */
        snow_combination(dz_soil[0], cell, snow);
        /* snow layer division */
        snow_division(snow);
    }

    /*******************************
      Snowpack hydrology processes
    *******************************/
    if (snow->swq == 0.0) {
        ice[0] += (snowfrost - snow_sublim) * step_dt / (dz_soil[0] * MM_PER_M);
        if (ice[0] < 0.0) {
            liq[0] += ice[0];
            ice[0] = 0.0;
        }
    }
    if (snow->Nsnow == 0 && snow->swq > 0.0) {
        double tmp_swq = snow->swq;
        snow->swq -= snow_sublim * step_dt + snowfrost * step_dt;      
        double ratio = snow->swq / tmp_swq;
        snow_depth = max(0.0, ratio * snow_depth);
        snow_depth = min(max(snow_depth, snow->swq / 500.0), snow->swq / 50.0);
        snow->snow_depth = snow_depth;
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
    if (snow->Nsnow > 0) {
        double tmp_liq = pack_ice[0] + pack_liq[0]; // top layer total snow water before sublimation
        double tmp_ice = pack_ice[0] - snow_sublim * step_dt + snowfrost * step_dt;
        pack_ice[0] = tmp_ice;
        if (tmp_ice < param.TOL_A && snow->Nsnow > 0) {
            snow_combination(dz_soil[0], cell, snow);
        }
        if (tmp_ice > param.TOL_A && snow->Nsnow > 0) {
            dz_snow[0] *= (pack_ice[0] + pack_liq[0]) / tmp_liq; 
        }
        if (snow->Nsnow > 0) {
            pack_liq[0] += rainfall * step_dt;
            pack_liq[0] = max(0.0, pack_liq[0]);
        }
    }

    /* Porosity and partial volume */
    for (i = 0; i < snow->Nsnow; i++) {
        theta_ice[i] = min(1.0, pack_ice[i] / (dz_snow[i] * CONST_RHOICE));
        porosity[i] = 1.0 - theta_ice[i];
    }
    /* compute inter-layer snow water flow */
    for (i = 0; i < snow->Nsnow; i++) {
        pack_liq[i] += snow_inflow;
        theta_liq[i] = pack_liq[i] / (dz_snow[i] * CONST_RHOFW);
        pack_outflow[i] = max(0.0, (theta_liq[i] - 
                            param.SNOW_LIQUID_WATER_CAPACITY * 
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
    for (i = 0; i < snow->Nsnow; i++) {
        dz_snow[i] = max(dz_snow[i], pack_liq[i] / CONST_RHOFW + pack_ice[i] / CONST_RHOICE);
    }
    /* Liquid water from snow bottom to soil [mm/s] */
    snow_outflow = pack_outflow[0] / step_dt;
    for (i = 0; i < snow->Nsnow; i++) {
        pack_outflow[i] /= step_dt;
    }
    // obtain equilibrium state of snow in glacier region
    if (snow->swq > param.SNOW_MAX_SURFACE_SWE) {
        snow_density = pack_ice[0] / dz_snow[0];
        excess_flux = snow->swq - param.SNOW_MAX_SURFACE_SWE;
        pack_ice[0] = pack_ice[0] - excess_flux;
        dz_snow[0] -= excess_flux / snow_density;
        excess_flux /= step_dt;
    }

    /* sum up snow mass for layered snow */
    if (snow->Nsnow > 0) {
        snow->swq = 0.0;
        for (i = 0; i < snow->Nsnow; i++) {
            snow->swq += pack_ice[i] + pack_liq[i];
        }
    }

    /* Update SnowDepth for multi-layer snow */
    if (snow->Nsnow > 0) {
        snow->snow_depth = 0.0;
        for (i = 0; i < snow->Nsnow; i++) {
            snow->snow_depth += dz_snow[i];
        }
    }
    // update snow quantity
    if (snow->snow_depth < param.TOL_A || snow->swq < param.TOL_A) {
        snow->swq = 0.0;
        snow->snow_depth = 0.0;
    }

    /* accumulate glacier excessive flow [mm] */
    snow->glac_excess += excess_flux * step_dt;

    /* frozen ground and soil */
    if (options.FROZEN_SOIL && energy->FrozenGrnd) {
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

    for (i = 0; i < options.Nroot; i++) {
        soil_transp[i] = cell->transp * transp_fact[i];
    }
    /* soil_inflow */
    soil_inflow = (snow->pack_melt + snow->pack_comb +
                            snow->pack_transp) / step_dt / MM_PER_M;
    if (snow->Nsnow == 0) {
        soil_inflow += (snow_outflow + dewsoil + rainfall) / MM_PER_M;
    }
    else {
        soil_inflow += (snow_outflow + dewsoil) / MM_PER_M;
    }
    
    ppt += soil_inflow;  // m/s
    cell->soil_inflow = soil_inflow;
    /********************************************************
       Compute Runoff, Baseflow, and Soil Moisture Transport
    ********************************************************/

    ErrorFlag = runoff(ppt, cell, soil_con);

    cell->baseflow += snow->glac_excess;

    return(ErrorFlag);
    
}
