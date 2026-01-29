/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine iterates to determine the temperature of the canopy, and solve
 * the resulting fluxes between the canopy and the atmosphere and the canopy
 * and the ground.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    Calculate the canopy energy balance.
 *****************************************************************************/
int
func_canopy_energy_bal(double             step_dt,
                       double             wind,
                       double             Qair,
                       double             vp,
                       double             longwave,
                       double             air_temp,
                       double             air_density,
                       double             pressure,
                       double             rh_grnd,
                       double             ref_height,
                       double            *roughness,
                       double            *displacement,
                       energy_bal_struct *energy,
                       cell_data_struct  *cell,
                       snow_data_struct  *snow,
                       soil_con_struct   *soil_con,
                       veg_var_struct    *veg_var,
                       veg_lib_struct    *veg_lib)
{
    extern parameters_struct param;

    /* General Model Parameters */
    int         ErrorFlag;
    size_t      lindex;
    double      layer_depth;
    double      SVP_liq;
    double      SVP_ice;
    double      esat_slope;
    double      liq_slope;
    double      ice_slope;
    double      esat_Tsub;
    double      esat_Tcanopy;
    double      Fmo_tmp;
    double      layer_T;
    double      coef_sensible;
    double      coef_latent;
    double      psi_h_grnd;
    double      Ra_leaf;
    double      Ra_canopy;
    double      SensibleSub;
    double      NetLongSub;
    double      LatentSub;
    double      GroundSub;

    /* Initialization */
    double NetLAI = veg_var->NetLAI;
    double NetSAI = veg_var->NetSAI;
    double VPcanopy = cell->VPcanopy;

    double Canopy_Upper = veg_lib->Canopy_Upper;
    double Tair = energy->Tair;
    double Tcanopy = energy->Tcanopy;
    double Tupper = energy->Tupper;
    double total_transp = cell->total_transp;
    double NetShortGrnd = energy->NetShortGrnd;
    double EmissLongSub = energy->EmissLongSub;
    double EmissLongGrnd = energy->EmissLongGrnd;
    double *T = energy->T;
    double *Ra_over = cell->Ra_over;
    double *Ra_sub = cell->Ra_sub;
    double *kappa_node = energy->kappa_node;

    /* saturation vapor pressure at Tupper */
    svp(Tupper, &SVP_liq, &SVP_ice);
    if (Tupper > CONST_TKFRZ) {
        esat_Tsub = SVP_liq;
    }
    else {
        esat_Tsub = SVP_ice;
    }
    /* wind speed at canopy height */
    double tmp_wind = wind * log((Canopy_Upper - displacement[1] +
                            roughness[1]) / roughness[1]) /
                                    log(param.REF_HEIGHT_WIND / roughness[1]);
    if (Canopy_Upper - displacement[1] < 0.0) {
        log_err("Canopy_Upper: < displacement");
    }

    /********************************
      Vegetated canopy energy flux
    ********************************/
    size_t iter = 0;
    bool last_iter = false;
    double zL_over = 1.0;
    double ustar = 0.06;
    double delta_over = param.TOL_A + 1.;
    double sensible_over = 0.0;
    double sensible_sub = 0.0;
    double psi_h_new = 1.;
    /* prepare for longwave rad */
    double coef_lw_atmos = -EmissLongSub * (1.0 + (1.0 - EmissLongSub) * 
                        (1.0 - EmissLongGrnd)) * longwave -
                            EmissLongSub * EmissLongGrnd * 
                                    CONST_BOLTZ * pow(Tupper, 4);
    double coef_lw_canopy = (2.0 - EmissLongSub * (1.0 - EmissLongGrnd)) *
                                       EmissLongSub * CONST_BOLTZ;

    /* begin stability iteration for Tcanopy */
    do {

        /* Compute canopy aero_resist */
        ErrorFlag = CalcAerodynamic(iter, Qair,
                                    air_temp, &zL_over,
                                    &ustar, Ra_over,
                                    &sensible_over, wind,
                                    air_density,
                                    displacement[1],
                                    ref_height,
                                    roughness[1]);
        if (ErrorFlag == ERROR) {
            return (ERROR);
        }

        /***************************
          Under-canopy aero_resist
        ***************************/
        double L_grnd = 0.0;
        double zL_grnd = 0.0;

        if (iter > 0) {
            Fmo_tmp = CONST_KARMAN * (CONST_G / Tair) * sensible_sub /
                                                (air_density * CONST_CPDAIR);
            if (fabs(Fmo_tmp) < param.TOL_A) {
                Fmo_tmp = param.TOL_A;
            }
            L_grnd = -1.0 * pow(ustar, 3) / Fmo_tmp;
            zL_grnd = min((displacement[1] - roughness[0]) / L_grnd, 1.0);
        }
        if (zL_grnd < 0.0) {
            psi_h_new = pow((1.0 - 15.0 * zL_grnd), -0.25);
        }
        else {
            psi_h_new = 1.0 + 4.7 * zL_grnd;
        }
        if (iter == 0) {
            psi_h_grnd = psi_h_new;
        }
        else {
            psi_h_grnd = 0.5 * (psi_h_grnd + psi_h_new);
        }

        /* wind attenuation within canopy */
        double alpha_wind = pow((veg_lib->alpha_canopy * (NetLAI + NetSAI) * 
                                            Canopy_Upper * psi_h_grnd), 0.5);
        double tmp_exp1 = exp(-alpha_wind * roughness[0] / Canopy_Upper);
        double tmp_exp2 = exp(-alpha_wind * (roughness[1] + displacement[1]) / Canopy_Upper);
        double FRA_h_tmp = Canopy_Upper * exp(alpha_wind) / alpha_wind * (tmp_exp1 - tmp_exp2);
        double Ch_turb = max(CONST_KARMAN * ustar * (Canopy_Upper - displacement[1]), param.TOL_A);
        Ra_sub[0] = 0.0;
        Ra_sub[1] = FRA_h_tmp / Ch_turb;
        Ra_sub[2] = Ra_sub[1];

        /* leaf boundary layer resistance */
        double tmp_side = alpha_wind * 50.0 / (1.0 - exp(-alpha_wind / 2.0));
        Ra_leaf = tmp_side * sqrt(veg_lib->d_leaf / tmp_wind);
        Ra_leaf = min(max(Ra_leaf, 5.0), 50.0);
        cell->Ra_leaf = Ra_leaf;

        /* Saturated vapor pressure at Tcanopy */
        svp(Tcanopy, &SVP_liq, &SVP_ice);
        svp_slope(Tcanopy, &liq_slope, &ice_slope);
        if (Tcanopy > CONST_TKFRZ) {
            esat_Tcanopy = SVP_liq;
            esat_slope = liq_slope;
        }
        else {
            esat_Tcanopy = SVP_ice;
            esat_slope = ice_slope;
        }

        /**********************
          Stomatal Resistance
        **********************/
        if (iter == 0) {
            // sunlit case
            calc_rc(0, VPcanopy,
                    pressure,
                    total_transp,
                    energy, 
                    veg_var, veg_lib);
            // shaded case
            calc_rc(1, VPcanopy,
                    pressure,
                    total_transp,
                    energy, 
                    veg_var, veg_lib);
        }
        /* Iteratively solve the Tcanopy */
        Tcanopy = calc_atmos_energy_bal(step_dt, air_temp,
                                        air_density,
                                        esat_Tsub, vp,
                                        coef_lw_atmos,
                                        coef_lw_canopy,
                                        veg_lib->c_biomass,                                   
                                        esat_Tcanopy,
                                        esat_slope, &delta_over,
                                        energy, cell, veg_var);
        // update Tair and VPcanopy.
        Tair = energy->Tair;
        VPcanopy = cell->VPcanopy;
        energy->Tcanopy = Tcanopy;
        sensible_over = air_density * CONST_CPDAIR * 
                            (Tair - air_temp) /
                                        Ra_over[1];
        sensible_sub = air_density * CONST_CPDAIR *
                            (Tupper - Tair) /
                                        Ra_sub[1];
        
        double Qair_surf = (0.622 * VPcanopy) / (pressure - 0.378 * VPcanopy);

        if (last_iter) {
            break;
        }
        if (iter >= 5 && fabs(delta_over) < 0.01 && !last_iter) {
            last_iter = true;
        }
        iter++;
    }
    while (iter < param.MAX_ITER_OVER && fabs(delta_over) > 0.01);
    
    /********************************
      Vegetated surface energy flux
    ********************************/
    coef_lw_atmos = -EmissLongGrnd * (1.0 - EmissLongSub) * 
                        longwave - EmissLongGrnd * EmissLongSub * 
                                     CONST_BOLTZ * pow(Tcanopy, 4);
    coef_lw_canopy = EmissLongGrnd * CONST_BOLTZ;
    coef_sensible = CONST_CPMAIR * air_density / Ra_sub[1];
    coef_latent = air_density * CONST_CPDAIR / (energy->PsyCh_grnd * 
                                        (Ra_sub[2] + cell->Ra_evap));

    iter = 0;
    if (snow->Nsnow > 0) {
        lindex = snow->Nsnow - 1; // 存在雪层
        layer_depth = snow->dz_snow[lindex];
        layer_T = snow->pack_T[0];
    }
    else {
        layer_depth = soil_con->dz_soil[0]; // 不存在雪层，kappa_node第一层为土壤热力学参数
        layer_T = cell->soil_T[0];
    }
    
    double coef_ground = 2.0 * kappa_node[0] / layer_depth;
    double delta_sub = param.TOL_A + 1.;
    /* begin stability iteration for Tupper */
    do {

        /* Saturated vapor pressure at Tgrnd */
        svp(Tupper, &SVP_liq, &SVP_ice);
        svp_slope(Tupper, &liq_slope, &ice_slope);
        if (Tupper > CONST_TKFRZ) {
            esat_Tsub = SVP_liq;
            esat_slope = liq_slope;
        }
        else {
            esat_Tsub = SVP_ice;
            esat_slope = ice_slope;
        }
        /* Calculate the sensible heat flux */
        SensibleSub = calc_sensible_heat(Tupper, Tair,
                                         coef_sensible);
        /* Calculate longwave exchange and net radiation */
        NetLongSub = calc_longwave(Tupper, coef_lw_canopy) +
                                            coef_lw_atmos;
        /* Calculate Ground Heat Flux */
        GroundSub = calc_ground_heat(Tupper, layer_T,
                                     coef_ground);
        /* Computes the latent heat flux */
        LatentSub = calc_latent_heat(esat_Tsub,
                                     rh_grnd,
                                     VPcanopy,
                                     coef_latent);
        /* Calculate energy balance error at the snowpack surface */
        double RestTerm = NetShortGrnd - NetLongSub - SensibleSub -
                         LatentSub - GroundSub + energy->AdvectSub;
        double coef_flux = 4.0 * coef_lw_canopy * pow(Tupper, 3) + coef_sensible + 
                                         coef_latent * esat_slope + coef_ground;
        delta_sub = RestTerm / coef_flux;
        NetLongSub += 4.0 * coef_lw_canopy * pow(Tupper, 3) * delta_sub;
        SensibleSub += coef_sensible * delta_sub;
        LatentSub += coef_latent * esat_slope * delta_sub;
        GroundSub += coef_ground * delta_sub;

        /* update Tgrnd */
        Tupper += delta_sub;
        iter++;
    }
    while (iter < param.MAX_ITER_MOST && fabs(delta_sub) > 0.01);

    /* if snow on ground and Tupper > freezing point: reset 
        Tupper = freezing point. reevaluate ground fluxes. */
    if (snow->snow_depth > 0.05 && Tupper > CONST_TKFRZ) {

        Tupper = CONST_TKFRZ;
        /* update ground fluxes */
        NetLongSub = coef_lw_canopy * pow(Tupper, 4) - 
                            EmissLongGrnd * (1.0 - EmissLongSub) *
                            longwave - EmissLongGrnd * EmissLongSub * 
                            CONST_BOLTZ * pow(Tcanopy, 4);
        SensibleSub = coef_sensible * (Tupper - Tair);
        LatentSub = coef_latent * (esat_Tsub * rh_grnd - VPcanopy);
        GroundSub = NetShortGrnd + energy->AdvectSub - 
                                NetLongSub - SensibleSub - 
                                                     LatentSub;
    }
    // 复制到energy结构体
    cell->Ra_leaf = Ra_leaf;
    energy->Tupper = Tupper;
    energy->Tcanopy = Tcanopy;
    energy->NetLongSub = NetLongSub;
    energy->SensibleSub = SensibleSub;
    energy->LatentSub = LatentSub;
    energy->GroundSub = GroundSub;

    return (0);
}
