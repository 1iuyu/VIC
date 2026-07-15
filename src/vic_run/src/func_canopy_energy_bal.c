/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine iterates to determine the temperature of the canopy, and solve
 * the resulting fluxes between the canopy and the atmosphere and the canopy
 * and the ground.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    Calculate the canopy energy balance.
 *****************************************************************************/
int
func_canopy_energy_bal(size_t             hidx,
                       double             step_dt,
                       double             air_temp,
                       force_data_struct *force,
                       energy_bal_struct *energy,
                       cell_data_struct  *cell,
                       snow_data_struct  *snow,
                       soil_con_struct   *soil_con,
                       veg_var_struct    *veg_var,
                       veg_lib_struct    *veg_lib)
{
    extern parameters_struct param;
    extern option_struct     options;
    /* General Model Parameters */
    int         ErrorFlag;
    double      Ra_leaf;
    /* Initialization */
    double Qair = force->Qair[hidx];
    double wind = force->wind[hidx];
    double pressure = force->pressure[hidx];
    double longwave = force->longwave[hidx];
    double air_density = force->density[hidx];
    double theta_pot = force->theta_pot[hidx];
    double theta_v = force->theta_v[hidx];
    double daylen = force->daylen[hidx];
    double NetLAI = veg_var->NetLAI;
    double NetSAI = veg_var->NetSAI;
    double fcanopy = veg_var->fcanopy;
    double wetFrac = veg_var->wetFrac;
    double dryFrac = veg_var->dryFrac;
    double canopy_swq = veg_var->canopy_swq;
    double *T = energy->T;
    double *soil_T = cell->soil_T;
    double *Ra_over = cell->Ra_over;
    double *Ra_sub = cell->Ra_sub;
    double *Z0m_grnd = cell->Z0m_grnd;
    double *Z0m_sub = cell->Z0m_sub;
    double *displacement = cell->displacement;
    double *ref_height = cell->ref_height;
    double Qair_grnd = cell->Qair_grnd;
    double Canopy_Upper = veg_lib->Canopy_Upper;
    double Tcanopy = energy->Tcanopy;
    double Tfoliage = energy->Tfoliage;
    double Tstem = energy->Tstem;
    double Tgrnd = energy->Tgrnd;
    double SensibleLeaf = energy->SensibleLeaf;
    double SensibleStem = energy->SensibleStem;
    double LatentLeaf = energy->LatentLeaf;
    double NetShortSub = energy->NetShortSub;
    double NetLongOver = energy->NetLongOver;
    double EmissLongSub = energy->EmissLongSub;
    double EmissLongGrnd = energy->EmissLongGrnd;
    /********************************
      Vegetated canopy energy flux
    ********************************/
    size_t iter = 0;
    double ustar = 0.06;
    double wstar = 0.0;
    double wind_sub = 0.0;
    double h_leaf = 0.0;
    double temp_profile = 0.0;
    double Qair_profile = 0.0;
    double wind_over = 0.0;
    double delt_T = 1.0;
    double L_disp = 0.0;
    double Qair_over = 0.0;
    double corr_wind = 0.0;
    double Cs_turb = 0.0;
    double conv_h = 1000.0;
    double f_abs_stem;
    double RestTerm, coef_flux;
    double sa_leaf, sa_stem;
    double cp_leaf, cp_stem, Ra_stem;
    double stem_biomass;
    size_t moz_signchg_count = 0;
    double L_old = 0.0;
    double coef_latent = 0.0;
    double coef_sensible = 0.0;
    // 设置临时变量
    double wtgq0, wtlq0, wtaq0, wtl0, wta0, wtal, wtg;
    double wtstem0, wtgq, wtg0, wtalq, wtgaq, wtshi;
    double wtair, wtl, wtstem, wtga, wtaq, wtlq, wtsqi;

    // 计算生物质的热容量
    ErrorFlag = calc_biomass_heat(&f_abs_stem, 
                                  &sa_leaf,
                                  &sa_stem, 
                                  &cp_leaf, &cp_stem,
                                  &stem_biomass,
                                  &Ra_stem, 
                                  veg_var, veg_lib);
    if (ErrorFlag == ERROR) {
        return (ERROR);
    }
    // 调整冠层的空气动力学参数
    double lt = 0.0;
    double NetVEG = NetLAI + NetSAI;
    if (options.AERO_RESIST == AR_ZENG) {
        lt = min(NetVEG, 2.0);
        double egvf =(1.0 - exp(-lt)) / (1.0 - exp(-2.0));
        displacement[1] = egvf * displacement[1];
        Z0m_sub[0] = exp(egvf * log(Z0m_sub[0]) + (1.0 - egvf) * log(Z0m_grnd[0]));
    }
    else if (options.AERO_RESIST == AR_MEIER) {
        lt = max(NetVEG, 1.e-5);
        displacement[1] = Canopy_Upper * (1.0 - (1.0 - exp(-pow(7.5 * lt, 0.5)) / pow(7.5 * lt, 0.5)));
        lt = min(lt, veg_lib->Z0sub_LAImax);
        double delt = 2.0;
        double ustar_ini = pow(veg_lib->Z0sub_Cs + veg_lib->Z0sub_Cr * lt * 0.5, -0.5) * 
                                    veg_lib->Z0sub_c * lt * 0.25;
        double ustar1 = ustar_ini;

        while (delt > 1.e-4) {
            double ustar_prev = ustar1;
            ustar1 = ustar_ini * exp(ustar_prev);
            delt = fabs(ustar1 - ustar_prev);
        }
        ustar1 = 4.0 * ustar1 / lt / veg_lib->Z0sub_c;
        Z0m_sub[0] = Canopy_Upper * (1.0 - displacement[1] / Canopy_Upper) * exp(-CONST_KARMAN * ustar1 +
                      log(veg_lib->Z0sub_cw) - 1.0 + pow(veg_lib->Z0sub_cw, -1.0));
    }
    else {
        log_err("Unknown AERO_RESIST option");
    }
    Z0m_sub[1] = Z0m_sub[0];
    Z0m_sub[2] = Z0m_sub[1];
    ref_height[0] = param.REF_HEIGHT_WIND + Z0m_sub[0] + displacement[1];
    ref_height[1] = param.REF_HEIGHT + Z0m_sub[1] + displacement[1];
    ref_height[2] = param.REF_HEIGHT + Z0m_sub[2] + displacement[1];
    double thm = air_temp + 0.0098 * ref_height[1];
    // 设置未被积雪覆盖的植被比例
    double f_snow_veg;
    if (NetVEG > 0.05) {
        f_snow_veg = 1.0;
    }
    else {
        f_snow_veg = 0.0;
    }
    // 树冠和地面净吸收的长波辐射
    double coef_lw_atmos = EmissLongSub * (1.0 + (1.0 - EmissLongSub) * 
                        (1.0 - EmissLongGrnd)) * longwave;
    double coef_lw_grnd = EmissLongSub * EmissLongGrnd * CONST_STEBOL;
    double coef_lw_canopy = -(2.0 - EmissLongSub * (1.0 - EmissLongGrnd)) *
                                       EmissLongSub * CONST_STEBOL;
    // 计算比湿和梯度
    double qsatdT = 0.0;
    double qsat_T = 0.0;
    double esat_T = 0.0;
    svp_flags(Tfoliage, pressure,
                &esat_T, &qsat_T, 
                NULL, &qsatdT,
                ESAT| QSAT | QSDT);
    // 初始化植被/地面表面温度、树冠内部空气温度
    Tcanopy = (Tgrnd + thm) / 2.0; // 初始猜测
    Qair_over = (Qair + Qair_grnd) / 2.0;
    wind = max(wind, 0.01);
    double dth = thm - Tcanopy;
    double dqh = Qair - Qair_over;
    double deltq = Qair_grnd - Qair_over;
    double dthv = dth * (1.0 + 0.61 * Qair) + 0.61 * dqh * theta_pot;  // 虚位温差
    double zL_over = ref_height[0] - displacement[1];

    /* Initialize Obukhov length scale */
    L_disp = initialize_MOST(wind, dthv,
                                zL_over,
                                theta_v, Z0m_sub[0],
                                &corr_wind);
    // 初始化迭代变量
    double Tfoliage_init = Tfoliage;
    // double Tstem_init = Tstem;
    double canopyevap = 0.0;
    // 迭代求解冠层温度和相关的能量通量
    do {
        // Perform stability iteration
        ErrorFlag = FrictionVelocity(L_disp, corr_wind, 
                                     ref_height, &ustar,
                                     &temp_profile,
                                     &Qair_profile,
                                     Z0m_sub, displacement[1]);
        if (ErrorFlag == ERROR) {
            return (ERROR);
        }
        // 设置空气阻力
        Ra_over[0] = 1.0 / (ustar * ustar / corr_wind);
        Ra_over[1] = 1.0 / (temp_profile * ustar);
        Ra_over[2] = 1.0 / (Qair_profile * ustar);
        // 叶片的整片叶层阻力
        wind_over = corr_wind * sqrt(1.0 / (Ra_over[0] / corr_wind));
        wind_sub = min(0.4, 0.03 * corr_wind / ustar);
        h_leaf = 0.01 / (sqrt(wind_over) * sqrt(veg_lib->d_leaf));
        Ra_leaf = 1.0 / (h_leaf * wind_over);
        cell->Ra_leaf = Ra_leaf;
        // 计算下垫面土壤与冠层空气之间的湍流传送系数
        double w_frac = exp(-(NetVEG));
        double Cs_bare = CONST_KARMAN / 0.13 * pow((Z0m_sub[1] * wind_over / 1.5e-5), -0.45);
        double ri = (CONST_G * Canopy_Upper * (Tcanopy - Tgrnd)) / pow(Tcanopy * wind_over, 2.0);
        double Cs_dense = 0.004;
        if (Tcanopy - Tgrnd > 0.0) {
            Cs_dense = 0.004 / (1.0 + 0.5 * min(ri, 10.0));
        }
        Cs_turb = Cs_bare * w_frac + Cs_dense * (1.0 - w_frac);

        if (options.BIOMASST) {
            Ra_sub[1] = 1.0 / (Cs_turb * wind_sub);
        }
        else {
            Ra_sub[1] = 1.0 / (Cs_turb * wind_over);
        }
        Ra_sub[2] = Ra_sub[1];

        double vp_over = pressure * Qair_over * PA_PER_KPA / param.SVP_RDAIR;

        /**********************
          Stomatal Resistance
        **********************/
        photosynth_hydrostress(thm, daylen, esat_T,
                               qsat_T, vp_over, 
                               Qair_over, pressure, 
                               air_density, Tfoliage,
                               cell, soil_con,
                               veg_var, veg_lib);

        // 计算感热通量传导
        wtair = 1.0 / Ra_over[1];
        wtl = sa_leaf / Ra_leaf;
        wtg = 1.0 / Ra_sub[1];
        wtstem = sa_stem / (Ra_stem + Ra_leaf);
        wtshi  = 1.0 / (wtair + wtl + wtstem + wtg);
        wtl0 = wtl * wtshi;
        wtg0 = wtg * wtshi;
        wta0 = wtair * wtshi;
        wtstem0 = wtstem * wtshi;
        wtga = wta0 + wtg0 + wtstem0;
        wtal = wta0 + wtl0 + wtstem0;

        // 计算叶片潜在蒸发量的比例
        double rppdry = 0.0;
        if (dryFrac > 0.0) {
            rppdry = dryFrac * Ra_leaf * (veg_var->LAI_sun / (Ra_leaf + veg_var->RS_sunlit) + 
                            veg_var->LAI_sha / (Ra_leaf + veg_var->RS_shade)) / NetLAI;
        }
        else {
            rppdry = 0.0;
        }
        // 计算潜在蒸发量的比例
        double efpot = air_density * (NetVEG / Ra_leaf) * (qsat_T - Qair_over);
        double rpp = 0.0;
        if (efpot > 0.0) {
            if (cell->transp_fact > 0.0) {
                rpp = rppdry + wetFrac;
            } else {
                rpp = wetFrac; // No transpiration if transp_fact below 1.e-10
            }
            // Check total evapotranspiration from leaves
            rpp = min(rpp, (cell->transp + canopy_swq / step_dt) / efpot);
        }
        else {
            // No transpiration if potential evaporation less than zero
            rpp = 1.0;
        }
        wtaq = f_snow_veg / Ra_over[2];
        wtlq = f_snow_veg * NetVEG / Ra_leaf * rpp;
        double fsno_dl = snow->snow_depth / 0.05;    // effective snow cover for (dry)plant litter
        double elai_dl = 0.5 * (1.0 - min(fsno_dl, 1.0)); // exposed (dry)litter area index
        double rdl = (1.0 - exp(-elai_dl)) / (0.004 * wind_over); // dry litter layer resistance
        if (deltq < 0.0) {
            wtgq = f_snow_veg / (Ra_sub[2] + rdl);
        }
        else {
            wtgq = f_snow_veg / (Ra_sub[2] + cell->Ra_evap);
        }
        wtsqi = 1.0 / (wtaq + wtlq + wtgq);

        wtgq0 = wtgq * wtsqi;      // ground
        wtlq0 = wtlq * wtsqi;      // leaf
        wtaq0 = wtaq * wtsqi;      // air
        wtgaq = wtaq0 + wtgq0;     // air + ground
        wtalq = wtaq0 + wtlq0;     // air + leaf

        coef_sensible = air_density * CONST_CPDAIR * wtl;
        coef_latent = CONST_LATVAP * air_density * wtlq;

        SensibleLeaf = coef_sensible * (wtga * Tfoliage - wtg0 * Tgrnd - wta0 * thm - wtstem0 * Tstem);
        SensibleStem = air_density * CONST_CPDAIR * wtstem * ((wta0 + wtg0 + wtl0) * 
                                    Tstem - wtg0 * Tgrnd - wta0 * thm - wtl0 * Tfoliage);
        
        LatentLeaf = coef_latent * (wtgaq * qsat_T - wtgq0 * Qair_grnd - wtaq0 * Qair);

        NetLongOver = coef_lw_atmos + coef_lw_canopy * pow(Tfoliage, 4.0);
        RestTerm = (1.0 - f_abs_stem) * (NetShortSub + NetLongOver + coef_lw_grnd * pow(Tgrnd, 4.0)) -
                    SensibleLeaf - LatentLeaf - (cp_leaf / step_dt) * (Tfoliage - Tfoliage_init);
        coef_flux = (1.0 - f_abs_stem) * (-4.0 * coef_lw_canopy * pow(Tfoliage, 3.0)) +
                    coef_sensible * wtga + coef_latent * wtgaq * qsatdT + cp_leaf / step_dt;
        delt_T = RestTerm / coef_flux;
        Tfoliage += delt_T;

        efpot = air_density * (NetVEG / Ra_sub[1]) *
                (wtgaq * (qsat_T + qsatdT * delt_T) -
                 wtgq0 * Qair_grnd - wtaq0 * Qair);
        canopyevap = rpp * efpot;
        double res_energy = max(0.0, canopyevap - cell->transp - canopy_swq / step_dt);
        canopyevap = min(canopyevap, cell->transp + canopy_swq / step_dt);

        // Update 
        SensibleLeaf += coef_sensible * wtga * delt_T + CONST_LATVAP * res_energy;
        SensibleStem += air_density * CONST_CPDAIR * wtstem * (-wtl0 * delt_T);

        svp_flags(Tfoliage, pressure, 
                  NULL, &qsat_T, NULL, 
                  &qsatdT, QSAT | QSDT);
        // 更新冠层温度和水汽含量的加权平均值
        Tcanopy = wtg0 * Tgrnd + wta0 * thm + wtl0 * Tfoliage + wtstem0 * Tstem;
        Qair_over = wtlq0 * qsat_T + Qair * wtaq0 + wtgq0 * Qair_grnd;

        dth = thm - Tcanopy;
        dqh = Qair - Qair_over;
        deltq = wtalq * Qair_grnd - wtlq0 * qsat_T - wtaq0 * Qair;

        double tstar = temp_profile * dth;
        double qstar = Qair_profile * dqh;

        double thvstar = tstar * (1.0 + 0.61 * Qair) + 0.61 * theta_pot * qstar;
        double zeta = zL_over * CONST_KARMAN * CONST_G * thvstar / (pow(ustar, 2.0) * theta_v);

        if (zeta >= 0.0) {
            zeta = min(2.0, max(zeta, 0.01));
            corr_wind = max(wind, 0.1);
        }
        else {
            zeta = max(-100.0, min(zeta, -0.01));
            if (ustar * thvstar > 0.0) {
                log_warn("ustar * thvstar is positive and has to be negative");
                wstar = 0.0;
            }
            else {
                wstar = pow(-CONST_G * ustar * thvstar * conv_h / theta_v, 0.333);
            }
            corr_wind = sqrt(wind * wind + wstar * wstar);
        }
        L_disp = zL_over / zeta;

        if (L_old * L_disp < 0.0) {
            moz_signchg_count++;
        }
        if (moz_signchg_count >= 4) {
            L_disp = zL_over / (-0.01);
        }
        L_old = L_disp;
        
        iter++;
    }
    while (iter < param.MAX_ITER_OVER && fabs(delt_T) > 0.01);
    
    /********************************
      Vegetated surface energy flux
    ********************************/
    energy->error = (1.0 - f_abs_stem) * (NetShortSub + coef_lw_atmos + 
            coef_lw_canopy * pow(Tfoliage, 3.0) * (Tfoliage + 4.0 * 
            delt_T) + coef_lw_grnd * pow(Tgrnd, 4.0)) - SensibleLeaf - CONST_LATVAP * 
            canopyevap - ((Tfoliage - Tfoliage_init) * cp_leaf / step_dt);
    double dt_stem = 0.0;
    if (options.BIOMASST) {
        if (stem_biomass > 0.0) {
            dt_stem = (f_abs_stem * (NetShortSub + coef_lw_atmos + coef_lw_canopy * 
                pow(Tfoliage_init, 4.0) +
                coef_lw_grnd * pow(Tgrnd, 4.0)) - SensibleStem) / (cp_stem / step_dt -
                f_abs_stem * coef_lw_canopy * 4.0 * pow(Tfoliage_init, 3.0));
        } else {
            dt_stem = 0.0;
        }
        Tstem += dt_stem;
    } else {
        dt_stem = 0.0;
    }
    // compute individual sensible heat fluxes
    double delt = wtal * Tgrnd - wtl0 * Tfoliage - wta0 * thm - wtstem0 * Tstem;
    double delt_snow = wtal * T[0] - wtl0 * Tfoliage - wta0 * thm - wtstem0 * Tstem;
    double delt_soil = wtal * soil_T[0] - wtl0 * Tfoliage - wta0 * thm - wtstem0 * Tstem;
    energy->sensible += fcanopy * (CONST_CPDAIR * air_density * wtg * delt - energy->sensible);
    energy->SensibleSnow += fcanopy * (CONST_CPDAIR * air_density * 
                            wtg * delt_snow - energy->SensibleSnow);
    energy->SensibleSoil += fcanopy * (CONST_CPDAIR * air_density * 
                            wtg * delt_soil - energy->SensibleSoil);
    // compute individual latent heat fluxes
    double delq_snow = wtalq * cell->Qair_snow - wtlq0 * qsat_T - wtaq0 * Qair;
    double delq_soil = wtalq * cell->Qair_soil - wtlq0 * qsat_T - wtaq0 * Qair;
    energy->latent += fcanopy * (air_density * wtgq * deltq * CONST_LATVAP - energy->latent);
    energy->LatentSnow += fcanopy * (air_density * wtgq * delq_snow * CONST_LATVAP - energy->LatentSnow);
    energy->LatentSoil += fcanopy * (air_density * wtgq * delq_soil * CONST_LATVAP - energy->LatentSoil);
    cell->esoil += fcanopy * (air_density * wtgq * deltq - cell->esoil);
    energy->deriv_evap += fcanopy * (air_density * wtgq * cell->Qair_grnd * 
                    CONST_G / (CONST_RWV * Tgrnd) / CONST_RHOFW - energy->deriv_evap);

    // 更新累积的露水 (kg/m2)
    if (Tfoliage > CONST_TKFRZ) {
        if ((canopyevap - cell->transp) * step_dt > veg_var->iter_intrain) {
            veg_var->int_snow = max(0.0, veg_var->iter_intsnow + veg_var->iter_intrain +
                     (cell->transp - canopyevap) * step_dt);
        }
        veg_var->int_rain = max(0.0, veg_var->iter_intrain + (cell->transp - canopyevap) * step_dt);
    }
    else if (Tfoliage <= CONST_TKFRZ) {
        if ((canopyevap - cell->transp) * step_dt > veg_var->iter_intsnow) {
            veg_var->int_rain = veg_var->iter_intrain + veg_var->iter_intsnow + (cell->transp - 
                                    canopyevap) * step_dt;
        }
        veg_var->int_snow = max(0.0, veg_var->iter_intsnow + (cell->transp - canopyevap) * step_dt);
    }
    veg_var->canopy_swq = veg_var->int_rain + veg_var->int_snow;
    
    // 复制到energy结构体
    cell->Ra_leaf = Ra_leaf;
    cell->Ra_stem = Ra_stem;
    energy->Tcanopy = Tcanopy;
    energy->Tstem = Tstem;
    energy->Tfoliage = Tfoliage;
    cell->canopyevap = canopyevap;
    energy->SensibleStem = SensibleStem;
    energy->SensibleLeaf = SensibleLeaf;
    energy->NetLongOver = NetLongOver;
    coef_sensible = air_density * CONST_CPDAIR * wtg * wtal;
    coef_latent = air_density * wtgq * wtalq * qsatdT * energy->LatentVapOver;
    double deriv_sub = -(coef_sensible + coef_latent + 4.0 * 
                        EmissLongSub * CONST_STEBOL * pow(Tgrnd, 3.0));
    energy->deriv_terms += fcanopy * (deriv_sub - energy->deriv_terms);
    energy->Tsurf = EmissLongSub * Tfoliage + (1.0 - EmissLongSub) * sqrt(sqrt(pow(Tgrnd, 4.0)));

    double LongSubIn = (1.0 - EmissLongSub) * EmissLongGrnd * longwave + EmissLongSub * 
                EmissLongGrnd * CONST_STEBOL * pow(Tfoliage, 4.0) * (1.0 - f_abs_stem) + 
                EmissLongSub * EmissLongGrnd * CONST_STEBOL * pow(Tstem, 4.0) * f_abs_stem;

    energy->longwave += fcanopy * (LongSubIn - EmissLongGrnd * CONST_STEBOL * 
                        pow(Tgrnd, 4.0) - energy->longwave);
    energy->NetLongSnow += fcanopy * (LongSubIn - EmissLongGrnd * CONST_STEBOL * 
                        pow(T[0], 4.0) - energy->NetLongSnow);
    energy->NetLongSoil += fcanopy * (LongSubIn - EmissLongGrnd * CONST_STEBOL * 
                        pow(soil_T[0], 4.0) - energy->NetLongSoil);

    return (0);
}
