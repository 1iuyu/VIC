/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine computes the surface energy balance for bare soil and
 * vegetation uncovered by snow.  It computes outgoing longwave, sensible heat
 * flux, ground heat flux, and storage of heat in the thin upper layer, based
 * on the given surface temperature.
 *
 * The Energy Balance Equation used comes from Xu Liang's Paper "Insights of
 * the Ground Heat Flux in Land Surface Parameterization Schemes."
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    Calculate the surface energy balance.
 *****************************************************************************/
int
func_surf_energy_bal(double             air_density,
                     double             Ra_evap,
                     double             longwave,
                     double             air_temp,
                     double             rh_grnd,
                     double             theta,          // potential temperature (Kelvin)
                     double             vp,
                     double             pressure,
                     double             Qair,
                     double             wind,
                     double             ref_height,
                     double            *roughness,
                     double            *displacement,
                     energy_bal_struct *energy,
                     cell_data_struct  *cell,
                     snow_data_struct  *snow,
                     soil_con_struct   *soil_con)
{
    extern parameters_struct param;
    extern option_struct     options;

    /* meteorological forcing terms */
    size_t       lindex;
    double       Error;
    double       ErrorFlag;
    double       SVP_liq;
    double       SVP_ice;
    double       liq_slope;
    double       ice_slope;
    double       esat_Tgrnd;
    double       esat_slope;
    double       coef_latent;
    double       layer_T;
    double       Qair_surf;
    double       coef_sensible;
    double       NetLongGrnd;
    double       SensibleGrnd;
    double       LatentGrnd;
    double       GroundGrnd;
    double       RestTerm;
    double       coef_flux;
    double       layer_depth;

    /***************
      MAIN ROUTINE
    ***************/
    Error = 0;
    double Tlower = energy->Tlower;
    double PsyCh_grnd = energy->PsyCh_grnd;
    double AdvectGrnd = energy->AdvectGrnd;
    double NetShortGrnd = energy->NetShortGrnd;
    double EmissLongGrnd = energy->EmissLongGrnd;
    double *Ra_grnd = cell->Ra_grnd;
    double *kappa_node = energy->kappa_node;

    if (snow->Nsnow > 0) {
        lindex = snow->Nsnow - 1; // 存在雪层
        layer_depth = snow->dz_snow[lindex];
        layer_T = snow->pack_T[0];
    }
    else {
        layer_depth = soil_con->dz_soil[0]; // 不存在雪层，kappa_node第一层为土壤热力学参数
        layer_T = cell->soil_T[0];
    }
    double coef_longwave = EmissLongGrnd * CONST_BOLTZ;
    double coef_ground = 2.0 * kappa_node[0] / layer_depth;
    size_t iter = 0;
    double delta_grnd = param.TOL_A + 1.;
    double corr_wind = 0.0;
    /* virtual potential temperature difference between ground and air */
    double dth = air_temp + 0.0098 * param.REF_HEIGHT_WIND - Tlower;
    double dqh = Qair / Qair_surf;
    double theta_v = theta * (1.0 + 0.61 * Qair);
    double dthv = dth + (1 + 0.61 * Qair) + 0.61 * dqh * theta;
    double zL_grnd = param.REF_HEIGHT_WIND;
    double ustar = 0.06;
    double tstar = 0.0;
    double qstar = 0.0;
    double thvstar = 0.0;
    double zeta = 0.0;
    double wstar = 0.0;
    double conv_h = 1000.0;  // convective boundary height [m]
    double conv_beta = 1.0;  // coefficient of convective velocity
    double Qair_profile = 0.;
    double temp_profile = 0.;
    
    double sensible_grnd = 0.0;
    double L_disp = initialize_MOST(wind, dthv,
                                    param.REF_HEIGHT_WIND, 
                                    theta, roughness[0],
                                    &corr_wind);

    do {
        ErrorFlag = FrictionVelocity(wind, L_disp, corr_wind, &ustar,
                                     &temp_profile, &Qair_profile,
                                     roughness[0], displacement[0]);
        if (ErrorFlag == ERROR) {
            return (ERROR);
        }
        tstar = temp_profile * dth;
        qstar = Qair_profile * dqh;
        thvstar = tstar * (1.0 + 0.61 * Qair) + 0.61 * theta * qstar;
        zeta = zL_grnd * CONST_KARMAN * CONST_G * thvstar / (pow(ustar, 2.0) * theta_v);

        if (zeta >= 0.0) {
            zeta = max(zeta, 0.01);
            corr_wind = max(wind, 0.1);
        }
        else {
            zeta = max(-100.0, min(zeta,-0.01));
            wstar = conv_beta * pow(-CONST_G * ustar * thvstar * conv_h / theta_v, 0.333);
            corr_wind = sqrt(wind * wind + wstar * wstar);
        }
        L_disp = zL_grnd / zeta;

        // Determine aerodynamic resistances
        Ra_grnd[0] = 1.0 / (ustar * ustar / corr_wind);
        Ra_grnd[1] = 1.0 / (temp_profile * ustar);
        Ra_grnd[2] = 1.0 / (Qair_profile * ustar);

        /* Saturated vapor pressure at Tgrnd */
        svp(Tlower, &SVP_liq, &SVP_ice);
        svp_slope(Tlower, &liq_slope, &ice_slope);
        if (Tlower > CONST_TKFRZ) {
            esat_Tgrnd = SVP_liq;
            esat_slope = liq_slope;
        }
        else {
            esat_Tgrnd = SVP_ice;
            esat_slope = ice_slope;
        }
        /* ground fluxes and temperature change */
        coef_sensible = CONST_CPDAIR * air_density / Ra_grnd[1];
        coef_latent = air_density * CONST_CPDAIR /
                            PsyCh_grnd / (Ra_grnd[2] + Ra_evap);
        NetLongGrnd = calc_longwave(Tlower, coef_longwave) - EmissLongGrnd * longwave;
        SensibleGrnd = calc_sensible_heat(Tlower, air_temp, coef_sensible);
        LatentGrnd = calc_latent_heat(esat_Tgrnd, rh_grnd, vp, coef_latent);
        GroundGrnd = calc_ground_heat(Tlower, layer_T, coef_ground);

        RestTerm = NetShortGrnd - NetLongGrnd - SensibleGrnd -
                         LatentGrnd - GroundGrnd + AdvectGrnd;
        coef_flux = 4.0 * coef_longwave * pow(Tlower, 3) + coef_sensible + 
                                         coef_latent * esat_slope + coef_ground;
        delta_grnd = RestTerm / coef_flux;

        // 更新通量（线性近似）
        NetLongGrnd += 4.0 * coef_longwave * pow(Tlower, 3) * delta_grnd;
        SensibleGrnd += coef_sensible * delta_grnd;
        LatentGrnd += coef_latent * esat_slope * delta_grnd;
        GroundGrnd += coef_ground * delta_grnd;

        /* update ground temperature */
        Tlower += delta_grnd;

        /* for computing M-O length */
        sensible_grnd = coef_sensible * (Tlower - air_temp);

        /* update specific humidity */
        svp(Tlower, &SVP_liq, &SVP_ice);
        if (Tlower > CONST_TKFRZ) {
            esat_Tgrnd = SVP_liq;
        }
        else {
            esat_Tgrnd = SVP_ice;
        }
        Qair_surf = 0.622 * (esat_Tgrnd * rh_grnd) / 
                            (pressure - 0.378 * (esat_Tgrnd * rh_grnd));
        double moisture_flux = (Qair_surf - Qair) * coef_latent *
                                         PsyCh_grnd / CONST_CPDAIR;
        iter++;
    }
    while (iter < param.MAX_ITER_MOST && fabs(delta_grnd) > 0.01);    // iter 5

    /* if snow on grnd and Tgrnd>0.0, reset Tgrnd=0.0 */
    if (snow->snow_depth > 0.05 && Tlower > CONST_TKFRZ) {

        Tlower = CONST_TKFRZ;
        /* 重新计算饱和水汽压 */
        svp(Tlower, &SVP_liq, &SVP_ice);
        if (Tlower > CONST_TKFRZ) {
            esat_Tgrnd = SVP_liq;
        }
        else {
            esat_Tgrnd = SVP_ice;
        }
        /* update ground fluxes */        
        NetLongGrnd = calc_longwave(Tlower, coef_longwave) - 
                                        EmissLongGrnd * longwave;
        SensibleGrnd = coef_sensible * (Tlower - air_temp);
        LatentGrnd = coef_latent * (esat_Tgrnd * rh_grnd - vp);
        GroundGrnd = NetShortGrnd + AdvectGrnd - 
                                NetLongGrnd - SensibleGrnd - 
                                                     LatentGrnd;
    }
    // 复制到energy结构体
    energy->Tlower = Tlower;
    energy->NetLongGrnd = NetLongGrnd;
    energy->SensibleGrnd = SensibleGrnd;
    energy->LatentGrnd = LatentGrnd;
    energy->GroundGrnd = GroundGrnd;

    return (Error);
}
