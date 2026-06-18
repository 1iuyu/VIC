/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine computes the surface energy balance for bare soil. It 
 * computes sensible heat and latent heat fluxes and their derivatives based 
 * on the given surface temperature.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    Calculate the surface energy balance.
 *****************************************************************************/
int
func_surf_energy_bal(size_t             hidx,
                     double             air_temp,
                     force_data_struct *force,
                     energy_bal_struct *energy,
                     cell_data_struct  *cell,
                     snow_data_struct  *snow)
{
    extern parameters_struct param;
    extern option_struct options;
    /***************
      MAIN ROUTINE
    ***************/
    int ErrorFlag;
    double Tgrnd = energy->Tgrnd;
    double *Ra_grnd = cell->Ra_grnd;
    double *Z0m_grnd = cell->Z0m_grnd;
    double *ref_height = cell->ref_height;
    double *displacement = cell->displacement;
    double Qair = force->Qair[hidx];
    double wind = force->wind[hidx];
    double pressure = force->pressure[hidx];
    double longwave = force->longwave[hidx];
    double air_density = force->density[hidx];
    double theta_pot = force->theta_pot[hidx];
    double theta_v = force->theta_v[hidx];
    double corr_wind = 0.0;
    double ustar = 0.06;
    double wstar = 0.0;
    double last_ustar = 0.0;
    double last_Ldisp = 0.0;
    double conv_h = 1000.0;  // convective boundary height [m]
    double conv_beta = 1.0;  // coefficient of convective velocity
    double Qair_profile = 0.;
    double temp_profile = 0.;
    double L_disp = 0.0;  // Monin-Obukhov长度
    size_t moz_signchg_count = 0;
    /* 计算地表比湿 */
    ErrorFlag = calc_surf_humidity(pressure, 
                                   Qair, Tgrnd, 
                                   snow, cell);
    if (ErrorFlag == ERROR) {
        return (ERROR);
    }
    // 初始设置参考高度
    ref_height[0] = param.REF_HEIGHT_WIND + Z0m_grnd[0] + displacement[0]; // 动量
    ref_height[1] = param.REF_HEIGHT + Z0m_grnd[1] + displacement[0];  // 感热
    ref_height[2] = param.REF_HEIGHT + Z0m_grnd[2] + displacement[0];  // 潜热
    double thm = air_temp + 0.0098 * ref_height[1];
    double dth = thm - Tgrnd;
    double dqh = Qair - cell->Qair_grnd;  // 比湿差
    double dthv = dth * (1 + 0.61 * Qair) + 0.61 * dqh * theta_pot;  // 虚位温差
    double zL_grnd = ref_height[0];
    /* Initialize Obukhov length scale */
    L_disp = initialize_MOST(wind, dthv,
                             zL_grnd,
                             theta_v, Z0m_grnd[0],
                             &corr_wind);

    int iterT = 0;
    do {

        last_ustar = ustar;
        // Perform stability iteration
        ErrorFlag = FrictionVelocity(L_disp, corr_wind,
                                     ref_height, &ustar,
                                     &temp_profile, 
                                     &Qair_profile,
                                     Z0m_grnd, displacement[0]);
        if (ErrorFlag == ERROR) {
            return (ERROR);
        }

        double tstar = temp_profile * dth;
        double qstar = Qair_profile * dqh;
        double thvstar = tstar * (1.0 + 0.61 * Qair) + 0.61 * theta_pot * qstar;
        double zeta = zL_grnd * CONST_KARMAN * CONST_G * thvstar / 
                                            (pow(ustar, 2.0) * theta_v);

        // 更新潜热和感热的粗糙度长度
        if (options.AERO_RESIST == AR_ZENG) {
            Z0m_grnd[1] = Z0m_grnd[0] / exp(0.13 * pow(ustar * Z0m_grnd[0] / 1.5e-5, 0.45));
        }
        else if (options.AERO_RESIST == AR_MEIER) {
            Z0m_grnd[1] = param.ROUGH3 * param.ROUGH_NU / ustar * 
                        exp(-param.ROUGH_BETA * pow(ustar, 0.5) * pow(fabs(tstar), 0.25));
        }
        else {
            log_err("Unknown AERO_RESIST option");
        }
        // 假设潜热和感热的粗糙度长度相同
        Z0m_grnd[2] = Z0m_grnd[1];
        // 更新参考高度
        ref_height[0] = param.REF_HEIGHT_WIND + Z0m_grnd[0] + displacement[0];
        ref_height[1] = param.REF_HEIGHT + Z0m_grnd[1] + displacement[0];
        ref_height[2] = param.REF_HEIGHT + Z0m_grnd[2] + displacement[0];

        if (zeta >= 0.0) {
            zeta = min(2.0, max(zeta, 0.01));
            corr_wind = max(wind, 0.1);
        }
        else {
            zeta = max(-100.0, min(zeta,-0.01));
            wstar = conv_beta * pow(-CONST_G * ustar * thvstar * conv_h / theta_v, 0.333);
            corr_wind = sqrt(wind * wind + wstar * wstar);
        }
        L_disp = zL_grnd / zeta;
        // 检查L_disp的符号变化以判断MOST迭代是否收敛
        if (last_Ldisp * L_disp < 0.0) {
            moz_signchg_count++;
        }
        if (moz_signchg_count >= 4) {
            log_err("MOST iteration failed to converge after 4 sign changes in L_disp");
        }
        last_Ldisp = L_disp;

        iterT++;
    } 
    while (iterT < param.MAX_ITER_MOST && fabs(ustar - last_ustar) > param.TOL_A);

    // Determine aerodynamic resistances
    Ra_grnd[0] = 1.0 / (ustar * ustar / corr_wind);
    Ra_grnd[1] = 1.0 / (temp_profile * ustar);
    Ra_grnd[2] = 1.0 / (Qair_profile * ustar);

    double esat_Tgrnd = 0.0;
    double dew_point = 0.0;
    double coef_latent = 0.0;
    double coef_longwave = energy->EmissLongGrnd * CONST_STEBOL;
    double coef_sensible = CONST_CPDAIR * air_density / Ra_grnd[1];
    // Get saturated vapor pressure at forcing height
    svp_flags(Tgrnd, pressure,
              &esat_Tgrnd, NULL,
              NULL, NULL, ESAT);
    double vp = max(Qair * pressure / (0.622 + Qair), esat_Tgrnd * 0.01);
    // 计算露点温度
    if (Tgrnd < CONST_TKFRZ) {
        dew_point = 273.86 * log(vp / 611.21) / (22.587 - log(vp / 611.21));
    }
    else {
        dew_point = 243.04 * log(vp / 610.94) / (17.625 - log(vp / 610.94));
    }
    dew_point += CONST_TKFRZ;
    if (dqh > 0.0) {
        if (Tgrnd > dew_point) {
            coef_latent = 0.0;
        } else {
            coef_latent = air_density / Ra_grnd[2];
        }
    } else {
        coef_latent = air_density / (Ra_grnd[2] + cell->Ra_evap);
    }
    // 计算分子扩散最小通量（静稳空气限制）
    double sensible_min = CONST_KDAIR / dth / (ref_height[1] - displacement[1]);
    double vapor_min = CONST_VAPDIFF * dqh * air_density / (ref_height[2] - displacement[2]);
    double SensibleGrnd = -coef_sensible * dth;
    double VaporGrnd = -coef_latent * dqh;

    if (fabs(SensibleGrnd) < sensible_min) {
        SensibleGrnd = sensible_min;
        coef_sensible = CONST_KDAIR / (ref_height[1] - displacement[1]);
    }
    if (fabs(VaporGrnd) < vapor_min) {
        VaporGrnd = vapor_min;
        coef_latent = CONST_VAPDIFF * air_density;
    }
    // 计算地表的潜热和长波辐射
    cell->esoil_grnd = VaporGrnd;
    energy->LatentGrnd = VaporGrnd * energy->LatentVapGrnd;
    energy->NetLongGrnd = energy->EmissLongGrnd * longwave - coef_longwave * pow(Tgrnd, 4.0);
    energy->SensibleGrnd = SensibleGrnd;
    energy->deriv_grnd = -(coef_sensible + 4.0 * coef_longwave * pow(Tgrnd, 3) +
                         coef_latent * energy->LatentVapGrnd * cell->Qair_deriv);
    energy->deriv_egrnd = -coef_latent * cell->Qair_grnd * CONST_G /
                                                (CONST_RWV * Tgrnd) / CONST_RHOFW;
                        
    return (0);
}
