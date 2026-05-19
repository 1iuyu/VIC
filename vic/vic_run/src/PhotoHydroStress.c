/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine computes canopy stomatal resistance and foliage photosynthesis
 * based on Ball-Berry  and Jarvis-Montanari scheme
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    Calculate canopy photosynthesis and stomatal resistance.
 *****************************************************************************/
int
PhotoHydroStress(double             air_temp,
                 double             esat_T,
                 double             qsat_T,
                 double             vp_over,
                 double             pressure,
                 energy_bal_struct *energy,
                 cell_data_struct  *cell,
                 soil_con_struct   *soil_con,
                 veg_var_struct    *veg_var,
                 veg_lib_struct    *veg_lib)
{   
    extern option_struct     options;
    extern parameters_struct param;
    size_t i, Nsoil, Ncanopy;
    double tmp_aPAR = 0.;
    double aPAR_sunlit = energy->aPAR_sunlit;
    double aPAR_shade = energy->aPAR_shade;
    double fcanopy = veg_var->fcanopy;
    double Tfoliage = energy->Tfoliage;
    double *dz_soil = soil_con->dz_soil;
    double *root = cell->root;
    double *LAI_z = veg_var->LAI_z;
    double kp_sun[MAX_CANOPYS];
    double kp_sha[MAX_CANOPYS];
    double par_sun[MAX_CANOPYS];
    double par_sha[MAX_CANOPYS];
    double lmr_sun[MAX_CANOPYS];
    double lmr_sha[MAX_CANOPYS];
    double vcmax_sun[MAX_CANOPYS];
    double jmax_sun[MAX_CANOPYS];
    double tpu_sun[MAX_CANOPYS];
    double vcmax_sha[MAX_CANOPYS];
    double jmax_sha[MAX_CANOPYS];
    double tpu_sha[MAX_CANOPYS];
    // 初始化变量
    for (i = 0; i < MAX_CANOPYS; i++) {
        kp_sun[i] = 0.0;
        kp_sha[i] = 0.0;
        par_sun[i] = 0.0;
        par_sha[i] = 0.0;
        lmr_sun[i] = 0.0;
        lmr_sha[i] = 0.0;
        vcmax_sun[i] = 0.0;
        jmax_sun[i] = 0.0;
        tpu_sun[i] = 0.0;
        jmax_sha[i] = 0.0;
        tpu_sha[i] = 0.0;
    }

    double RGL = veg_lib->RGL;
    double rmax = veg_lib->rmax;
    double rmin = veg_lib->rmin;
    double m_vpd = veg_lib->m_vpd;
    double T_opt_trans = veg_lib->T_opt_trans;

    // initialize variables
    Nsoil = cell->Nsoil;
    Ncanopy = cell->Ncanopy;
    double f_N = 0.0;
    double CF = pressure / (8.314 * air_temp) * param.MAX_LIMIT;
    double RS_tmp = 1.0 / cell->Ra_leaf * CF;
    double tmp_photoleaf = 0.;
    double atmosO2 = pressure * 0.209;
    double atmosCO2 = pressure * 395.0e-06;
    double carbonylmax = 0.0;
    double indexgrow = 0.0;
    double rho_root = 0.0;
    double froot_carbon = 20.0; // kg/m2
    double lmrhd = 150650.0; // J/mol
    double lmrse = 490.0; // J/mol/K
    double froot_leaf = 0.0;
    double NetLAI = veg_var->NetLAI;
    double NetSAI = veg_var->NetSAI;
    double *matric = cell->matric;
    double *Zsum_soil = soil_con->Zsum_soil;
    double *conductivity = cell->conductivity;
    double lmrc = fth25(lmrhd, lmrse);
    double root_radius = 0.29e-03; // m
    double root_density = 0.31e06; // (g biomass / m3 root)
    double hk_total = 0.0;
    double matric_50 = 0.0; // ???
    double krmax = 0.0; // ???
    double daylen_fact = 0.0;
    // 计算常数：根横截面积 (m2)
    double xsec_root = CONST_PI * root_radius * root_radius;
    
    for (i = 0; i < Nsoil; i++) {
        // 根生物量密度：g biomass / m^3 soil
        rho_root = 2.0 * froot_carbon * root[i] / dz_soil[i];
        // 确保最小根生物量 (1 gC/m^2)
        rho_root = max(2.0, rho_root);
        // 根长密度：m root / m^3 soil
        double rlen_dens = rho_root / (xsec_root * root_density);
        // 根面积指数 (RAI)
        double rai = (NetLAI + NetSAI) * froot_leaf * root[i];
        // 粗根长度 - 使用指定的侧向长度
        double croot_length = 0.25; // specified lateral coarse root length [m]
        double r_soil = sqrt(1.0 / (CONST_PI * rlen_dens));
        // 土壤导水率 (m/s)
        double hk_soil = conductivity[i] / r_soil;
        // 使用植被PLC函数调整根区导水率
        // fs: 由于根水势降低（更负）导致的电导减少因子
        double fs = plc(matric[i], matric_50);
        // 根电导：单位面积单位长度的电导 (m/s)
        double hk_root = (fs * rai * krmax) / (croot_length + Zsum_soil[i]);
        hk_soil = max(hk_soil, 1.e-16);
        hk_root = max(hk_root, 1.e-16);
        // 计算土壤和根的总阻力，然后取倒数得到总电导
        double Ra_total = 1.0 / hk_soil + 1.0 / hk_root;
        // 电导是阻力的倒数
        // 对于表层土壤，明确设置电导为0
        if (rai * root[i] > 0.0 && i > 0) {
            hk_total = 1.0 / Ra_total;
        } 
        else {
            hk_total = 0.0;
        }
    }
    // Miscellaneous parameters
    double kc25_coef = 404.9e-6; // mol/mol
    double ko25_coef = 278.4e-3; // mol/mol
    double cp25_yr = 42.75e-6; // mol/mol
    double cpha = 37830.0; // Activation energy for cp [J/mol]
    double koha = 36380.0; // Activation energy for ko [J/mol]
    double kcha = 79430.0; // Activation energy for kc [J/mol]
    double lmrha = 46390.0; // Activation energy for lmr [J/mol]
    double sco  = 0.5 * 0.209 / cp25_yr;
    double cp25 = 0.5 * atmosO2 / sco;
    double KC = kc25_coef * pressure * ft(Tfoliage, kcha);
    double KO = ko25_coef * pressure * ft(Tfoliage, koha);
    double CP = cp25 * ft(Tfoliage, cpha);
    // 从veg_lib参数获取
    double leaf_CN = 0.0; // ??? // leaf C:N (gC/gN)
    double SLA_top = 0.0; // ??? // specific leaf area at top of canopy [m2/gC]
    double flnr = 0.0; // ??? // fraction of leaf N in Rubisco enzyme (gN Rubisco/gN leaf)
    double fnr = 7.16; // Mass ratio: gRubisco/gN in Rubisco
    double tpu25ratio = 0.167;
    double kp25ratio = 20000.0;
    double vcmaxha = 72000.0; // Activation energy for vcmax [J/mol]
    double jmaxha = 50000.0; // Activation energy for jmax [J/mol]
    double tpuha = 72000.0; // Activation energy for tpu [J/mol]
    double vcmaxhd = 200000.0; // Deactivation energy for vcmax [J/mol]
    double jmaxhd = 200000.0; // Deactivation energy for jmax [J/mol]
    double tpuhd = 200000.0; // Deactivation energy for tpu [J/mol]
    double vcmaxse_sf = 1.0; // Scale factor for vcmaxse
    double jmaxse_sf = 1.0; // Scale factor for jmaxse
    double tpuse_sf = 1.0; // Scale factor for tpuse
    double fnps = 0.15; // Fraction of light absorbed by non-photosynthetic pigment
    // 叶片氮浓度
    double lnc = 1.0 / (leaf_CN * SLA_top);
    lnc = min(lnc, 10.0);
    double vcmax25 = lnc * flnr * 7.16 * 60.0 * daylen_fact;
    double jmax25 = ((2.59 - 0.035 * min(max((Tfoliage - CONST_TKFRZ), 11.0), 35.0)) * vcmax25); // ???
    double tpu25 = tpu25ratio * vcmax25;
    double kp25 = kp25ratio * vcmax25;
    double luvcmax25 = vcmax25;
    double lujmax25t = jmax25;
    double lutpu25 = tpu25;
    // 氮衰减系数
    if (daylen_fact == 0.0) {
        f_N = 0.0;
    }
    else {
        f_N = exp(0.00963 * vcmax25 / daylen_fact - 2.43);
    }
    // Leaf maintenance respiration in proportion to vcmax25
    double lmr25top = 0.0;
    if (veg_lib->Ctype == 0) {
        lmr25top = vcmax25 * 0.015; // C3_veg = 0.015;
    }
    else {
        lmr25top = vcmax25 * 0.025; // C4_veg = 0.025;
    }
    // 遍历冠层
    double LAIcanopy = 0.0;
    double nscaler_sun = 0.0;
    double nscaler_sha = 0.0;
    for (i = 0; i < Ncanopy; i++) {
        if (i == 0) {
            LAIcanopy = 0.5 * NetLAI;
        }
        else {
            LAIcanopy = 0.5 * NetLAI;
        }
        if (Ncanopy == 1) {
            nscaler_sun = vcmaxcint_sun;
            nscaler_sha = vcmaxcint_sha;
        }
        else {
            nscaler_sun = exp(-f_N * LAIcanopy);
            nscaler_sha = exp(-f_N * LAIcanopy);
        }
        // 维护呼吸
        double lmr25_sun = lmr25top * nscaler_sun;
        double lmr25_sha = lmr25top * nscaler_sha;

        if (veg_lib->Ctype == 0) {
            lmr_sun[i] = lmr25_sun * ft(Tfoliage, lmrha) * fth(Tfoliage, lmrhd, lmrse, lmrc);
            lmr_sha[i] = lmr25_sha * ft(Tfoliage, lmrha) * fth(Tfoliage, lmrhd, lmrse, lmrc);
        } else {
            lmr_sun[i] = lmr25_sun * pow(2.0, (Tfoliage - (CONST_TKFRZ + 25.0)) / 10.0);
            lmr_sun[i] = lmr_sun[i] / (1.0 + exp(1.3 * (Tfoliage - (CONST_TKFRZ + 55.0))));
            lmr_sha[i] = lmr25_sha * pow(2.0, (Tfoliage - (CONST_TKFRZ + 25.0)) / 10.0);
            lmr_sha[i] = lmr_sha[i] / (1.0 + exp(1.3 * (Tfoliage - (CONST_TKFRZ + 55.0))));
        }
        // 低LAI时减少lmr
        lmr_sun[i] = lmr_sun[i] * min((0.2 * exp(3.218 * LAI_z[i])), 1.0);
        lmr_sha[i] = lmr_sha[i] * min((0.2 * exp(3.218 * LAI_z[i])), 1.0);
        // 光抑制
        if (light_inhibit && par_sun[i] > 0.0) {
            lmr_sun[i] = lmr_sun[i] * 0.67;
        }
        if (light_inhibit && par_sha[i] > 0.0) {
            lmr_sha[i] = lmr_sha[i] * 0.67;
        }
        // 白天计算光合作用
        if (par_sun[i] > 0.0) {
            // Vcmax, Jmax, TPU 的温度调整
            double vcmaxse = (668.39 - 1.07 * min(max((Tfoliage - CONST_TKFRZ), 11.0), 35.0)) * vcmaxse_sf;
            double jmaxse = (659.70 - 0.75 * min(max((Tfoliage - CONST_TKFRZ), 11.0), 35.0)) * jmaxse_sf;
            double tpuse = (668.39 - 1.07 * min(max((Tfoliage - CONST_TKFRZ), 11.0), 35.0)) * tpuse_sf;
            // 高温度抑制的缩放因子
            double lmrc = fth25(lmrhd, lmrse);
            double vcmaxc = fth25(vcmaxhd, vcmaxse);
            double jmaxc = fth25(jmaxhd, jmaxse);
            double tpuc = fth25(tpuhd, tpuse);
            
            vcmax_sun[i] = vcmax25 * nscaler_sun * ft(Tfoliage, vcmaxha) * fth(Tfoliage, vcmaxhd, vcmaxse, vcmaxc);
            jmax_sun[i] = jmax25 * nscaler_sun * ft(Tfoliage, jmaxha) * fth(Tfoliage, jmaxhd, jmaxse, jmaxc);
            tpu_sun[i] = tpu25 * nscaler_sun * ft(Tfoliage, tpuha) * fth(Tfoliage, tpuhd, tpuse, tpuc);
            
            vcmax_sha[i] = vcmax25 * nscaler_sha * ft(Tfoliage, vcmaxha) * fth(Tfoliage, vcmaxhd, vcmaxse, vcmaxc);
            jmax_sha[i] = jmax25 * nscaler_sha * ft(Tfoliage, jmaxha) * fth(Tfoliage, jmaxhd, jmaxse, jmaxc);
            tpu_sha[i] = tpu25 * nscaler_sha * ft(Tfoliage, tpuha) * fth(Tfoliage, tpuhd, tpuse, tpuc);
            
            if (veg_lib->Ctype == 1) {
                // C4植物的温度响应
                double temp_factor_sun = pow(2.0, (Tfoliage - (CONST_TKFRZ + 25.0)) / 10.0);
                temp_factor_sun = temp_factor_sun / (1.0 + exp(0.2 * ((CONST_TKFRZ + 15.0) - Tfoliage)));
                temp_factor_sun = temp_factor_sun / (1.0 + exp(0.3 * (Tfoliage - (CONST_TKFRZ + 40.0))));
                vcmax_sun[i] = vcmax25 * nscaler_sun * temp_factor_sun;
                
                double temp_factor_sha = pow(2.0, (Tfoliage - (CONST_TKFRZ + 25.0)) / 10.0);
                temp_factor_sha = temp_factor_sha / (1.0 + exp(0.2 * ((CONST_TKFRZ + 15.0) - Tfoliage)));
                temp_factor_sha = temp_factor_sha / (1.0 + exp(0.3 * (Tfoliage - (CONST_TKFRZ + 40.0))));
                vcmax_sha[i] = vcmax25 * nscaler_sha * temp_factor_sha;
            }
            
            kp_sun[i] = kp25 * nscaler_sun * pow(2.0, (Tfoliage - (CONST_TKFRZ + 25.0)) / 10.0);
            kp_sha[i] = kp25 * nscaler_sha * pow(2.0, (Tfoliage - (CONST_TKFRZ + 25.0)) / 10.0);

            if (energy->aPAR_sunlit <= 0.0) {
                double vegwp_sun = 1.0;
            }
            
            // 计算电子传输速率
            double qabs_sun = 0.5 * (1.0 - fnps) * par_sun[i] * 4.6;
            double qabs_sha = 0.5 * (1.0 - fnps) * par_sha[i] * 4.6;
            
            // 求解二次方程获取je
            double theta_psii = 0.7;  // 经验曲率参数
            double a = theta_psii;
            double b_sun = -(qabs_sun + jmax_sun[i]);
            double c_sun = qabs_sun * jmax_sun[i];
            double je_sun = solve_quadratic_min(a, b_sun, c_sun);
            
            double b_sha = -(qabs_sha + jmax_sha[i]);
            double c_sha = qabs_sha * jmax_sha[i];
            double je_sha = solve_quadratic_min(a, b_sha, c_sha);
            
            // 调用混合求解器计算Ci和gs
            double ceair = min(vp_over, esat_T);
            double rh_can;
            
            // 阳叶
            double qabs = 0.5 * (1.0 - fnps) * par_sun[i] * 4.6;
            double a = theta_psii;
            double b = -(qabs + jmax_sun[i]);
            double c_quad = qabs * jmax_sun[i];
            double r1, r2;
            solve_quadratic(a, b, c_quad, &r1, &r2);
            double je_sun = fmin(r1, r2);
            
            // 阴叶
            qabs = 0.5 * (1.0 - fnps) * par_sha[i] * 4.6;
            b = -(qabs + jmax_sha[i]);
            c_quad = qabs * jmax_sha[i];
            solve_quadratic(a, b, c_quad, &r1, &r2);
            double je_sha = fmin(r1, r2);
            
            // Ci初始猜测
            if (veg_lib->Ctype == 0) {
                ci_z_sun[i] = 0.7 * atmosCO2;
                ci_z_sha[i] = 0.7 * atmosCO2;
            } else {
                ci_z_sun[i] = 0.4 * atmosCO2;
                ci_z_sha[i] = 0.4 * atmosCO2;
            }
            
            // 寻找Ci和气孔导度
            int iter1, iter2;
            hybrid_PHS();
            
            double gsminsun, gsminsha, gs_slope_sun, gs_slope_sha;
            gsminsun = medlynintercept[i];
            gsminsha = medlynintercept[i];
            gs_slope_sun = medlynslope[i];
            gs_slope_sha = medlynslope[i];
            
            // 检查an < 0的情况
            if (an_sun[i] < 0.0) {
                gs_mol_sun[i] = fmax(bsun * gsminsun, 1.0);
            }
            if (an_sha[i] < 0.0) {
                gs_mol_sha[i] = fmax(bsha * gsminsha, 1.0);
            }
            
            // 午间气孔导度记录 (11AM-1PM)
            if (is_near_local_noon(londeg[g], 3600)) {
                gs_mol_sun_ln[i] = gs_mol_sun[i];
                gs_mol_sha_ln[i] = gs_mol_sha[i];
            } else {
                gs_mol_sun_ln[i] = SPVAL;
                gs_mol_sha_ln[i] = SPVAL;
            }
            
            // 最终计算cs和ci
            double cs_sun = atmosCO2 - 1.4 / RS_tmp * an_sun[i] * pressure;
            cs_sun = max(cs_sun, MAX_CS);
            ci_z_sun[i] = atmosCO2 - an_sun[i] * pressure *
                              (1.4 * gs_mol_sun[i] + 1.6 * RS_tmp) /
                              (RS_tmp * gs_mol_sun[i]);
            ci_z_sun[i] = max(ci_z_sun[i], param.TOL_A);
            
            double cs_sha = atmosCO2 - 1.4 / RS_tmp * an_sha[i] * pressure;
            cs_sha = max(cs_sha, MAX_CS);
            ci_z_sha[i] = atmosCO2 - an_sha[i] * pressure *
                              (1.4 * gs_mol_sha[i] + 1.6 * RS_tmp) /
                              (RS_tmp * gs_mol_sha[i]);
            ci_z_sha[i] = max(ci_z_sha[i], param.TOL_A);
            
            // 将gs_mol转换为gs (m/s)，再转换为rs (s/m)
            double gs = gs_mol_sun[i] / cf;
            rs_z_sun[i] = fmin(1.0 / gs, rsmax0);
            rs_z_sun[i] = rs_z_sun[p][iv] / o3coefg_sun;
            
            gs = gs_mol_sha[i] / cf;
            rs_z_sha[i] = fmin(1.0 / gs, rsmax0);
            rs_z_sha[i] = rs_z_sha[i] / o3coefg_sha;
            
            // 光合作用 - 保存限制性光合作用
            psn_z_sun[i] = ag_sun[i] * o3coefv_sun;
            
            
            // 检查迭代解的正确性
            if (gs_mol_sun[i] < 0.0 || gs_mol_sha[i] < 0.0) {
                printf("Negative stomatal conductance: p=%d, iv=%d, gs_mol_sun=%f, gs_mol_sha=%f\n",
                       p, iv, gs_mol_sun[p][iv], gs_mol_sha[p][iv]);
                exit(1);
            }
            
            // Ball-Berry模型校验
            double hs = (gb_mol * ceair + gs_mol_sun[i] * esat_T) /
                       ((gb_mol + gs_mol_sun[i]) * esat_T);
            rh_leaf_sun = hs;
            
            double gs_mol_err = gs_slope_sun * fmax(an_sun[i], 0.0) * hs / cs_sun * pressure +
                               fmax(bsun * gsminsun, 1.0);
            
            hs = (gb_mol * ceair + gs_mol_sha[i] * esat_T) /
                ((gb_mol + gs_mol_sha[i]) * esat_T);
            rh_leaf_sha = hs;
            
            gs_mol_err = gs_slope_sha * fmax(an_sha[i], 0.0) * hs / cs_sha * pressure +
                        max(bsha * gsminsha, 1.0);
        } else {
            // 夜间设置
            an_sun[i] = -bsun * lmr_sun[i];
            an_sha[i] = -bsha * lmr_sha[i];
            gs_mol_sun[i] = max(bsun * bbb, 1.0);
            gs_mol_sha[i] = max(bsha * bbb, 1.0);
        }
    }

    return (0);
}

/******************************************************************************
 * @brief    Calculate temperature scaling factor.
 *****************************************************************************/
double ft(double leaf_T, double ha)
{
    return exp(ha / (CONST_RGAS * (CONST_TKFRZ + 25.0)) * 
                    (1.0 - (CONST_TKFRZ + 25.0) / leaf_T));
}

/******************************************************************************
 * @brief    Calculate temperature scaling factor.
 *****************************************************************************/
double fth(double leaf_T, double hd, double se, double scalef)
{
    return scalef / (1.0 + exp((-hd + se * leaf_T) / (CONST_RGAS * leaf_T)));
}

/******************************************************************************
 * @brief    Calculate temperature scaling factor.
 *****************************************************************************/
double fth25(double hd, double se)
{
    return 1.0 + exp((-hd + se * (CONST_TKFRZ + 25.0)) / 
                    (CONST_RGAS * (CONST_TKFRZ + 25.0)));
}

/******************************************************************************
 * @brief    Calculate temperature scaling factor.
 *****************************************************************************/
double plc(double matric, double matric_50)
{
    double plc_val = 1.0;
    
    plc_val = pow(2.0, -pow(matric / matric_50, 4.0));
    if (plc_val < 0.005) {
        plc_val = 0.0;
    }
    return (plc_val);
}
