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
PhotoHydroStress(double             thm,
                 double             esat_T,
                 double             qsat_T,
                 double             vp_over,
                 double             Qair_over,
                 double             pressure,
                 double             air_density,
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
    double max_cs = 1.e-6; // Max CO2 partial pressure at leaf surface (Pa) for PHS
    double aPAR_sunlit = energy->aPAR_sunlit;
    double aPAR_shade = energy->aPAR_shade;
    double fcanopy = veg_var->fcanopy;
    double Tfoliage = energy->Tfoliage;
    double *dz_soil = soil_con->dz_soil;
    double *root = cell->root;
    double *LAI_z = veg_var->LAI_z;
    double *hksr_int = cell->hksr_int;
    double vegwp[4]; // 水势：0-阳叶，1-阴叶，2-木质部，3-根部
    // initialize variables
    Nsoil = cell->Nsoil;
    Ncanopy = cell->Ncanopy;
    double f_N = 0.0;
    double g_min = 2.0e3; // Ball-Berry minimum leaf conductance (mol/m^2/s)
    double CF = pressure / (CONST_RGAS * thm) * 1.0e6;
    double RS_mol = 1.0 / cell->Ra_leaf * CF;
    double tmp_photoleaf = 0.;
    double atmosO2 = pressure * 0.209;
    double atmosCO2 = pressure * 395.0e-06;
    double carbonylmax = 0.0;
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
        // fs: 由于根水势降低（更负）导致的导水率减少因子
        double fs = plc(matric[i], matric_50);
        // 根导水率：单位面积单位长度的导水率 (m/s)
        double hk_root = (fs * rai * krmax) / (croot_length + Zsum_soil[i]);
        hk_soil = max(hk_soil, 1.e-16);
        hk_root = max(hk_root, 1.e-16);
        // 计算土壤和根的总阻力，然后取倒数得到总导水率
        double Ra_total = 1.0 / hk_soil + 1.0 / hk_root;
        // 导水率是阻力的倒数
        // 对于表层土壤，明确设置导水率为0
        if (rai * root[i] > 0.0 && i > 0) {
            hksr_int[i] = 1.0 / Ra_total;
        } 
        else {
            hksr_int[i] = 0.0;
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
    double jmax25 = ((2.59 - 0.035 * min(max((Tfoliage - CONST_TKFRZ), 
                        11.0), 35.0)) * vcmax25) * jmaxse_sf;
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
    if (veg_lib->Ctype == 0) { // C3植物
        lmr25top = vcmax25 * 0.015; // C3_veg = 0.015;
    }
    else {  // C4植物
        lmr25top = vcmax25 * 0.025; // C4_veg = 0.025;
    }
    // 遍历冠层
    bool light_inhibit = false; // 是否存在光抑制
    double rsmax0 = 2.0e4; // Maximum stomatal resistance (s/m)
    double bsun, bsha;
    double gsun, gsha;
    double an_sun, an_sha;
    double ci_sun, ci_sha;
    double gs_sun, gs_sha;
    double ag_sun, ag_sha;
    double ac_sun, ac_sha;
    double aj_sun, aj_sha;
    double ap_sun, ap_sha;
    double lmr_sun, lmr_sha;
    double vcmax_sun, vcmax_sha;
    double tpu_sun, tpu_sha;
    double kp_sun, kp_sha;
    double jmax_sun, jmax_sha;
    double hs = 0.0;
    double gs = 0.0;
    double cs_sun = 0.0;
    double cs_sha = 0.0;
    double gs_mol_sun = 0.0;
    double gs_mol_sha = 0.0;
    double rs_sun = 0.0;
    double rs_sha = 0.0;
    double psn_sun[MAX_CANOPYS];
    double psn_sha[MAX_CANOPYS];
    double psn_wc_sun[MAX_CANOPYS];
    double psn_wj_sun[MAX_CANOPYS];
    double psn_wp_sun[MAX_CANOPYS];
    double psn_wc_sha[MAX_CANOPYS];
    double psn_wj_sha[MAX_CANOPYS];
    double psn_wp_sha[MAX_CANOPYS];
    // 初始化光合作用和水分胁迫变量
    for (i = 0; i < Ncanopy; i++) {
        psn_sun[i] = 0.0;
        psn_sha[i] = 0.0;
        psn_wc_sun[i] = 0.0;
        psn_wj_sun[i] = 0.0;
        psn_wp_sun[i] = 0.0;
        psn_wc_sha[i] = 0.0;
        psn_wj_sha[i] = 0.0;
        psn_wp_sha[i] = 0.0;
    }
    double LAIcanopy = 0.0;
    double nscaler_sun = 0.0;
    double nscaler_sha = 0.0;
    double gs_mol_err = 0.0;
    for (i = 0; i < Ncanopy; i++) {
        if (i == 0) {
            LAIcanopy = 0.5 * NetLAI;
        }
        else {
            LAIcanopy += 0.5 * (LAI_z[i-1] + LAI_z[i]);
        }
        if (Ncanopy == 1) {
            nscaler_sun = veg_var->ksun_vcmax;
            nscaler_sha = veg_var->ksha_vcmax;
        }
        else {
            nscaler_sun = exp(-f_N * LAIcanopy);
            nscaler_sha = exp(-f_N * LAIcanopy);
        }
        // 维护呼吸
        double lmr25_sun = lmr25top * nscaler_sun;
        double lmr25_sha = lmr25top * nscaler_sha;

        if (veg_lib->Ctype == 0) {
            lmr_sun = lmr25_sun * ft(Tfoliage, lmrha) * fth(Tfoliage, lmrhd, lmrse, lmrc);
            lmr_sha = lmr25_sha * ft(Tfoliage, lmrha) * fth(Tfoliage, lmrhd, lmrse, lmrc);
        } 
        else {
            lmr_sun = lmr25_sun * pow(2.0, (Tfoliage - (CONST_TKFRZ + 25.0)) / 10.0);
            lmr_sun = lmr_sun / (1.0 + exp(1.3 * (Tfoliage - (CONST_TKFRZ + 55.0))));
            lmr_sha = lmr25_sha * pow(2.0, (Tfoliage - (CONST_TKFRZ + 25.0)) / 10.0);
            lmr_sha = lmr_sha / (1.0 + exp(1.3 * (Tfoliage - (CONST_TKFRZ + 55.0))));
        }
        // 低LAI时减少lmr
        lmr_sun = lmr_sun * min((0.2 * exp(3.218 * LAI_z[i])), 1.0);
        lmr_sha = lmr_sha * min((0.2 * exp(3.218 * LAI_z[i])), 1.0);
        // 光抑制
        if (light_inhibit && energy->aPAR_sunlit > 0.0) {
            lmr_sun = lmr_sun * 0.67;
        }
        if (light_inhibit && energy->aPAR_shade > 0.0) {
            lmr_sha = lmr_sha * 0.67;
        }
        // 白天计算光合作用
        if (veg_var->leaf_sun > 0.0) {
            // Vcmax, Jmax, TPU 的温度调整
            double vcmaxse = (668.39 - 1.07 * min(max((Tfoliage - CONST_TKFRZ), 11.0), 35.0)) * vcmaxse_sf;
            double jmaxse = (659.70 - 0.75 * min(max((Tfoliage - CONST_TKFRZ), 11.0), 35.0)) * jmaxse_sf;
            double tpuse = (668.39 - 1.07 * min(max((Tfoliage - CONST_TKFRZ), 11.0), 35.0)) * tpuse_sf;
            // 高温度抑制的缩放因子
            double lmrc = fth25(lmrhd, lmrse);
            double vcmaxc = fth25(vcmaxhd, vcmaxse);
            double jmaxc = fth25(jmaxhd, jmaxse);
            double tpuc = fth25(tpuhd, tpuse);
            
            vcmax_sun = vcmax25 * nscaler_sun * ft(Tfoliage, vcmaxha) * fth(Tfoliage, vcmaxhd, vcmaxse, vcmaxc);
            jmax_sun = jmax25 * nscaler_sun * ft(Tfoliage, jmaxha) * fth(Tfoliage, jmaxhd, jmaxse, jmaxc);
            tpu_sun = tpu25 * nscaler_sun * ft(Tfoliage, tpuha) * fth(Tfoliage, tpuhd, tpuse, tpuc);
            
            vcmax_sha = vcmax25 * nscaler_sha * ft(Tfoliage, vcmaxha) * fth(Tfoliage, vcmaxhd, vcmaxse, vcmaxc);
            jmax_sha = jmax25 * nscaler_sha * ft(Tfoliage, jmaxha) * fth(Tfoliage, jmaxhd, jmaxse, jmaxc);
            tpu_sha = tpu25 * nscaler_sha * ft(Tfoliage, tpuha) * fth(Tfoliage, tpuhd, tpuse, tpuc);
            
            if (veg_lib->Ctype == 1) {
                // C4植物的温度响应
                double temp_factor_sun = pow(2.0, (Tfoliage - (CONST_TKFRZ + 25.0)) / 10.0);
                temp_factor_sun = temp_factor_sun / (1.0 + exp(0.2 * ((CONST_TKFRZ + 15.0) - Tfoliage)));
                temp_factor_sun = temp_factor_sun / (1.0 + exp(0.3 * (Tfoliage - (CONST_TKFRZ + 40.0))));
                vcmax_sun = vcmax25 * nscaler_sun * temp_factor_sun;
                
                double temp_factor_sha = pow(2.0, (Tfoliage - (CONST_TKFRZ + 25.0)) / 10.0);
                temp_factor_sha = temp_factor_sha / (1.0 + exp(0.2 * ((CONST_TKFRZ + 15.0) - Tfoliage)));
                temp_factor_sha = temp_factor_sha / (1.0 + exp(0.3 * (Tfoliage - (CONST_TKFRZ + 40.0))));
                vcmax_sha = vcmax25 * nscaler_sha * temp_factor_sha;
            }
            
            kp_sun = kp25 * nscaler_sun * pow(2.0, (Tfoliage - (CONST_TKFRZ + 25.0)) / 10.0);
            kp_sha = kp25 * nscaler_sha * pow(2.0, (Tfoliage - (CONST_TKFRZ + 25.0)) / 10.0);
            
            // 调用混合求解器计算Ci和gs
            double ceair = min(vp_over, esat_T);
            double rh_can = ceair / esat_T;
            
            // 阳叶
            double theta_psii = 0.7;  // 经验曲率参数
            double r1, r2;
            double aquad, bquad, cquad;
            double qabs = 0.5 * (1.0 - fnps) * energy->aPAR_sunlit * 4.6;
            aquad = theta_psii;
            bquad = -(qabs + jmax_sun);
            cquad = qabs * jmax_sun;
            solve_quadratic(aquad, bquad, cquad, &r1, &r2);
            double je_sun = min(r1, r2);
            
            // 阴叶
            qabs = 0.5 * (1.0 - fnps) * energy->aPAR_shade * 4.6;
            bquad = -(qabs + jmax_sha);
            cquad = qabs * jmax_sha;
            solve_quadratic(aquad, bquad, cquad, &r1, &r2);
            double je_sha = min(r1, r2);
            // 初始猜测Ci和气孔导度，后续迭代更新
            double ci_sun, ci_sha;
            if (veg_lib->Ctype == 0) {
                ci_sun = 0.7 * atmosCO2; // 初始猜测Ci为大气CO2的70%
                ci_sha = 0.7 * atmosCO2;
            }
            else {
                ci_sun = 0.4 * atmosCO2; // C4植物初始猜测Ci为大气CO2的40%
                ci_sha = 0.4 * atmosCO2;
            }
            // 计算ci和气孔导度的混合求解器
            hybrid_PHS(&ci_sun, &ci_sha,
                        vegwp, &bsun, &bsha,
                        je_sun, je_sha,
                        thm, RS_mol,
                        qsat_T, Qair_over,
                        pressure, air_density,
                        atmosCO2, atmosO2,
                        lmr_sun, lmr_sha,
                        energy->aPAR_sunlit, 
                        energy->aPAR_shade,
                        rh_can,
                        &kp_sun, &kp_sha,
                        &gsun, &gsha,
                        energy, cell, 
                        soil_con, veg_var, veg_lib);
            if (an_sun < 0.0) {
                gs_mol_sun = max(bsun * veg_lib->medlynint, 1.0);
            }
            if (an_sha < 0.0) {
                gs_mol_sha = max(bsha * veg_lib->medlynint, 1.0);
            }
            // ========== 计算叶面CO₂分压cs和最终ci ==========
            // 晴叶
            cs_sun = atmosCO2 - 1.4 / RS_mol * an_sun * pressure;
            cs_sun = max(cs_sun, max_cs);
            
            ci_sun = atmosCO2 - an_sun * pressure *
                (1.4 * gs_mol_sun + 1.6 * RS_mol) /
                (RS_mol * gs_mol_sun);
            ci_sun = max(ci_sun, 1.0e-6);
            
            // 阴叶
            cs_sha = atmosCO2 - 1.4 / RS_mol * an_sha * pressure;
            cs_sha = max(cs_sha, max_cs);
            
            ci_sha = atmosCO2 - an_sha * pressure *
                (1.4 * gs_mol_sha + 1.6 * RS_mol) /
                (RS_mol * gs_mol_sha);
            ci_sha = max(ci_sha, 1.0e-6);
            
            // ========== 转换为气孔阻力 ==========
            gs = gs_mol_sun / CF;
            rs_sun = min(1.0 / gs, rsmax0);
            
            gs = gs_mol_sha / CF;
            rs_sha = min(1.0 / gs, rsmax0);
            
            // ========== 记录光合速率和限制因子 ==========
            // 阳叶
            psn_sun[i] = ag_sun;
            psn_wc_sun[i] = 0.0;
            psn_wj_sun[i] = 0.0;
            psn_wp_sun[i] = 0.0;
            // 根据限制因子分配光合速率
            if (ac_sun <= aj_sun && ac_sun <= ap_sun) {
                psn_wc_sun[i] = psn_sun[i];
            } 
            else if (aj_sun < ac_sun && aj_sun <= ap_sun) {
                psn_wj_sun[i] = psn_sun[i];
            } 
            else if (ap_sun < ac_sun && ap_sun < aj_sun) {
                psn_wp_sun[i] = psn_sun[i];
            } 
            else if (ap_sun < ac_sun && ap_sun < aj_sun) {
                psn_wp_sun[i] = psn_sun[i];
            }
            // 阴叶
            psn_sha[i] = ag_sha;
            psn_wc_sha[i] = 0.0;
            psn_wj_sha[i] = 0.0;
            psn_wp_sha[i] = 0.0;
            
            if (ac_sha <= aj_sha && ac_sha <= ap_sha) {
                psn_wc_sha[i] = psn_sha[i];
            } 
            else if (aj_sha < ac_sha && aj_sha <= ap_sha) {
                psn_wj_sha[i] = psn_sha[i];
            } 
            else if (ap_sha < ac_sha && ap_sha < aj_sha) {
                psn_wp_sha[i] = psn_sha[i];
            }
            
            // ========== Ball-Berry 一致性校验 ==========
            // 晴叶
            hs = (RS_mol * ceair + gs_mol_sun * Qair_over) / ((RS_mol + gs_mol_sun) * Qair_over);
            
            gs_mol_err = veg_lib->medlynslope * max(an_sun, 0.0) * hs / cs_sun *
                            pressure + max(bsun * veg_lib->medlynint, 1.0);
            
            // 阴叶
            hs = (RS_mol * ceair + gs_mol_sha * Qair_over) / ((RS_mol + gs_mol_sha) * Qair_over);
            
            gs_mol_err = veg_lib->medlynslope * max(an_sha, 0.0) * hs / cs_sha *
                            pressure + max(bsha * veg_lib->medlynint, 1.0);
                
        }
        else {
            // 夜间或无光照条件下，光合作用为0，维护呼吸仍然存在
            an_sun = 0.0;
            an_sha = 0.0;
            ci_sun = 0.0;
            ci_sha = 0.0;
            gs_mol_sun = 0.0;
            gs_mol_sha = 0.0;
            rs_sun = rsmax0;
            rs_sha = rsmax0;
            psn_sun[i] = 0.0;
            psn_sha[i] = 0.0;
            psn_wc_sun[i] = 0.0;
            psn_wj_sun[i] = 0.0;
            psn_wp_sun[i] = 0.0;
            psn_wc_sha[i] = 0.0;
            psn_wj_sha[i] = 0.0;
            psn_wp_sha[i] = 0.0;
            // 调用calc_stress函数计算水分胁迫因子
            calc_stress(&bsun, &bsha, vegwp, thm, RS_mol, 
                        qsat_T, Qair_over,
                        pressure, air_density,
                        energy, cell, 
                        soil_con, veg_var, veg_lib);

        }
    }
    // 计算冠层阻力
    double laican_sun = 0.0;
    double gscan_sun = 0.0;
    double RS_sunlit = 0.0;
    double RS_shade = 0.0;
    for (i = 0; i < Ncanopy; i++) {
        gscan_sun += LAI_z[i] / (rs_sun + cell->Ra_leaf);
        laican_sun += LAI_z[i];
    }
    if (laican_sun > 0.0) {
        RS_sunlit = laican_sun / gscan_sun - cell->Ra_leaf;
    }
    else {
        RS_sunlit = 0.0;
    }
    double laican_sha = 0.0;
    double gscan_sha = 0.0;
    for (i = 0; i < Ncanopy; i++) {
        gscan_sha += LAI_z[i] / (rs_sha + cell->Ra_leaf);
        laican_sha += LAI_z[i];
    }
    if (laican_sha > 0.0) {
        RS_shade = laican_sha / gscan_sha - cell->Ra_leaf;
    }
    else {
        RS_shade = 0.0;
    }
    if (laican_sun + laican_sha > 0.0) {
        cell->transp_fact = bsun * (laican_sun / (laican_sun + laican_sha)) +
                            bsha * (laican_sha / (laican_sun + laican_sha));
    }
    else {
        cell->transp_fact = bsun;
    }
    veg_var->RS_sunlit = RS_sunlit;
    veg_var->RS_shade = RS_shade;

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
/******************************************************************************
 * @brief    Calculate temperature scaling factor.
 *****************************************************************************/
double d1plc(double matric, double matric_50, double ck)
{
    double ratio_pow = pow(matric / matric_50, ck);
    double exp_term = pow(2.0, -ratio_pow);
    double dplc_dmatric = -ck * log(2.0) * exp_term * ratio_pow / matric;
    return (dplc_dmatric);
}
