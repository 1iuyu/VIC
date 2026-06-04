/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine computes canopy stomatal resistance and foliage photosynthesis
 * based on plant hydraulic stress scheme.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    Calculate canopy photosynthesis and stomatal resistance.
 *****************************************************************************/
int
photosynth_hydrostress(double            thm,
                       double            daylen,
                       double            esat_T,
                       double            qsat_T,
                       double            vp_over,
                       double            Qair_over,
                       double            pressure,
                       double            air_density,
                       double            Tfoliage,
                       cell_data_struct *cell,
                       soil_con_struct  *soil_con,
                       veg_var_struct   *veg_var,
                       veg_lib_struct   *veg_lib)
{   
    extern parameters_struct param;
    size_t i, Ncanopy;
    double aPAR_sun = veg_var->aPAR_sun;
    double aPAR_sha = veg_var->aPAR_sha;
    double *dz_soil = soil_con->dz_soil;
    double *root = cell->root;
    double *LAI_z = veg_var->LAI_z;
    double *hksr_int = cell->hksr_int; // soil-root interface conductance (mm/s)
    double *mat_VEG = veg_var->mat_VEG; // 水势：0-阳叶，1-阴叶，2-木质部，3-根部
    // initialize variables
    Ncanopy = cell->Ncanopy;
    double f_N = 0.0;
    double CF = pressure / (CONST_RGAS * thm) * 1.0e6;
    double gb_mol = 1.0 / cell->Ra_leaf * CF;
    double atmosO2 = pressure * param.PHOTO_OX;
    double atmosCO2 = pressure * param.PHOTO_CX;
    double rho_root = 0.0;
    double froot_leaf = veg_lib->froot_leaf;
    double NetLAI = veg_var->NetLAI;
    double NetSAI = veg_var->NetSAI;
    double *matric = cell->matric;
    double *zc_soil = soil_con->zc_soil;
    double *conductivity = cell->conductivity;
    double lmrc = fth25(param.PHOTO_LMRHD, param.PHOTO_LMRSE);
    double max_daylen = cell->max_daylen;
    double daylen_fact = min(1.0, max(0.01, (daylen * daylen) / (max_daylen * max_daylen)));
    // 计算常数：根横截面积 (m2)
    double xsec_root = CONST_PI * pow(param.SOIL_RROOT, 2.0);
    
    for (i = 0; i < cell->Nroot; i++) {
        // 根生物量密度：g biomass/m3 soil
        rho_root = 2.0 * param.PHOTO_CROOT * root[i] / dz_soil[i];
        rho_root = max(2.0, rho_root);
        double rlen_dens = rho_root / (xsec_root * param.SOIL_RHOROOT);
        // 根面积指数 (RAI)
        double rai = (NetLAI + NetSAI) * froot_leaf * root[i];
        // 粗根长度 - 使用指定的侧向长度
        double croot_length = 0.25; // specified lateral coarse root length [m]
        double r_soil = sqrt(1.0 / (CONST_PI * rlen_dens));
        double hk_soil = conductivity[i] / r_soil;
        double fs = plc(matric[i], veg_lib->matric50);
        double hk_root = (fs * rai * veg_lib->kroot_max) / (croot_length + zc_soil[i]);
        hk_soil = max(hk_soil, param.PHOTO_MINCONDUCT);
        hk_root = max(hk_root, param.PHOTO_MINCONDUCT);
        // 计算土壤和根的总阻力，然后取倒数得到总导水率
        double Ra_total = 1.0 / hk_soil + 1.0 / hk_root;
        // 对于表层土壤，明确设置导水率为0
        if (rai * root[i] > 0.0 && i > 0) {
            hksr_int[i] = 1.0 / Ra_total;
        } 
        else {
            hksr_int[i] = 0.0;
        }
    }
    // Miscellaneous parameters
    double sco  = 0.5 * param.PHOTO_OX / param.PHOTO_CP;
    double cp25 = 0.5 * atmosO2 / sco;
    double KC = param.PHOTO_KC * pressure * ft(Tfoliage, param.PHOTO_EC);
    double KO = param.PHOTO_KO * pressure * ft(Tfoliage, param.PHOTO_EO);
    double CP = cp25 * ft(Tfoliage, param.PHOTO_EP);
    // 叶片氮浓度
    double lnc = min(1.0 / (veg_lib->leaf_CN * veg_lib->SLA_top), 10.0);
    double vcmax25 = lnc * veg_lib->fN_rub * param.PHOTO_FNR * param.PHOTO_SACT * daylen_fact;
    double jmax25 = ((2.59 - 0.035 * min(max((Tfoliage - CONST_TKFRZ), 
                        11.0), 35.0)) * vcmax25);
    double tpu25 = param.PHOTO_FTPU * vcmax25;
    double kp25 = param.PHOTO_FKP * vcmax25;
    // 氮衰减系数
    if (daylen_fact == 0.0) {
        f_N = 0.0;
    }
    else {
        f_N = exp(0.00963 * vcmax25 / daylen_fact - 2.43);
    }
    // Leaf maintenance respiration in proportion to vcmax25
    double lmr25top = 0.0;
    if (veg_lib->Ctype == PHOTO_C3) {
        lmr25top = vcmax25 * param.PHOTO_LRESC3; // C3_veg = 0.015;
    }
    else {
        lmr25top = vcmax25 * param.PHOTO_LRESC4; // C4_veg = 0.025;
    }
    // 遍历冠层
    bool light_inhibit = false; // 是否存在光抑制
    double bsun, bsha;
    double ci_sun, ci_sha;
    double lmr_sun, lmr_sha;
    double vcmax_sun, vcmax_sha;
    double tpu_sun, tpu_sha;
    double kp_sun, kp_sha;
    double jmax_sun, jmax_sha;
    double hs, gs;
    double cs_sun = 0.0;
    double cs_sha = 0.0;
    double gs_mol_sun = 0.0;
    double gs_mol_sha = 0.0;
    double Ra_sun[MAX_CANOPYS];
    double Ra_sha[MAX_CANOPYS];
    // 初始化光合作用和水分胁迫变量
    for (i = 0; i < MAX_CANOPYS; i++) {
        Ra_sun[i] = 0.0;
        Ra_sha[i] = 0.0;
    }
    double LAIcanopy = 0.0;
    double nscaler_sun = 0.0;
    double nscaler_sha = 0.0;
    for (i = 0; i < Ncanopy; i++) {
        if (i == 0) {
            LAIcanopy = 0.5 * LAI_z[i];
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

        if (veg_lib->Ctype == PHOTO_C3) {
            lmr_sun = lmr25_sun * ft(Tfoliage, param.PHOTO_EL) * 
                        fth(Tfoliage, param.PHOTO_LMRHD, param.PHOTO_LMRSE, lmrc);
            lmr_sha = lmr25_sha * ft(Tfoliage, param.PHOTO_EL) * 
                        fth(Tfoliage, param.PHOTO_LMRHD, param.PHOTO_LMRSE, lmrc);
        } 
        else {
            lmr_sun = lmr25_sun * pow(2.0, (Tfoliage - 298.15) / 10.0); // Q10 = 2.0, 298.15K = 25C
            lmr_sun = lmr_sun / (1.0 + exp(1.3 * (Tfoliage - 328.15))); // 高温抑制，328.15K = 55C
            lmr_sha = lmr25_sha * pow(2.0, (Tfoliage - 298.15) / 10.0);
            lmr_sha = lmr_sha / (1.0 + exp(1.3 * (Tfoliage - 328.15)));
        }
        // 低LAI时减少lmr
        lmr_sun = lmr_sun * min((0.2 * exp(3.218 * LAI_z[i])), 1.0);
        lmr_sha = lmr_sha * min((0.2 * exp(3.218 * LAI_z[i])), 1.0);
        // 光抑制
        if (light_inhibit && aPAR_sun > 0.0) {
            lmr_sun = lmr_sun * 0.67;
        }
        if (light_inhibit && aPAR_sha > 0.0) {
            lmr_sha = lmr_sha * 0.67;
        }
        // 白天计算光合作用
        if (veg_var->aPAR_sun > 0.0) {
            // Vcmax, Jmax, TPU 的温度调整
            double vcmaxse = (668.39 - 1.07 * min(max((Tfoliage - CONST_TKFRZ), 11.0), 35.0));
            double jmaxse = (659.70 - 0.75 * min(max((Tfoliage - CONST_TKFRZ), 11.0), 35.0));
            double tpuse = (668.39 - 1.07 * min(max((Tfoliage - CONST_TKFRZ), 11.0), 35.0));
            // 高温度抑制的缩放因子
            double vcmaxc = fth25(param.PHOTO_DV, vcmaxse);
            double jmaxc = fth25(param.PHOTO_DJ, jmaxse);
            double tpuc = fth25(param.PHOTO_DT, tpuse);
            
            vcmax_sun = vcmax25 * nscaler_sun * ft(Tfoliage, param.PHOTO_EV) * fth(Tfoliage, param.PHOTO_DV, vcmaxse, vcmaxc);
            jmax_sun = jmax25 * nscaler_sun * ft(Tfoliage, param.PHOTO_EJ) * fth(Tfoliage, param.PHOTO_DJ, jmaxse, jmaxc);
            tpu_sun = tpu25 * nscaler_sun * ft(Tfoliage, param.PHOTO_ET) * fth(Tfoliage, param.PHOTO_DT, tpuse, tpuc);
            
            vcmax_sha = vcmax25 * nscaler_sha * ft(Tfoliage, param.PHOTO_EV) * fth(Tfoliage, param.PHOTO_DV, vcmaxse, vcmaxc);
            jmax_sha = jmax25 * nscaler_sha * ft(Tfoliage, param.PHOTO_EJ) * fth(Tfoliage, param.PHOTO_DJ, jmaxse, jmaxc);
            tpu_sha = tpu25 * nscaler_sha * ft(Tfoliage, param.PHOTO_ET) * fth(Tfoliage, param.PHOTO_DT, tpuse, tpuc);
            
            if (veg_lib->Ctype == PHOTO_C4) {
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
            double rh_over = max((esat_T - ceair), 50.0) * 0.001;
            
            // 阳叶
            double theta_psii = 0.7;  // 经验曲率参数
            double r1, r2;
            double aquad, bquad, cquad;
            double qabs = 0.5 * (1.0 - param.PHOTO_FNPS) * aPAR_sun * 4.6;
            aquad = theta_psii;
            bquad = -(qabs + jmax_sun);
            cquad = qabs * jmax_sun;
            solve_quadratic(aquad, bquad, cquad, &r1, &r2);
            double je_sun = min(r1, r2);
            
            // 阴叶
            qabs = 0.5 * (1.0 - param.PHOTO_FNPS) * aPAR_sha * 4.6;
            bquad = -(qabs + jmax_sha);
            cquad = qabs * jmax_sha;
            solve_quadratic(aquad, bquad, cquad, &r1, &r2);
            double je_sha = min(r1, r2);
            // 初始猜测Ci和气孔导度，后续迭代更新
            if (veg_lib->Ctype == PHOTO_C3) {
                ci_sun = 0.7 * atmosCO2; // 初始猜测Ci为大气CO2的70%
                ci_sha = 0.7 * atmosCO2;
            }
            else {
                ci_sun = 0.4 * atmosCO2; // C4植物初始猜测Ci为大气CO2的40%
                ci_sha = 0.4 * atmosCO2;
            }
            // 计算ci和气孔导度
            hybrid_PHS(&ci_sun, &ci_sha,
                        mat_VEG, &bsun, &bsha,
                        je_sun, je_sha,
                        CP, KC, KO,
                        thm, gb_mol,
                        qsat_T, Qair_over,
                        pressure, air_density,
                        atmosCO2, atmosO2,
                        lmr_sun, lmr_sha,
                        rh_over,
                        vcmax_sun, vcmax_sha,
                        tpu_sun, tpu_sha,
                        kp_sun, kp_sha,
                        &gs_mol_sun, &gs_mol_sha,
                        cell, soil_con,
                        veg_var, veg_lib);

            if (veg_var->an_sun < 0.0) {
                gs_mol_sun = max(bsun * veg_lib->medlynint, 1.0);
            }
            if (veg_var->an_sha < 0.0) {
                gs_mol_sha = max(bsha * veg_lib->medlynint, 1.0);
            }
            // 计算叶面CO₂分压cs和最终ci
            // 阳叶
            cs_sun = atmosCO2 - 1.4 / gb_mol * veg_var->an_sun * pressure;
            cs_sun = max(cs_sun, param.PHOTO_MAXCS);
            
            ci_sun = atmosCO2 - veg_var->an_sun * pressure *
                (1.4 * gs_mol_sun + 1.6 * gb_mol) /
                (gb_mol * gs_mol_sun);
            ci_sun = max(ci_sun, 1.0e-6);
            
            // 阴叶
            cs_sha = atmosCO2 - 1.4 / gb_mol * veg_var->an_sha * pressure;
            cs_sha = max(cs_sha, param.PHOTO_MAXCS);
            
            ci_sha = atmosCO2 - veg_var->an_sha * pressure *
                (1.4 * gs_mol_sha + 1.6 * gb_mol) /
                (gb_mol * gs_mol_sha);
            ci_sha = max(ci_sha, 1.0e-6);
            
            // ========== 转换为气孔阻力 ==========
            gs = gs_mol_sun / CF;
            Ra_sun[i] = min(1.0 / gs, param.PHOTO_RSMAX);
            
            gs = gs_mol_sha / CF;
            Ra_sha[i] = min(1.0 / gs, param.PHOTO_RSMAX);
            
            // ========== Ball-Berry 一致性校验 ==========
            // 阳叶
            hs = (gb_mol * ceair + gs_mol_sun * Qair_over) / ((gb_mol + gs_mol_sun) * Qair_over);
            
            veg_var->PhotoError[0] = veg_lib->medlynslope * max(veg_var->an_sun, 0.0) * hs / cs_sun *
                            pressure + max(bsun * veg_lib->medlynint, 1.0);
            
            // 阴叶
            hs = (gb_mol * ceair + gs_mol_sha * Qair_over) / ((gb_mol + gs_mol_sha) * Qair_over);
            
            veg_var->PhotoError[1] = veg_lib->medlynslope * max(veg_var->an_sha, 0.0) * hs / cs_sha *
                            pressure + max(bsha * veg_lib->medlynint, 1.0);
                
        }
        else {
            // 夜间或无光照条件下，光合作用为0，维护呼吸仍然存在
            mat_VEG[0] = 1.0; // temporary signal for night time
            double gsminsun = veg_lib->medlynint;
            double gsminsha = veg_lib->medlynint;
            // 调用calc_stress函数计算水分胁迫因子
            calc_stress(&bsun, &bsha, 
                        mat_VEG, thm, gb_mol, 
                        qsat_T, Qair_over,
                        pressure, 
                        air_density,
                        gsminsun, gsminsha, 
                        cell, soil_con, 
                        veg_var, veg_lib);

            veg_var->ac_sun = 0.0;
            veg_var->ac_sha = 0.0;
            veg_var->ag_sun = 0.0;
            veg_var->ag_sha = 0.0;
            veg_var->aj_sun = 0.0;
            veg_var->aj_sha = 0.0;
            veg_var->an_sun = veg_var->ag_sun - bsun * lmr_sun;
            veg_var->an_sha = veg_var->ag_sha - bsha * lmr_sha;
            veg_var->ap_sun = 0.0;
            veg_var->ap_sha = 0.0;
            Ra_sun[i] = min(param.PHOTO_RSMAX, 1.0 / (max(bsun * gsminsun, 1.0)) * CF);
            Ra_sha[i] = min(param.PHOTO_RSMAX, 1.0 / (max(bsha * gsminsha, 1.0)) * CF);
        }
    }
    // 计算冠层阻力
    double laican_sun = 0.0;
    double gscan_sun = 0.0;
    double RS_sunlit = 0.0;
    double RS_shade = 0.0;
    for (i = 0; i < Ncanopy; i++) {
        gscan_sun += LAI_z[i] / (Ra_sun[i] + cell->Ra_leaf);
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
        gscan_sha += LAI_z[i] / (Ra_sha[i] + cell->Ra_leaf);
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
