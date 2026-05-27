/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate the transpiration stress using a plant hydraulics approach.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief  Calculate the transpiration stress using a plant hydraulics approach.
 *****************************************************************************/
void ci_func_PHS(bool              bflag,
                 double            cisun, 
                 double            cisha,
                 double           *fvalsun, 
                 double           *fvalsha,
                 double           *bsun, 
                 double           *bsha,
                 double           *gs_mol_sun, 
                 double           *gs_mol_sha,
                 double           *vegwp,
                 double            gsminsun,
                 double            gsminsha,
                 double            vcmax_sun, 
                 double            vcmax_sha,
                 double            tpu_sun, 
                 double            tpu_sha,
                 double            kp_sun, 
                 double            kp_sha,
                 double            CP, 
                 double            KC, 
                 double            KO,
                 double            thm, 
                 double            gb_mol,
                 double            qsat_T, 
                 double            Qair_over,
                 double            pressure, 
                 double            air_density,
                 double            jesun, 
                 double            jesha,
                 double            atmosCO2, 
                 double            atmosO2,
                 double            lmr_sun, 
                 double            lmr_sha,
                 double            rh_canopy,
                 cell_data_struct *cell,
                 soil_con_struct  *soil_con,
                 veg_var_struct   *veg_var,
                 veg_lib_struct   *veg_lib)
{
    double ai;                  // 中间协同限制光合速率
    double cs_sun, cs_sha;      // 叶面CO₂分压 (Pa)
    double aquad, bquad, cquad; // 二次方程系数
    double r1, r2;              // 二次方程根
    double term;                // Medlyn模型中间变量
    double an_sun, an_sha;      // 净光合速率
    double ac_sun, ac_sha;      // Rubisco限制光合
    double aj_sun, aj_sha;      // RuBP限制光合
    double ap_sun, ap_sha;      // 产物限制光合
    double ag_sun, ag_sha;      // 协同限制总光合
    double max_cs = 50000.0;    // 最大叶面CO₂分压 (Pa)
    double aPAR_sun = veg_var->aPAR_sun;
    double aPAR_sha = veg_var->aPAR_sha;
    double qe = 0.0; // Quantum efficiency (mol CO2 / mol photons) C3;
    if (veg_lib->Ctype == 1) {
        qe = 0.05;  // C4
    }
    double medlynslope = veg_lib->medlynslope;
    double medlynint = veg_lib->medlynint;
    double theta_cj = veg_lib->theta_cj; // AC-AJ耦合系数
    double theta_ip = 0.95; // AI-AP耦合系数
    
    // 如果需要重新计算水分胁迫因子
    if (bflag) {
        calc_stress(&bsun, &bsha, vegwp, thm, gb_mol, 
                    qsat_T, Qair_over,
                    pressure, air_density, 
                    gsminsun, gsminsha,
                    cell, soil_con, veg_var, veg_lib);
    }
    
    // ========== 计算三种限制下的光合速率 ==========
    if (veg_lib->Ctype == 0) {
        // C3植物
        // Rubisco限制
        ac_sun = (*bsun) * vcmax_sun * max(cisun - CP, 0.0) / (cisun + KC * (1.0 + atmosO2 / KO));
        ac_sha = (*bsha) * vcmax_sha * max(cisha - CP, 0.0) / (cisha + KC * (1.0 + atmosO2 / KO));
        
        // RuBP限制
        aj_sun = jesun * max(cisun - CP, 0.0) / (4.0 * cisun + 8.0 * CP);
        aj_sha = jesha * max(cisha - CP, 0.0) / (4.0 * cisha + 8.0 * CP);
        
        // 产物限制
        ap_sun = 3.0 * tpu_sun;
        ap_sha = 3.0 * tpu_sha;
    } else {
        // C4植物
        // Rubisco限制
        ac_sun = (*bsun) * vcmax_sun;
        ac_sha = (*bsha) * vcmax_sha;
        
        // RuBP限制（光响应）
        aj_sun = qe * aPAR_sun * 4.6;
        aj_sha = qe * aPAR_sha * 4.6;
        
        // PEP羧化酶限制
        ap_sun = kp_sun * max(cisun, 0.0) / pressure;
        ap_sha = kp_sha * max(cisha, 0.0) / pressure;
    }
    
    // ========== 协同限制计算总光合 ==========
    // 阳叶：先耦合AC和AJ
    aquad = theta_cj;
    bquad = -(ac_sun + aj_sun);
    cquad = ac_sun * aj_sun;
    solve_quadratic(aquad, bquad, cquad, &r1, &r2);
    ai = min(r1, r2);
    
    // 再耦合AI和AP
    aquad = theta_ip;
    bquad = -(ai + ap_sun);
    cquad = ai * ap_sun;
    solve_quadratic(aquad, bquad, cquad, &r1, &r2);
    ag_sun = max(0.0, min(r1, r2));
    
    // 阴叶
    aquad = theta_cj;
    bquad = -(ac_sha + aj_sha);
    cquad = ac_sha * aj_sha;
    solve_quadratic(aquad, bquad, cquad, &r1, &r2);
    ai = min(r1, r2);
    
    aquad = theta_ip;
    bquad = -(ai + ap_sha);
    cquad = ai * ap_sha;
    solve_quadratic(aquad, bquad, cquad, &r1, &r2);
    ag_sha = max(0.0, min(r1, r2));
    
    // ========== 净光合速率 ==========
    an_sun = ag_sun - (*bsun) * lmr_sun;
    an_sha = ag_sha - (*bsha) * lmr_sha;
    
    // ========== 负光合处理 ==========
    if (an_sun < 0.0) {
        *gs_mol_sun = medlynint;
        *gs_mol_sun = max((*bsun) * (*gs_mol_sun), 1.0);
        *fvalsun = 0.0;
    }
    
    if (an_sha < 0.0) {
        *gs_mol_sha = medlynint;
        *gs_mol_sha = max((*bsha) * (*gs_mol_sha), 1.0);
        *fvalsha = 0.0;
    }
    
    // 如果都为负，直接返回
    if (an_sun < 0.0 && an_sha < 0.0) {
        return;
    }
    
    // ========== 基于气孔模型计算导度（仅当an>=0时） ==========
    // 阳叶：计算叶面CO₂分压
    if (an_sun >= 0.0) {
        cs_sun = atmosCO2 - 1.4 / gb_mol * an_sun * pressure;
        cs_sun = max(cs_sun, max_cs);
    }
    
    // Medlyn模型
    if (an_sun >= 0.0) {
        term = 1.6 * an_sun / (cs_sun / pressure * 1.0e6);
        aquad = 1.0;
        bquad = -(2.0 * (medlynint * 1.0e-6 + term) + 
                    (medlynslope * term) * (medlynslope * term) / (gb_mol * 1.0e-6 * rh_canopy));
        cquad = medlynint * medlynint * 1.0e-12 +
                (2.0 * medlynint * 1.0e-6 + term * 
                (1.0 - medlynslope * medlynslope / rh_canopy)) * term;
        solve_quadratic(aquad, bquad, cquad, &r1, &r2);
        *gs_mol_sun = max(r1, r2) * 1.0e6;
    }
    
    // 阴叶
    if (an_sha >= 0.0) {
        cs_sha = atmosCO2 - 1.4 / gb_mol * an_sha * pressure;
        cs_sha = max(cs_sha, max_cs);
        
        term = 1.6 * an_sha / (cs_sha / pressure * 1.0e6);
        aquad = 1.0;
        bquad = -(2.0 * (medlynint * 1.0e-6 + term) + 
                    (medlynslope * term) * (medlynslope * term) / (gb_mol * 1.0e-6 * rh_canopy));
        cquad = medlynint * medlynint * 1.0e-12 +
                (2.0 * medlynint * 1.0e-6 + term * 
                (1.0 - medlynslope * medlynslope / rh_canopy)) * term;
        solve_quadratic(aquad, bquad, cquad, &r1, &r2);
        *gs_mol_sha = max(r1, r2) * 1.0e6;
    }
    
    // ========== 计算残差函数值 ==========
    if (an_sun >= 0.0) {
        if (*gs_mol_sun > 0.0) {
            *fvalsun = cisun - atmosCO2 + an_sun * pressure * 
                      (1.4 * (*gs_mol_sun) + 1.6 * gb_mol) / (gb_mol * (*gs_mol_sun));
        } else {
            *fvalsun = cisun - atmosCO2;
        }
    }
    
    if (an_sha >= 0.0) {
        if (*gs_mol_sha > 0.0) {
            *fvalsha = cisha - atmosCO2 + an_sha * pressure * 
                      (1.4 * (*gs_mol_sha) + 1.6 * gb_mol) / (gb_mol * (*gs_mol_sha));
        } else {
            *fvalsha = cisha - atmosCO2;
        }
    }
    // 写回结构体
    veg_var->ac_sun = ac_sun;
    veg_var->ac_sha = ac_sha;
    veg_var->ag_sun = ag_sun;
    veg_var->ag_sha = ag_sha;
    veg_var->aj_sun = aj_sun;
    veg_var->aj_sha = aj_sha;
    veg_var->an_sun = an_sun;
    veg_var->an_sha = an_sha;
    veg_var->ap_sun = ap_sun;
    veg_var->ap_sha = ap_sha;
}

/******************************************************************************
 * @brief    Solve the quadratic equation ax² + bx + c = 0 
 *           and return the two real roots.
 *****************************************************************************/
void solve_quadratic(double a, double b, double c, double *r1, double *r2) {
    double delta = b * b - 4.0 * a * c;
    if (delta < 0.0) {
        // 复数根情况返回实部
        *r1 = -b / (2.0 * a);
        *r2 = *r1;
    } else {
        double sqrt_delta = sqrt(delta);
        *r1 = (-b - sqrt_delta) / (2.0 * a);
        *r2 = (-b + sqrt_delta) / (2.0 * a);
    }
}

/******************************************************************************
 * @brief    Solve the quadratic equation ax² + bx + c = 0 
 *           and return the smaller root.
 *****************************************************************************/
double solve_quadratic_min(double a, double b, double c) {
    double r1, r2;
    solve_quadratic(a, b, c, &r1, &r2);
    return fmin(r1, r2);
}