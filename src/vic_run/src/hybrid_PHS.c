/******************************************************************************
 * @section DESCRIPTION
 *
 * use a hybrid solver to find the root of the ci_func equation for sunlit
 * and shaded leaves.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief  use a hybrid solver to find the root of the ci_func equation 
 *         for sunlit and shaded leaves
 *****************************************************************************/
void hybrid_PHS(double           *x0sun, 
                double           *x0sha,
                double           *vegwp,
                double           *bsun, 
                double           *bsha,
                double            jesun, 
                double            jesha,
                double            CP,
                double            KC,
                double            KO,
                double            thm, 
                double            gb_mol,
                double            qsat_T, 
                double            Qair_over,
                double            pressure, 
                double            air_density,
                double            atmosCO2, 
                double            atmosO2,
                double            lmr_sun, 
                double            lmr_sha,
                double            rh_canopy,
                double            vcmax_sun,
                double            vcmax_sha,
                double            tpu_sun,
                double            tpu_sha,
                double            kp_sun,
                double            kp_sha,
                double           *gs_mol_sun, 
                double           *gs_mol_sha,
                cell_data_struct *cell,
                soil_con_struct  *soil_con,
                veg_var_struct   *veg_var,
                veg_lib_struct   *veg_lib)
{
    // 局部变量
    size_t i, iter1, iter2;
    bool bflag;                     // 是否重新计算水分胁迫
    double x[4];
    double gs0sun, gs0sha;          // 无水分胁迫的气孔导度
    double x1sun, x1sha;            // ci的第二个猜测值
    double f0sun, f0sha;            // f(x0)残差
    double f1sun, f1sha;            // f(x1)残差
    double dxsun, dxsha;            // ci增量
    double b0sun, b0sha;            // 上一次的bsun/bsha
    double dbsun, dbsha;            // bsun/bsha变化量
    double tolsun, tolsha;          // ci收敛容差 (Pa)
    double minf;                    // 最小残差记录
    double minxsun, minxsha;        // 最小残差对应的ci解
    double soilflux;                // 土壤蒸腾通量 (mm/s)
    double xsun, xsha;              // Brent方法返回的解
    // 参数常量
    double toldb = 1e-2;      // bsun/bsha收敛容差
    double eps = 1e-2;        // 相对精度
    double eps1 = 1e-4;       // 绝对精度阈值
    size_t itmax = 3;            // 最大迭代次数
    
    // 初始化
    x1sun = *x0sun;
    x1sha = *x0sha;
    bflag = false;
    b0sun = -1.0;
    b0sha = -1.0;
    gs0sun = 0.0;
    gs0sha = 0.0;
    *bsun = 1.0;
    *bsha = 1.0;
    iter1 = 0;
    iter2 = 0;
    
    // 获取vegwp副本
    for (i = 0; i < 4; i++) {
        x[i] = vegwp[i];
    }
    
    // ========== 外层循环：更新水分胁迫因子 ==========
    while (1) {
        iter1++;
        iter2 = 0;
        
        x1sun = *x0sun;
        *x0sun = max(0.1, x1sun);
        x1sun = 0.99 * x1sun;
        
        x1sha = *x0sha;
        *x0sha = max(0.1, x1sha);
        x1sha = 0.99 * x1sha;
        
        // 设置容差
        tolsun = fabs(x1sun) * eps;
        tolsha = fabs(x1sha) * eps;
        
        // 第一次ci_func调用
        ci_func_PHS(bflag, *x0sun, *x0sha, &f0sun, &f0sha,
                    bsun, bsha, gs_mol_sun, gs_mol_sha, 
                    vegwp, gs0sun, gs0sha, vcmax_sun, vcmax_sha,
                    tpu_sun, tpu_sha, kp_sun, kp_sha, 
                    CP, KC, KO, thm, gb_mol,
                    qsat_T, Qair_over, pressure, air_density,
                    jesun, jesha, atmosCO2, atmosO2,
                    lmr_sun, lmr_sha, rh_canopy,
                    cell, soil_con, veg_var, veg_lib);
        
        // 更新水分胁迫收敛检查变量
        dbsun = b0sun - *bsun;
        dbsha = b0sha - *bsha;
        b0sun = *bsun;
        b0sha = *bsha;
        bflag = false;
        
        // 第二次ci_func调用（计算f(x1)）
        ci_func_PHS(bflag, x1sun, x1sha, &f1sun, &f1sha,
                    bsun, bsha, gs_mol_sun, gs_mol_sha, 
                    vegwp, gs0sun, gs0sha, vcmax_sun, vcmax_sha,
                    tpu_sun, tpu_sha, kp_sun, kp_sha, 
                    CP, KC, KO, thm, gb_mol,
                    qsat_T, Qair_over, pressure, air_density,
                    jesun, jesha, atmosCO2, atmosO2,
                    lmr_sun, lmr_sha, rh_canopy,
                    cell, soil_con, veg_var, veg_lib);
        
        // ========== 内层循环：割线法求解ci ==========
        while (1) {
            // 检查是否已收敛
            if (fabs(f0sun) < eps1 && fabs(f0sha) < eps1) {
                x1sun = *x0sun;
                x1sha = *x0sha;
                break;
            }
            if (fabs(f1sun) < eps1 && fabs(f1sha) < eps1) {
                break;
            }
            
            iter2++;
            
            // 阳叶：割线法更新
            if (f1sun - f0sun == 0.0) {
                dxsun = 0.5 * (x1sun + *x0sun) - x1sun;
            } else {
                dxsun = -f1sun * (x1sun - *x0sun) / (f1sun - f0sun);
            }
            
            // 阴叶：割线法更新
            if (f1sha - f0sha == 0.0) {
                dxsha = 0.5 * (x1sha + *x0sha) - x1sha;
            } else {
                dxsha = -f1sha * (x1sha - *x0sha) / (f1sha - f0sha);
            }
            
            // 更新变量
            *x0sun = x1sun;
            x1sun = x1sun + dxsun;
            *x0sha = x1sha;
            x1sha = x1sha + dxsha;
            
            // 调用ci_func计算新的f值
            ci_func_PHS(bflag, x1sun, x1sha, &f1sun, &f1sha,
                        bsun, bsha, gs_mol_sun, gs_mol_sha, 
                        vegwp, gs0sun, gs0sha, vcmax_sun, vcmax_sha,
                        tpu_sun, tpu_sha, kp_sun, kp_sha, 
                        CP, KC, KO, thm, gb_mol,
                        qsat_T, Qair_over, pressure, air_density,
                        jesun, jesha, atmosCO2, atmosO2,
                        lmr_sun, lmr_sha, rh_canopy,
                        cell, soil_con, veg_var, veg_lib);
            
            // 检查增量收敛
            if (fabs(dxsun) < tolsun && fabs(dxsha) < tolsha) {
                *x0sun = x1sun;
                *x0sha = x1sha;
                break;
            }
            
            // 记录最优解（残差最小的点）
            if (iter2 == 1) {
                minf = fabs(f1sun + f1sha);
                minxsun = x1sun;
                minxsha = x1sha;
            } else {
                if (fabs(f1sun + f1sha) < minf) {
                    minf = fabs(f1sun + f1sha);
                    minxsun = x1sun;
                    minxsha = x1sha;
                }
            }
            
            // 检查函数值收敛
            if (fabs(f1sun) < eps1 && fabs(f1sha) < eps1) {
                break;
            }
            
            // 如果函数值异号（根被包围），切换到Brent方法
            if ((f1sun * f0sun < 0.0) && (f1sha * f0sha < 0.0)) {
                brent_PHS(*x0sun, x1sun, f0sun, f1sun,
                          *x0sha, x1sha, f0sha, f1sha,
                          tolsun, vegwp, 
                          vcmax_sun, vcmax_sha,
                          tpu_sun, tpu_sha, kp_sun, kp_sha,
                          CP, KC, KO, thm, 
                          gb_mol, jesun, jesha, 
                          atmosCO2, atmosO2,
                          lmr_sun, lmr_sha,
                          rh_canopy, qsat_T, Qair_over,
                          pressure, air_density, 
                          &xsun, &xsha, bsun, bsha,
                          gs_mol_sun, gs_mol_sha,
                          cell, soil_con, veg_var, veg_lib);

                *x0sun = xsun;
                *x0sha = xsha;
                break;
            }
            
            // 最大迭代次数保护
            if (iter2 > itmax) {
                x1sun = minxsun;
                x1sha = minxsha;
                ci_func_PHS(bflag, x1sun, x1sha, &f1sun, &f1sha,
                            bsun, bsha, gs_mol_sun, gs_mol_sha, 
                            vegwp, gs0sun, gs0sha, vcmax_sun, vcmax_sha,
                            tpu_sun, tpu_sha, kp_sun, kp_sha, 
                            CP, KC, KO, thm, gb_mol,
                            qsat_T, Qair_over, pressure, air_density,
                            jesun, jesha, atmosCO2, atmosO2,
                            lmr_sun, lmr_sha, rh_canopy,
                            cell, soil_con, veg_var, veg_lib);
                break;
            }
        } // 内层循环结束
        
        // 更新无水分胁迫的气孔导度
        if (*bsun > 0.01) {
            gs0sun = *gs_mol_sun / *bsun;
        }
        if (*bsha > 0.01) {
            gs0sha = *gs_mol_sha / *bsha;
        }
        
        bflag = true;  // true，下次外层循环重新计算水分胁迫
        
        // 检查外层循环收敛
        if (fabs(dbsun) < toldb && fabs(dbsha) < toldb) {
            break;
        }
        
        // 最大迭代次数保护
        if (iter1 > itmax) {
            break;
        }
        
    }
    
    *x0sun = x1sun;
    *x0sha = x1sha;
    
    // 根据最终气孔导度更新植被水势
    getvegwp(vegwp, gb_mol, &gs0sun, &gs0sha, 
             qsat_T, Qair_over, &soilflux,
             pressure, air_density, thm,
             cell, soil_con, veg_var, veg_lib);
    
    // 保存vegwp
    for (i = 0; i < 4; i++) {
        vegwp[i] = x[i];
    }
    
    // 确保蒸腾通量非负
    if (soilflux < 0.0) {
        soilflux = 0.0;
    }
    cell->transp = soilflux;
}