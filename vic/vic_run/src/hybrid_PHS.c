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
void hybrid_PHS(double *x0sun, double *x0sha,
                double *vegwp,
                double gb_mol,
                double *bsun, double *bsha,
                double jesun, double jesha,
                double thm, double RS_mol,
                double qsat_T, double Qair_over,
                double pressure, double air_density,
                double Canopy_Upper, double matric50, 
                double conduct_max,
                double atmosCO2, double atmosO2,
                double lmr_sun, double lmr_sha,
                double par_sun, double par_sha,
                double rh_can,
                double *gs_mol_sun, double *gs_mol_sha,
                double qsatl, double qaf,
                int *iter1, int *iter2,
                energy_bal_struct *energy,
                cell_data_struct  *cell,
                soil_con_struct   *soil_con,
                veg_var_struct    *veg_var,
                veg_lib_struct    *veg_lib)
{
    // 局部变量
    size_t i;
    double x[4];
    double gs0sun, gs0sha;          // 无水分胁迫的气孔导度
    int bflag;                      // 是否重新计算水分胁迫
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
    const double toldb = 1e-2;      // bsun/bsha收敛容差
    const double eps = 1e-2;        // 相对精度
    const double eps1 = 1e-4;       // 绝对精度阈值
    const int itmax = 3;            // 最大迭代次数
    
    // 初始化
    x1sun = *x0sun;
    x1sha = *x0sha;
    bflag = 0;
    b0sun = -1.0;
    b0sha = -1.0;
    gs0sun = 0.0;
    gs0sha = 0.0;
    *bsun = 1.0;
    *bsha = 1.0;
    *iter1 = 0;
    
    // 获取vegwp副本（用于calcstress）
    for (i = 0; i < 4; i++) {
        x[i] = vegwp[i];
    }
    
    // ========== 外层循环：更新水分胁迫因子 ==========
    while (1) {
        (*iter1)++;
        *iter2 = 0;
        
        x1sun = *x0sun;
        *x0sun = max(0.1, x1sun);
        x1sun = 0.99 * x1sun;
        
        x1sha = *x0sha;
        *x0sha = max(0.1, x1sha);
        x1sha = 0.99 * x1sha;
        
        // 设置容差
        tolsun = fabs(x1sun) * eps;
        tolsha = fabs(x1sha) * eps;
        
        // 第一次ci_func调用（计算f(x0)）
        ci_func_PHS(x, *x0sun, *x0sha, &f0sun, &f0sha,
                    bsun, bsha, bflag,
                    gb_mol, gs0sun, gs0sha,
                    gs_mol_sun, gs_mol_sha,
                    jesun, jesha, atmosCO2, atmosO2,
                    lmr_sun, lmr_sha,
                    par_sun, par_sha, rh_can,
                    qsat_T, Qair_over);
        
        // 更新水分胁迫收敛检查变量
        dbsun = b0sun - *bsun;
        dbsha = b0sha - *bsha;
        b0sun = *bsun;
        b0sha = *bsha;
        bflag = 0;  // false
        
        // 第二次ci_func调用（计算f(x1)）
        ci_func_PHS(x, x1sun, x1sha, &f1sun, &f1sha,
                    bsun, bsha, bflag,
                    gb_mol, gs0sun, gs0sha,
                    gs_mol_sun, gs_mol_sha,
                    jesun, jesha, atmosCO2, atmosO2,
                    lmr_sun, lmr_sha,
                    par_sun, par_sha, rh_can,
                    qsat_T, Qair_over);
        
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
            
            (*iter2)++;
            
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
            ci_func_PHS(x, x1sun, x1sha, &f1sun, &f1sha,
                        bsun, bsha, bflag,
                        gb_mol, gs0sun, gs0sha,
                        gs_mol_sun, gs_mol_sha,
                        jesun, jesha, atmosCO2, atmosO2,
                        lmr_sun, lmr_sha,
                        par_sun, par_sha, rh_can,
                        qsat_T, Qair_over);
            
            // 检查增量收敛
            if (fabs(dxsun) < tolsun && fabs(dxsha) < tolsha) {
                *x0sun = x1sun;
                *x0sha = x1sha;
                break;
            }
            
            // 记录最优解（残差最小的点）
            if (*iter2 == 1) {
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
                brent_PHS(&xsun, *x0sun, x1sun, f0sun, f1sun,
                          &xsha, *x0sha, x1sha, f0sha, f1sha,
                          tolsun, gb_mol,
                          jesun, jesha, atmosCO2, atmosO2,
                          lmr_sun, lmr_sha,
                          par_sun, par_sha, rh_can,
                          gs_mol_sun, gs_mol_sha,
                          bsun, bsha, qsat_T, Qair_over);
                *x0sun = xsun;
                *x0sha = xsha;
                break;
            }
            
            // 最大迭代次数保护
            if (*iter2 > itmax) {
                x1sun = minxsun;
                x1sha = minxsha;
                ci_func_PHS(x, x1sun, x1sha, &f1sun, &f1sha,
                            bsun, bsha, bflag,
                            gb_mol, gs0sun, gs0sha,
                            gs_mol_sun, gs_mol_sha,
                            jesun, jesha, atmosCO2, atmosO2,
                            lmr_sun, lmr_sha,
                            par_sun, par_sha, rh_can,
                            qsat_T, Qair_over);
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
        
        bflag = 1;  // true，下次外层循环重新计算水分胁迫
        
        // 检查外层循环收敛
        if (fabs(dbsun) < toldb && fabs(dbsha) < toldb) {
            break;
        }
        
        // 最大迭代次数保护
        if (*iter1 > itmax) {
            break;
        }
        
    }
    
    *x0sun = x1sun;
    *x0sha = x1sha;
    
    // 根据最终气孔导度更新植被水势
    getvegwp(vegwp, RS_mol, &gs0sun, &gs0sha, 
                qsat_T, Qair_over, &soilflux, 
                Canopy_Upper, matric50, 
                pressure, conduct_max, air_density,
                thm, cell, soil_con, veg_var);
    
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