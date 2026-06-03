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
void brent_PHS(double            x1sun, 
               double            x2sun, 
               double            f1sun, 
               double            f2sun, 
               double            x1sha, 
               double            x2sha, 
               double            f1sha, 
               double            f2sha,
               double            tol,
               double           *mat_VEG,
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
               double            jesun, 
               double            jesha,
               double            atmosCO2, 
               double            atmosO2,
               double            lmr_sun, 
               double            lmr_sha,
               double            rh_canopy,
               double            qsat_T, 
               double            Qair_over,
               double            pressure, 
               double            air_density,
               double           *xsun,
               double           *xsha,
               double           *bsun, 
               double           *bsha,
               double           *gs_mol_sun, 
               double           *gs_mol_sha,
               cell_data_struct *cell,
               soil_con_struct  *soil_con,
               veg_var_struct   *veg_var,
               veg_lib_struct   *veg_lib)
{
    // 局部变量
    double gs0sun, gs0sha;              // 无水分胁迫的气孔导度副本
    
    // Brent方法变量（分别对应晴叶和阴叶）
    double a[2], b[2], c[2], d[2], e[2];    // 区间和插值变量
    double fa[2], fb[2], fc[2];             // 函数值
    double pv[2], qv[2], rv[2], sv[2];      // 反二次插值系数
    double tol1[2], xm[2];                  // 收敛容差和中点
    const int nphs = 2;                     // 晴叶(sun=0)和阴叶(sha=1)
    const int itmax = 20;                   // 最大迭代次数
    const double eps = 1e-4;                // 相对误差容限
    bool bflag = false;                     // 是否重新计算水分胁迫
    
    int iter;
    int phase;
    
    // 初始化晴叶变量
    a[0] = x1sun;
    b[0] = x2sun;
    fa[0] = f1sun;
    fb[0] = f2sun;
    
    // 初始化阴叶变量
    a[1] = x1sha;
    b[1] = x2sha;
    fa[1] = f1sha;
    fb[1] = f2sha;
    
    // 检查根是否被包围（函数值必须异号）
    for (phase = 0; phase < nphs; phase++) {
        if ((fa[phase] > 0.0 && fb[phase] > 0.0) || 
            (fa[phase] < 0.0 && fb[phase] < 0.0)) {
            log_err("brent_PHS: root must be bracketed for phase %d (fa=%f, fb=%f)", 
                    phase, fa[phase], fb[phase]);
        }
    }
    
    // 初始化c和fc为b和fb
    c[0] = b[0];
    c[1] = b[1];
    fc[0] = fb[0];
    fc[1] = fb[1];
    
    iter = 0;
    
    while (1) {
        if (iter >= itmax) {
            log_warn("brent_PHS: exceeding maximum iterations (%d)", itmax);
            break;
        }
        iter++;
        
        for (phase = 0; phase < nphs; phase++) {
            // 重新排列a, b, c，确保fb是最小的函数值
            if ((fb[phase] > 0.0 && fc[phase] > 0.0) || 
                (fb[phase] < 0.0 && fc[phase] < 0.0)) {
                c[phase] = a[phase];
                fc[phase] = fa[phase];
                d[phase] = b[phase] - a[phase];
                e[phase] = d[phase];
            }
            
            if (fabs(fc[phase]) < fabs(fb[phase])) {
                a[phase] = b[phase];
                b[phase] = c[phase];
                c[phase] = a[phase];
                fa[phase] = fb[phase];
                fb[phase] = fc[phase];
                fc[phase] = fa[phase];
            }
        }
        
        // 收敛检查
        tol1[0] = 2.0 * eps * fabs(b[0]) + 0.5 * tol;
        tol1[1] = 2.0 * eps * fabs(b[1]) + 0.5 * tol;
        xm[0] = 0.5 * (c[0] - b[0]);
        xm[1] = 0.5 * (c[1] - b[1]);
        
        // 检查是否收敛（两个相位都要收敛）
        if ((fabs(xm[0]) <= tol1[0] || fb[0] == 0.0) &&
            (fabs(xm[1]) <= tol1[1] || fb[1] == 0.0)) {
            *xsun = b[0];
            *xsha = b[1];
            return;
        }
        
        for (phase = 0; phase < nphs; phase++) {
            // 尝试反二次插值或割线法
            if (fabs(e[phase]) >= tol1[phase] && fabs(fa[phase]) > fabs(fb[phase])) {
                // 反二次插值
                sv[phase] = fb[phase] / fa[phase];
                
                if (a[phase] == c[phase]) {
                    // 割线法退化为线性插值
                    pv[phase] = 2.0 * xm[phase] * sv[phase];
                    qv[phase] = 1.0 - sv[phase];
                } 
                else {
                    // 完整反二次插值
                    qv[phase] = fa[phase] / fc[phase];
                    rv[phase] = fb[phase] / fc[phase];
                    pv[phase] = sv[phase] * (2.0 * xm[phase] * qv[phase] * (qv[phase] - rv[phase]) - 
                                            (b[phase] - a[phase]) * (rv[phase] - 1.0));
                    qv[phase] = (qv[phase] - 1.0) * (rv[phase] - 1.0) * (sv[phase] - 1.0);
                }
                
                // 检查是否在边界内
                if (pv[phase] > 0.0) {
                    qv[phase] = -qv[phase];
                }
                pv[phase] = fabs(pv[phase]);
                
                if (2.0 * pv[phase] < min(3.0 * xm[phase] * qv[phase] - fabs(tol1[phase] * qv[phase]),
                                           fabs(e[phase] * qv[phase]))) {
                    // 接受插值结果
                    e[phase] = d[phase];
                    d[phase] = pv[phase] / qv[phase];
                } else {
                    // 插值失败，使用二分法
                    d[phase] = xm[phase];
                    e[phase] = d[phase];
                }
            } else {
                // 区间收缩太慢，使用二分法
                d[phase] = xm[phase];
                e[phase] = d[phase];
            }
            
            // 更新a和fa（将当前b移到a位置）
            a[phase] = b[phase];
            fa[phase] = fb[phase];
            
            // 计算新的试探点
            if (fabs(d[phase]) > tol1[phase]) {
                b[phase] += d[phase];
            } else {
                if (xm[phase] > 0.0) {
                    b[phase] += tol1[phase];
                } else {
                    b[phase] -= tol1[phase];
                }
            }
        }
        
        // 保存当前气孔导度作为无胁迫导度的副本
        gs0sun = *gs_mol_sun;
        gs0sha = *gs_mol_sha;
        
        // 调用ci_func_PHS计算新点的函数值
        ci_func_PHS(bflag, b[0], b[1], 
                    &fb[0], &fb[1],
                    bsun, bsha, 
                    gs_mol_sun, gs_mol_sha,
                    mat_VEG, gs0sun, gs0sha, 
                    vcmax_sun, vcmax_sha, 
                    tpu_sun, tpu_sha,
                    kp_sun, kp_sha, 
                    CP, KC, KO, thm,
                    gb_mol, qsat_T, 
                    Qair_over, pressure, 
                    air_density, jesun, jesha, 
                    atmosCO2, atmosO2, 
                    lmr_sun, lmr_sha,
                    rh_canopy, cell,
                    soil_con, 
                    veg_var, veg_lib);
        
        // 如果两个相位的函数值都为零，收敛
        if (fb[0] == 0.0 && fb[1] == 0.0) {
            break;
        }
    }
    
    // 达到最大迭代次数时输出警告
    if (iter >= itmax) {
        log_warn("brent_PHS: exceeding maximum iterations (%d)", itmax);
    }
    
    // 返回结果
    *xsun = b[0];
    *xsha = b[1];
}