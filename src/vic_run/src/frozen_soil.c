/******************************************************************************
 * @section DESCRIPTION
 * 
 *    This subroutine computes the maximum amount of unfrozen water that
 *    can exist at the current temperature.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    This subroutine computes the maximum amount of unfrozen water that
 *           can exist at the current temperature.
 *****************************************************************************/
int
frozen_soil(size_t     nidx,
            double     T,
            double    *liq,
            double    *ice,
            double    *matric,
            soil_con_struct  *soil_con)
{
    double  tmp_liq, tmp_ice;
    double  error = 1.0;
    double  delt_liq = 1.0;
    double  tmp_matric = 0.0;
    double  tmp_deriv = 0.0;
    double  total_potent = 0.0;
    double *Wsat_node = soil_con->Wsat_node;
    double *Wpwp_node = soil_con->Wpwp_node;

    size_t iter = 0;
    size_t count = 0;
    if (T > CONST_TKFRZ - 1.0e-3) {
        tmp_liq = liq[nidx];
        tmp_ice = 0.0;
    }
    else {
        // ---- 冻结情况：迭代求解平衡态未冻水含量 ----
        // 计算目标总水势（Clausius-Clapeyron 方程）
        total_potent = CONST_LATICE * ((T - CONST_TKFRZ) / T) / CONST_G;
        // 初始猜测：假设纯水在该温度下的未冻水含量
        tmp_liq = SoilWaterRetentionCurve(MOIST_FLAG, nidx, 0.0,
                                         total_potent, soil_con);

        while (iter < 30 && count == 0) {
            // 防止未冻水含量降至残余含水量以下
            if (tmp_liq < Wpwp_node[nidx]) {
                tmp_liq += 10.0e-7;
            }
            // 计算当前未冻水含量对应的基质势
            tmp_matric = SoilWaterRetentionCurve(MATRIC_FLAG, nidx,
                                                 tmp_liq, 0.0, soil_con);
            // 计算基质势对未冻水含量的导数 d(ψ_m)/dθ_l
            tmp_deriv = SoilWaterRetentionCurve(DERIV_FLAG, nidx, tmp_liq,
                                                tmp_matric, soil_con);
            // 如果导数为零（如饱和状态），使用目标水势重新计算
            if (tmp_deriv == 0.0) {
                tmp_deriv = SoilWaterRetentionCurve(DERIV_FLAG, nidx, tmp_liq,
                                                    total_potent, soil_con);
            }
            // 计算误差
            error = total_potent - tmp_matric;
            // Newton-Raphson 迭代
            delt_liq = error / (1.0 / tmp_deriv);
            tmp_liq += delt_liq;

            // ---- 物理约束 ----
            // 不能超过饱和含水量
            if (tmp_liq > Wsat_node[nidx]) {
                tmp_liq = Wsat_node[nidx];
            }
            // 不能越过残余含水量的一半（防止过度干燥）
            double tol_wpwp = (tmp_liq - delt_liq + Wpwp_node[nidx]) / 2.0;
            if (tmp_liq < tol_wpwp) {
                tmp_liq = tol_wpwp;
            }
            // 检查收敛
            if (fabs(delt_liq) < 1.0e-5) {
                count++;
                break;
            }
            iter++;
        }
    }
    // ---- 冻融状态判定与相变量分配 ----
    if (liq[nidx] + ice[nidx] * CONST_RHOICE / CONST_RHOFW < tmp_liq) {
        // 情况 A：实际总水量 < 最大可能液态水量
        // 所有水都能保持液态，无冰生成
        ice[nidx] = 0.0;
        liq[nidx] += ice[nidx] * CONST_RHOICE / CONST_RHOFW;  // 实际液态水就是全部水量
        // 状态：未冻结
    }
    else {
        // 情况 B：实际总水量 >= 最大液态水量
        // 存在冰
        double prev_liq = liq[nidx];  // 保存之前的液态水含量
        tmp_ice = ice[nidx] + (prev_liq - tmp_liq) * CONST_RHOFW / CONST_RHOICE;
        // tmp_liq 保持为平衡态未冻水含量不变
        if (tmp_ice < 0.0) {
            // 数值误差导致冰为负，设为零
            tmp_liq = liq[nidx] + ice[nidx] * CONST_RHOICE / CONST_RHOFW;
            tmp_ice = 0.0;
        }
        liq[nidx] = tmp_liq;
        ice[nidx] = tmp_ice;
    }

    // ---- 重新计算基质势（确保与最终液态水含量一致） ----
    if (tmp_liq > Wpwp_node[nidx]) {
        tmp_matric = SoilWaterRetentionCurve(MATRIC_FLAG, nidx,
                                             tmp_liq, 0.0, soil_con);
        if (tmp_matric >= 0.0) {
            // 数值误差导致非负基质势，强制修正为理论值
            tmp_matric = total_potent;
            if (tmp_matric >= 0.0) {
                tmp_matric = 0.0;
            }
        }
    }
    else {
        // 液态水在残余含水量处，基质势由温度决定
        tmp_matric = total_potent;
    }
    matric[nidx] = tmp_matric;
    
    return (0);
}
