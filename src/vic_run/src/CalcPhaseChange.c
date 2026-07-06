/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate the phase change (melting/freezing) of snow water and soil water.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    Calculate the phase change (melting/freezing) of snow and soil 
 *           water, re-adjust the ice and liquid mass, and layer temperature.
 *****************************************************************************/
int
CalcPhaseChange(size_t             nidx,
                energy_bal_struct *energy,
                cell_data_struct  *cell,
                soil_con_struct   *soil_con)
{
    /* Initialize variables */
    double tmp_tkfrz = 0.0;
    double eff_Cs = 0.0;
    double EnergyRes = 0.0;
    double liq_deriv = 0.0;
    double equil_liq = 0.0;
    double tmp_matric = 0.0;
    double fusion_flux = 0.0;
    double *ice = cell->ice;
    double *liq = cell->liq; 
    double *T = energy->T;
    double *matric = cell->matric;
    double *last_Cs = energy->last_Cs;
    double *Wpwp_node = soil_con->Wpwp_node;
    double *Wsat_node = soil_con->Wsat_node;
    // 保存临时变量
    double tmp_ice = ice[nidx];
    double tmp_liq = liq[nidx];
    double tmp_T = T[nidx];
    // 计算冰点和等效热容量（假设所有水为液态）
    double total_liq = liq[nidx] + ice[nidx] * CONST_RHOICE / CONST_RHOFW;
    if (total_liq > Wsat_node[nidx]) {
        total_liq = Wsat_node[nidx];
    }
    tmp_matric = SoilWaterRetentionCurve(MATRIC_FLAG, nidx,
                                         total_liq, 0.0, soil_con);
    if (tmp_matric < 0.0) {
        tmp_tkfrz = CONST_TKTRIP * tmp_matric /
                    (CONST_LATICE / CONST_G - tmp_matric) + CONST_TKTRIP;
    } 
    else {
        tmp_tkfrz = CONST_TKTRIP;
    }
    // 计算等效热容量
    liq_deriv = water_curve_deriv(nidx, tmp_tkfrz, 
                                  total_liq, tmp_matric, soil_con);
    eff_Cs = last_Cs[nidx] + CONST_RHOFW * CONST_LATICE * liq_deriv;

    // 分情况处理
    if (T[nidx] <= tmp_tkfrz) {
        equil_liq = frozen_soil(nidx, tmp_tkfrz, T, liq, ice, soil_con);
        // 热力学上所有水都能保持液态，但需要检查是否有足够冷量来冻结
        if (total_liq < equil_liq) {
            liq[nidx] = total_liq;
            ice[nidx] = 0.0;
        }
        else {
            // 存在冰
            double new_ice = (total_liq - equil_liq) * CONST_RHOFW / CONST_RHOICE;
            double new_liq = equil_liq;
            if (tmp_ice == 0.0 && new_ice > 0.0) {
                EnergyRes = last_Cs[nidx] * (tmp_tkfrz - T[nidx]);
                T[nidx] = tmp_tkfrz - EnergyRes / eff_Cs;
                // 用修正后的温度重新计算平衡态
                equil_liq = frozen_soil(nidx, tmp_tkfrz, T, liq, ice, soil_con);
                liq[nidx] = equil_liq;
                ice[nidx] = (total_liq - equil_liq) * CONST_RHOFW / CONST_RHOICE;
            }
            else {
                liq[nidx] = new_liq;
                ice[nidx] = new_ice;
            }
        }
    }
    else {
        if (ice[nidx] == 0.0) {
            liq[nidx] = total_liq;
            ice[nidx] = 0.0;
        }
        else {
            EnergyRes = last_Cs[nidx] * (tmp_T - tmp_tkfrz);
            fusion_flux = CONST_RHOICE * CONST_LATICE * ice[nidx];
            if (EnergyRes >= fusion_flux) {
                liq[nidx] = total_liq;
                ice[nidx] = 0.0;
                T[nidx] = tmp_tkfrz + (EnergyRes - fusion_flux) / eff_Cs;
            }
            else {
                double melted_ice = EnergyRes / (CONST_RHOICE * CONST_LATICE);
                liq[nidx] = tmp_liq + melted_ice * CONST_RHOICE / CONST_RHOFW;
                ice[nidx] = tmp_ice - melted_ice;
                T[nidx] = tmp_tkfrz;
                // 用修正后的温度确认热力学一致性
                equil_liq = frozen_soil(nidx, tmp_tkfrz, T, liq, ice, soil_con);
                if (equil_liq < total_liq) {
                    liq[nidx] = equil_liq;
                    ice[nidx] = (total_liq - equil_liq) * CONST_RHOFW / CONST_RHOICE;
                }
            }
        }
    }
    // 重新计算基质势
    if (liq[nidx] > Wpwp_node[nidx]) {
        matric[nidx] = SoilWaterRetentionCurve(MATRIC_FLAG, nidx,
                                                liq[nidx], 0.0, soil_con);
    }

    return (0);
}