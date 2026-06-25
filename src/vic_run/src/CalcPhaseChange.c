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
                snow_data_struct  *snow,
                soil_con_struct   *soil_con)
{
    double tmp_tkfrz = 0.0;
    double eff_Cs = 0.0;
    double EnergyRes = 0.0;
    double liq_deriv = 0.0;
    double tmp_matric = 0.0;
    double fusion_flux = 0.0;
    double *ice = cell->ice;
    double *liq = cell->liq; 
    double *T = energy->T;
    double *matric = cell->matric;
    double *Cs_node = energy->Cs_node;
    double *Wsat_node = soil_con->Wsat_node;
    // 保存临时变量
    double tmp_ice = ice[nidx];
    double tmp_liq = liq[nidx];
    double total_liq = liq[nidx] + ice[nidx] * CONST_RHOICE / CONST_RHOFW;
    if (total_liq > Wsat_node[nidx]) {
        total_liq = Wsat_node[nidx];
    }
    tmp_matric = SoilWaterRetentionCurve(MATRIC_FLAG, nidx,
                                         total_liq, 0.0, soil_con);
    if (tmp_matric < 0.0) {
        tmp_tkfrz = CONST_TKFRZ * tmp_matric /
                    (CONST_LATICE / CONST_G - tmp_matric) + CONST_TKFRZ;
    } 
    else {
        tmp_tkfrz = CONST_TKFRZ;
    }
    liq_deriv = water_curve_deriv(nidx, tmp_tkfrz, total_liq, tmp_matric, soil_con);
    eff_Cs = Cs_node[nidx] + CONST_RHOFW * CONST_LATICE * liq_deriv;

    if (T[nidx] > tmp_tkfrz && ice[nidx] > 0.0) {
        EnergyRes = Cs_node[nidx] * (T[nidx] - tmp_tkfrz);
        fusion_flux = CONST_RHOICE * CONST_LATICE * ice[nidx];
        if (EnergyRes >= fusion_flux) {
            liq[nidx] = total_liq;
            ice[nidx] = 0.0;
            T[nidx] = tmp_tkfrz + (EnergyRes - fusion_flux) / Cs_node[nidx];
            // 重新计算基质势
            if (liq[nidx] > soil_con->Wpwp_node[nidx]) {
                matric[nidx] = SoilWaterRetentionCurve(MATRIC_FLAG, nidx,
                                                       liq[nidx], 0.0, soil_con);
            }
        }
        else {
            double melted_ice = max(EnergyRes, 0.0) / (CONST_RHOICE * CONST_LATICE);
            liq[nidx] = tmp_liq + melted_ice * CONST_RHOICE / CONST_RHOFW;
            ice[nidx] = tmp_ice - melted_ice;
            T[nidx] = tmp_tkfrz;
            // 用修正后的温度计算冻结状态
            frozen_soil(nidx, T[nidx], liq, ice, matric, soil_con);
        }
    }
    else if (T[nidx] < tmp_tkfrz && ice[nidx] > 0.0 && ice[nidx] > tmp_ice) {
        EnergyRes = Cs_node[nidx] * (tmp_tkfrz - T[nidx]);
        T[nidx] = tmp_tkfrz - EnergyRes / eff_Cs;
        // 用修正后的温度计算冻结状态
        frozen_soil(nidx, T[nidx], liq, ice, matric, soil_con);
    }

    return (0);
}