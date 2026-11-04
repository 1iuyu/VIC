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
double
frozen_soil(size_t             nidx,    // 节点索引
            double             tmp_tkfrz,
            double             soil_T,
            cell_data_struct  *cell,
            soil_con_struct   *soil_con)
{
    extern parameters_struct param;
    double equil_liq;
    double *moist = cell->moist;
    double *Wsat_node = soil_con->Wsat_node;
    double *Wpwp_node = soil_con->Wpwp_node;

    // 通过 Clausius-Clapeyron 方程计算
    double total_potent = CONST_LATICE * ((soil_T - tmp_tkfrz) / soil_T) / CONST_G;
    equil_liq = SoilWaterRetentionCurve(MOIST_FLAG, nidx, 0.0,
                                        total_potent, 0.0, soil_con);

    // Newton-Raphson 迭代
    size_t iter = 0;
    bool converged = false;
    while (iter < 30 && !converged) {
        if (equil_liq < Wpwp_node[nidx]) {
            equil_liq = Wpwp_node[nidx] + param.TOL_A;
        }
        double tmp_mat = SoilWaterRetentionCurve(MATRIC_FLAG, nidx, equil_liq, 
                                                 0.0, 0.0, soil_con);
        double dtheta_dpsi = SoilWaterRetentionCurve(DW_DPSI_FLAG, nidx, equil_liq,
                                                     tmp_mat, 0.0, soil_con);
        if (dtheta_dpsi == 0.0) {
            dtheta_dpsi = SoilWaterRetentionCurve(DW_DPSI_FLAG, nidx, equil_liq,
                                                  total_potent, 0.0, soil_con);
        }
        double error = total_potent - tmp_mat;
        double delta = error / (1.0 / dtheta_dpsi);
        equil_liq += delta;

        if (equil_liq > Wsat_node[nidx]) {
            equil_liq = Wsat_node[nidx];
        }
        double halfway = (equil_liq - delta + Wpwp_node[nidx]) / 2.0;
        if (equil_liq < halfway) {
            equil_liq = halfway;
        }

        if (fabs(delta) < 1.0e-5) {
            converged = true;
        }

        iter++;
    }
    // 平衡态液态水不能超过总水量
    if (equil_liq > moist[nidx]) {
        equil_liq = moist[nidx];
    }
    
    return (equil_liq);
}
