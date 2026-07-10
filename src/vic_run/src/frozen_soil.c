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
frozen_soil(size_t           nidx,
            double           tmp_tkfrz,
            double           T,
            double          *liq,
            double          *ice,
            soil_con_struct *soil_con)
{
    extern parameters_struct param;
    double equil_liq;
    double *Wsat_node = soil_con->Wsat_node;
    double *Wpwp_node = soil_con->Wpwp_node;

    // 热力学诊断
    double total_liq = liq[nidx] + ice[nidx] * CONST_RHOICE / CONST_RHOFW;

    // 通过 Clausius-Clapeyron 方程计算
    double total_potent = CONST_LATICE * ((T - tmp_tkfrz) / T) / CONST_G;
    equil_liq = SoilWaterRetentionCurve(MOIST_FLAG, nidx, 0.0,
                                        total_potent, soil_con);

    // Newton-Raphson 迭代
    int iter = 0, converged = 0;
    while (iter < 30 && !converged) {
        if (equil_liq < Wpwp_node[nidx]) {
            equil_liq = Wpwp_node[nidx] + param.TOL_A;
        }
        double tmp_mat = SoilWaterRetentionCurve(MATRIC_FLAG, nidx,
                                                    equil_liq, 0.0, soil_con);
        double tmp_der = SoilWaterRetentionCurve(DERIV_FLAG, nidx, equil_liq,
                                                    tmp_mat, soil_con);
        if (tmp_der == 0.0) {
            tmp_der = SoilWaterRetentionCurve(DERIV_FLAG, nidx, equil_liq,
                                                total_potent, soil_con);
        }
        double error = total_potent - tmp_mat;
        double delta = error / (1.0 / tmp_der);
        equil_liq += delta;

        if (equil_liq > Wsat_node[nidx]) {
            equil_liq = Wsat_node[nidx];
        }
        double halfway = (equil_liq - delta + Wpwp_node[nidx]) / 2.0;
        if (equil_liq < halfway) {
            equil_liq = halfway;
        }

        if (fabs(delta) < 1.0e-5) {
            converged = 1;
        }

        iter++;
    }
    // 平衡态液态水不能超过总水量
    if (equil_liq > total_liq) {
        equil_liq = total_liq;
    }
    
    return (equil_liq);
}
