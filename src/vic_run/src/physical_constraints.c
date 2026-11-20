/******************************************************************************
* @section DESCRIPTION
*
* This routine refines the Newton increments using a backtracking line search
* and imposes physical constraints on the trial state variables.
******************************************************************************/

#include "vic_run.h"

/******************************************************************************
* @brief Refine Newton increments using physical constraints and line search.
******************************************************************************/
int
physical_constraints(size_t             Nstate,
                     double             step_dt,
                     double             fOld,
                     double             slopeInit,
                     double            *fScale,
                     double            *mat_RHS,
                     double            *mat_state,
                     force_data_struct *force,
                     energy_bal_struct *energy,
                     cell_data_struct  *cell,
                     snow_data_struct  *snow,
                     soil_con_struct   *soil_con,
                     veg_var_struct    *veg_var,
                     veg_lib_struct    *veg_lib)
{
    int ErrorFlag;
    const double maxTempInc = 20.0; // Maximum temperature increment (K)
    const double maxPsiInc = 1.0;   // Maximum matric(ψ) increment (m)
    const double maxSnowLiqInc = 0.05;
    size_t i, j, iLine, nidx;
	size_t Nsnow = snow->Nsnow;
	size_t Nnode = cell->Nnode;
    double lambda = 1.0;
    double fNew;
    double delt_t, delt_q;
    double frac_liq, frac_ice;
    double alpha = 1.0e-4;
    double lambda_min = 0.1;
    double increment;
    double delta[MAX_STATES] = {0};
    double residual[MAX_STATES] = {0};
    double *pack_T = snow->pack_T;
    double *soil_T = cell->soil_T;
    double *theta_moist = snow->theta_moist;
    double *moist = cell->moist;
    double *total_matric = cell->total_matric;
    double *Wpwp_node = soil_con->Wpwp_node;
    double *Wsat_node = soil_con->Wsat_node;

    // 保存临时数组
    for (i = 0; i < Nstate; i++) {
        delta[i] = mat_RHS[i];
        residual[i] = 0.0;
    }
    // 线搜索和物理约束
    for (iLine = 0; iLine < 5; iLine++) {

        /* Impose physical constraints */
        for (i = 0, j = 0; i < Nnode; i++, j+=2) {
            // 温度物理约束
            delt_t = lambda * mat_RHS[j];
            if (fabs(delt_t) > maxTempInc) {
                delt_t = copysign(maxTempInc, delt_t);
                if (i < Nsnow) {
                    if (mat_state[j] + delt_t > CONST_TKFRZ) {
                        delt_t = (CONST_TKFRZ - mat_state[j]) * 0.5;
                    }
                }
                else {
                    nidx = i - Nsnow;
                    double Tfrz = calc_tkfrz(nidx, NULL, cell, soil_con);
                    double diff_Tfrz = Tfrz - mat_state[j];
                    if (diff_Tfrz > 0.0) {
                        if (delt_t > diff_Tfrz) {
                            delt_t = diff_Tfrz + 1.0e-7;
                        }
                    }
                    else {
                        if (delt_t < diff_Tfrz) {
                            delt_t = diff_Tfrz - 1.0e-7;
                        }
                    }
                }
            }
            // 水量物理约束
            delt_q = lambda * mat_RHS[j+1];
            if (i < Nsnow) {
                if (fabs(delt_q) > maxSnowLiqInc) {
                    delt_q = copysign(maxSnowLiqInc, delt_q);
                    frac_liq = SnowFreezingCurve(SNOWLIQ_FLAG, pack_T[i], 0.0);
                    frac_ice = 0.0;
                    if (-delt_q > frac_liq) {
                        delt_q = -0.5 * frac_liq;
                    }
                    else if (delt_q > 1.0 - frac_liq - frac_ice) {
                        delt_q = 0.5 * (1.0 - frac_liq - frac_ice);
                    }
                }
                // 更新状态变量
                pack_T[i] = mat_state[j] + delt_t;
                theta_moist[i] = mat_state[j+1] + delt_q;
            }
            else {
                nidx = i - Nsnow;
                if (fabs(delt_q) > maxPsiInc) {
                    delt_q = copysign(maxPsiInc, delt_q);
                    double Tfrz = calc_tkfrz(nidx, NULL, cell, soil_con);
                    if (soil_T[nidx] < Tfrz) {
                        double xConst = CONST_LATICE / (CONST_G * CONST_TKFRZ);
                        double psiliq = xConst * (soil_T[nidx] - CONST_TKFRZ);
                        frac_liq = SoilWaterRetentionCurve(MOIST_FLAG, nidx, 0.0, psiliq, 0.0, soil_con);
                    }
                    else {
                        frac_liq = mat_state[j+1];
                    }
                    frac_ice = (mat_state[j+1] - frac_liq) * (CONST_RHOFW / CONST_RHOICE);

                    // 检查排水量是否超过可用容量，或增加超出可能范围的水量
                    if (-delt_q > frac_liq - Wpwp_node[nidx]) {
                        delt_q = -0.5 * (frac_liq - Wpwp_node[nidx]);
                    }
                    else if (delt_q > Wsat_node[nidx] - frac_liq - frac_ice) {
                        delt_q = 0.5 * (Wsat_node[nidx] - frac_liq - frac_ice);
                    }
                }
                // 更新状态变量
                soil_T[nidx] = mat_state[j] + delt_t;
                total_matric[nidx] = mat_state[j+1] + delt_q;              
            }
        }

        // 更新雪层和土层固液分数比例并调整温度
        CalcPhaseChange(energy, cell, snow, soil_con);

        /* Recalculate residual at trial state */
        update_fluxes(step_dt, force, 
                      energy, 
                      cell, snow,
                      soil_con,
                      veg_var, veg_lib);

        // 计算残差向量
        ErrorFlag = calc_resid_vec(step_dt, residual,
                                   energy, cell, 
                                   snow, soil_con);

        if (ErrorFlag == ERROR) {
            return (ERROR);
        }

        /* Calculate scaled residual norm */
        fNew = 0.0;
        for (i = 0; i < Nstate; i++) {
            increment = fScale[i] * residual[i];
            fNew += 0.5 * increment * increment;
        }

        /* Check Armijo condition */
        if (fNew <= fOld + alpha * lambda * slopeInit) {
            break;
        }

        /* Reduce step length */
        lambda = max(lambda_min, 0.5 * lambda);
    }
    
    // 线搜索失败返回 ERROR，由外层做时间步减半
    if (iLine >= 5) {
        log_err("physical_constraints: line search failed, "
                "fNew=%.6e fOld=%.6e slopeInit=%.6e lambda=%.6e",
                fNew, fOld, slopeInit, lambda);
        return (ERROR);
    }

    return (0);
}