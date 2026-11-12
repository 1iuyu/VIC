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
CalcPhaseChange(energy_bal_struct *energy,
                cell_data_struct  *cell,
                snow_data_struct  *snow,
                soil_con_struct   *soil_con)
{
    /* Initialize variables */
    size_t i, nidx, iter;
    size_t Nsnow = snow->Nsnow;
    size_t Nnode = cell->Nnode;
    double tmp_ice, tmp_T, tmp_tkfrz;
    double EnergyRes = 0.0;
    double equil_liq = 0.0;
    double tempInc = 0.0;
    double frac_liq = 0.0;
    double critDiff = 0.0;
    double ice_trial = 0.0;
    double *ice = cell->ice;
    double *liq = cell->liq;
    double *soil_T = cell->soil_T;
    double *pack_T = snow->pack_T;
    double *theta_ice = snow->theta_ice;
    double *theta_liq = snow->theta_liq;
    double *theta_moist = snow->theta_moist;
    double *matric = cell->matric;
    double *dfliq_dT = energy->dfliq_dT;
    double *Cs_node = energy->Cs_node;
    double *moist = cell->moist;
    double *Wpwp_node = soil_con->Wpwp_node;
    double *Wsat_node = soil_con->Wsat_node;
    for (i = 0; i < Nnode; i++) {
        // 保存临时变量
        if (i < Nsnow) {
            tmp_ice = theta_ice[i];
            tmp_T = pack_T[i];
            tmp_tkfrz = CONST_TKFRZ;
        }
        else {
            nidx = i - Nsnow;
            tmp_ice = ice[i];
            tmp_T = soil_T[i];
            tmp_tkfrz = calc_tkfrz(nidx, NULL, cell, soil_con);
        }

        /* Initialize temperature iteration */
        double xTemp = tmp_T;
        double tempMin = tmp_T - 10.0;
        double tempMax = tmp_T + 10.0;

        /* ============================================================
        *     Energy-conserving phase change iteration (SUMMA method)
        * ============================================================ */
        for (iter = 0; iter < 100; iter++) {

            /* Restrict temperature using the iteration bracket */
            if (xTemp <= tempMin || xTemp >= tempMax) {
                xTemp = 0.5 * (tempMin + tempMax);
            }
            if (i < Nsnow) {
                frac_liq = SnowFreezingCurve(SNOWLIQ_FLAG, pack_T[i], 0.0);
                theta_liq[i] = theta_moist[i] * frac_liq;
                theta_ice[i] = theta_moist[i] * (1.0 - frac_liq) * (CONST_RHOFW / CONST_RHOICE);
            }
            else {
                /* Calculate equilibrium liquid water content */
                equil_liq = frozen_soil(nidx, tmp_tkfrz, xTemp, cell, soil_con);
                equil_liq = max(Wpwp_node[nidx], min(equil_liq, Wsat_node[nidx]));

                /* Calculate ice content (water equivalent) */
                ice_trial = (moist[nidx] - equil_liq) * (CONST_RHOFW / CONST_RHOICE);
            }

            /* Calculate energy residual */
            EnergyRes = -Cs_node[i] * (xTemp - tmp_T)
                        + CONST_RHOFW * CONST_LATICE * (ice_trial - tmp_ice);

            /* Update temperature bracket */
            if (EnergyRes < 0.0) {
                tempMax = min(xTemp, tempMax);
            }
            else {
                tempMin = max(tempMin, xTemp);
            }

            /* Calculate Newton iteration increment */
            tempInc = EnergyRes / (Cs_node[i] + CONST_RHOFW * CONST_LATICE * dfliq_dT[i]);
            
            /* Check convergence */
            if (fabs(EnergyRes) < 1.0e-4 || fabs(tempInc) < 1.0e-6) {
                break;
            }

            /* Handle the discontinuity at the critical freezing temperature */
            critDiff = tmp_tkfrz - xTemp;

            if (critDiff > 0.0) {
                if (tempInc > critDiff) {
                    tempInc = critDiff + 1.0e-7;
                }
            }
            else {
                if (tempInc < critDiff) {
                    tempInc = critDiff - 1.0e-7;
                }
            }

            /* Update temperature */
            xTemp += tempInc;
        }

        /* ============================================================
        *     Update final liquid water and ice contents
        * ============================================================ */

        /* Calculate equilibrium liquid water at the final temperature */
        if (i < Nsnow) {
            frac_liq = SnowFreezingCurve(SNOWLIQ_FLAG, xTemp, 0.0);
            pack_T[i] = xTemp;
            theta_liq[i] = theta_moist[i] * frac_liq;
            theta_ice[i] = theta_moist[i] * (1.0 - frac_liq) * (CONST_RHOFW / CONST_RHOICE);
        }
        else {
            equil_liq = frozen_soil(nidx, tmp_tkfrz, xTemp, cell, soil_con);
            equil_liq = max(Wpwp_node[nidx], min(equil_liq, Wsat_node[nidx]));

            /* Update temperature */
            soil_T[nidx] = xTemp;

            /* Update liquid water and ice */
            liq[nidx] = equil_liq;
            ice[nidx] = (moist[nidx] - equil_liq) * CONST_RHOFW / CONST_RHOICE;

            /* Recalculate matric potential */
            if (liq[nidx] > Wpwp_node[nidx]) {
                matric[nidx] = SoilWaterRetentionCurve(MATRIC_FLAG, nidx, liq[nidx], 0.0, 0.0, soil_con);
            }
        }
    }

    return (0);
}