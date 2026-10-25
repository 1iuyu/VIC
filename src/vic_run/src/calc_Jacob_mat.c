/******************************************************************************
* @section DESCRIPTION
*
* This routine was written to handle the various calls and data
* handling needed to solve the various components of the new VIC
* snow code for both the full_energy and water_balance models.
******************************************************************************/

#include "vic_run.h"

/******************************************************************************
* @brief        This routine was written to handle the various calls and data
*               handling needed to solve the various components of the new VIC
*               snow code for both the full_energy and water_balance models.
******************************************************************************/
int
calc_Jacob_mat(double             step_dt,
               double            *mat_LHS,
               force_data_struct *force,
               energy_bal_struct *energy,
               cell_data_struct  *cell,
               snow_data_struct  *snow,
               soil_con_struct   *soil_con)
{
    extern parameters_struct param;

    /* Initialize variables */
	size_t i, j, lidx, idx;
	size_t Nsnow = snow->Nsnow;
    size_t Nsoil = cell->Nsoil;
	size_t Nnode = cell->Nnode;
    double coverage = snow->coverage;
    double zc_node[MAX_SNOWS+1] = {0};
    double fact[MAX_NODES] = {0};
	double *T = energy->T;
    double *last_T = energy->last_T;
    double *liq_deriv = snow->liq_deriv;
    double *dfliq_dT = energy->dfliq_dT;
    double *dCs_dT = energy->dCs_dT;
    double *dCm_dT = energy->dCm_dT;
	double *dz_soil = soil_con->dz_soil;
    double *zc_soil = soil_con->zc_soil;
    double *zc_snow = snow->zc_snow;
    double *dz_snow = snow->dz_snow;
    double *Zsum_snow = snow->Zsum_snow;
    double *Cs_node = energy->Cs_node;
    double *Cm_node = energy->Cm_node;
    double *dWdQ_surf = energy->dWdQ_surf;
    double *dWdE_surf = energy->dWdE_surf;
    double *dEdW_left = energy->dEdW_left;
    double *dEdW_right = energy->dEdW_right;
    double *dEdT_left = energy->dEdT_left;
    double *dEdT_right = energy->dEdT_right;
    double *dWdQ_left = energy->dWdQ_left;
    double *dWdQ_right = energy->dWdQ_right;
    double *dWdE_left = energy->dWdE_left;
    double *dWdE_right = energy->dWdE_right;
    double *dCs_dtheta = energy->dCs_dtheta;
    double *dtheta_dpsi = cell->dtheta_dpsi;
    double *dtheta_dT = energy->dtheta_dT;
    double *dCs_dpsi = energy->dCs_dpsi;
    double *dCm_dpsi = energy->dCm_dpsi;
    double *theta_ice = snow->theta_ice;
    double *theta_liq = snow->theta_liq;
    double deriv_snow = energy->deriv_snow;
    double deriv_soil = energy->deriv_soil;
    double deriv_terms = energy->deriv_terms;
	/* initialization */
    for (i = 0; i < Nsnow; i++) {
        zc_node[i] = zc_snow[i];
        if (cell->IS_VEG) {
            zc_node[Nsnow] = zc_soil[0] + Zsum_snow[Nsnow-1];
        }
        else if (cell->IS_GLAC) {
            zc_node[Nsnow] = cell->h2osfc * 0.25 + Zsum_snow[Nsnow-1];
        }
    }
    // 计算地表水分通量限制
    double *moist = energy->moist;
    double *last_moist = energy->last_moist;
    double deriv_evap = energy->deriv_evap;
    double *dcompress_dpsi = cell->dcompress_dpsi;
    // 计算固定的时间步长系数
    size_t tmp_Nsnow = Nsnow;
    double capr = 0.34;  // heat capacity ratio
    for (i = 0; i < Nnode; i++) {
        if (Nsnow > 0 && i < Nsnow) {
            if (i == 0) {
                fact[i] = 0.5 * (zc_node[i] + capr * zc_node[i+1]) / step_dt;
            }
            else {
                fact[i] = dz_snow[i] / step_dt;
            }
        }
        else if (i == Nsnow && cell->h2osfc > param.TOL_A) {
            fact[i] = 0.5 * cell->h2osfc / step_dt;
            tmp_Nsnow++;
        }
        else {
            if (i == 0) {
                fact[i] = step_dt / 0.5 / (zc_soil[i] + capr * zc_soil[i+1]);
            }
            else {
                lidx = i - tmp_Nsnow;
                fact[i] = step_dt / dz_soil[lidx];
            }
        }
    }

    // ============================================================
    //     构建Jacobian矩阵
    // ============================================================
    for (i = 0, j = 0; i < Nnode; i++, j+=2) {
        if (i < Nsnow) {
            if (i == 0) {
                // 主对角项(本层能量对温度的导数) dE(i)/dT(i)
                idx = get_band_idx(j, j);
                mat_LHS[idx] = (dEdT_left[i] - deriv_snow) * fact[i] + Cs_node[i] + CONST_LATICE * 
                                CONST_RHOFW * dtheta_dT[i] + dCs_dT[i] * (T[i] - last_T[i]) + dCm_dT[i] * 
                                (moist[i] - last_moist[i]) + CONST_LATICE * CONST_RHOFW * dfliq_dT[i] * 
                                (moist[i] - last_moist[i]);
                // 主对角项(本层能量对水量的导数) dE(i)/dθ(i)
                idx = get_band_idx(j, j+1);
                mat_LHS[idx] = (-1.0 + theta_liq[i]) * CONST_LATICE * CONST_RHOFW + dCs_dtheta[i] * 
                                (T[i] - last_T[i]) + Cm_node[i] + fact[i] * dEdW_left[i] + 
                                CONST_LATICE * CONST_RHOFW * (T[i] - last_T[i]) * dfliq_dT[i];
            }
            else {
                // 主对角项(本层能量对温度的导数) dE(i)/dT(i)
                idx = get_band_idx(j, j);
                mat_LHS[idx] = (dEdT_left[i] - dEdT_right[i-1]) * fact[i] + Cs_node[i] + CONST_LATICE * 
                                CONST_RHOFW * dtheta_dT[i] + dCs_dT[i] * (T[i] - last_T[i]) + dCm_dT[i] * 
                                (moist[i] - last_moist[i]) + CONST_LATICE * CONST_RHOFW * dfliq_dT[i] * 
                                (moist[i] - last_moist[i]);
                // 主对角项(本层能量对水量的导数) dE(i)/dθ(i)
                idx = get_band_idx(j, j+1);
                mat_LHS[idx] = (-1.0 + theta_liq[i]) * CONST_LATICE * CONST_RHOFW + dCs_dtheta[i] * (T[i] - 
                                last_T[i]) + Cm_node[i] + fact[i] * (dEdW_left[i] - dEdW_right[i-1]) + 
                                CONST_LATICE * CONST_RHOFW * (T[i] - last_T[i]) * dfliq_dT[i];
            }

            // 上对角项(上层温度对本层温度的导数) dE(i-1)/dT(i)
            if (i > 0) {
                idx = get_band_idx(j-2, j);
                mat_LHS[idx] = fact[i-1] * dEdT_right[i-1];
            }
            // 下对角项(下层温度对本层温度的导数) dE(i+1)/dT(i)
            if (i < Nsnow - 1) {
                idx = get_band_idx(j+2, j);
                mat_LHS[idx] = -fact[i+1] * dEdT_left[i];
            }
            // 上对角项(上层温度对本层水量的导数) dE(i-1)/dθ(i)
            if (i > 0) {
                idx = get_band_idx(j-2, j+1);
                mat_LHS[idx] = fact[i-1] * dEdW_right[i-1];
            }
            // 下对角项(下层温度对本层水量的导数) dE(i+1)/dθ(i)
            if (i < Nsnow - 1) {
                idx = get_band_idx(j+2, j+1);
                mat_LHS[idx] = -fact[i+1] * dEdW_left[i];
            }

            // 主对角项(本层水量对液态水的导数) dW(i)/dθ(i)
            idx = get_band_idx(j+1, j+1);
            mat_LHS[idx] = fact[i] * liq_deriv[i] * theta_ice[i] + 1.0;
            // 下对角项(下层水量对本层水量的导数) dW(i)/dθ(i+1)
            if (i < Nsnow - 1) {
                if (theta_ice[i+1] < param.SNOW_MAX_LIQUID_FRAC) {
                    idx = get_band_idx(j+3, j+1);
                    mat_LHS[idx] = -fact[i+1] * liq_deriv[i] * theta_ice[i+1];
                }
            }
            // 主对角项(本层水量对温度的导数) dW(i)/dT(i)
            idx = get_band_idx(j+1, j);
            mat_LHS[idx] = fact[i] * liq_deriv[i] * dtheta_dT[i];
            // 下对角项(下层水量对本层温度的导数) dW(i)/dT(i+1)
            if (i < Nsnow - 1) {
                if (theta_ice[i+1] < param.SNOW_MAX_LIQUID_FRAC) {
                    idx = get_band_idx(j+3, j);
                    mat_LHS[idx] = -fact[i+1] * liq_deriv[i] * dtheta_dT[i];
                }
            }
        }
        else {
            lidx = i - Nsnow;
            if (lidx == 0) {
                if (cell->h2osfc > param.TOL_A || Nsnow > 0) {
                    // 主对角项(本层能量对温度的导数) dE(i)/dT(i)
                    idx = get_band_idx(j, j);
                    mat_LHS[idx] = (dEdT_left[i] - deriv_soil * (1.0 - coverage)) * fact[i] + Cs_node[i] + 
                                    CONST_LATICE * CONST_RHOFW * dtheta_dT[i] + dCs_dT[i] * (T[i] - last_T[i]) + 
                                    dCm_dT[i] * (moist[i] - last_moist[i]) + CONST_LATICE * CONST_RHOFW * 
                                    dfliq_dT[i] * (moist[i] - last_moist[i]);
                    // 主对角项(本层水量对液态水的导数) dW(i)/dθ(i)
                    idx = get_band_idx(j+1, j+1);
                    mat_LHS[idx] = fact[i] * dWdQ_left[lidx] + 
                                dcompress_dpsi[lidx] * dtheta_dpsi[i];
                    // 主对角项(本层水量对温度的导数) dW(i)/dT(i)
                    idx = get_band_idx(j+1, j);
                    mat_LHS[idx] = fact[i] * (dWdE_left[lidx] - deriv_evap * (1.0 - coverage));
                }
                else {
                    // 主对角项(本层能量对温度的导数) dE(i)/dT(i) 
                    idx = get_band_idx(j, j); // 表土层裸露于大气
                    mat_LHS[idx] = (dEdT_left[i] - deriv_terms) * fact[i] + Cs_node[i] + CONST_LATICE * 
                                    CONST_RHOFW * dtheta_dT[i] + dCs_dT[i] * (T[i] - last_T[i]) + dCm_dT[i] * 
                                    (moist[i] - last_moist[i]) + CONST_LATICE * CONST_RHOFW * dfliq_dT[i] * 
                                    (moist[i] - last_moist[i]);
                    // 主对角项(本层水量对液态水的导数) dW(i)/dθ(i)
                    idx = get_band_idx(j+1, j+1);
                    mat_LHS[idx] = fact[i] * dWdQ_left[lidx] + 
                                dcompress_dpsi[lidx] * dtheta_dpsi[i];
                    // 主对角项(本层水量对温度的导数) dW(i)/dT(i)
                    idx = get_band_idx(j+1, j);
                    mat_LHS[idx] = fact[i] * (dWdE_left[lidx] - deriv_evap);                 
                }
                // 主对角项(本层能量对水量的导数) dE(i)/dθ(i)
                idx = get_band_idx(j, j+1);
                mat_LHS[idx] = Cm_node[i] * dtheta_dpsi[i] + dCs_dpsi[i] * (T[i] - last_T[i]) + dCm_dpsi[i] * 
                                (moist[i] - last_moist[i]) + step_dt * fact[i] * dEdW_left[i];
                // 如果土层存在冰
                if ((moist[i] - last_moist[i]) > 1.0e-12) {
                    mat_LHS[idx] -= CONST_LATICE * CONST_RHOFW * dtheta_dpsi[i];
                }
            }
            else {
                // 主对角项(本层能量对温度的导数) dE(i)/dT(i) 
                idx = get_band_idx(j, j);
                mat_LHS[idx] = (dEdT_left[i] - dEdT_right[i-1]) * fact[i] + Cs_node[i] + CONST_LATICE * 
                                CONST_RHOFW * dtheta_dT[i] + dCs_dT[i] * (T[i] - last_T[i]) + dCm_dT[i] * 
                               (moist[i] - last_moist[i]) + CONST_LATICE * CONST_RHOFW * dfliq_dT[i] * 
                               (moist[i] - last_moist[i]);
                // 主对角项(本层能量对水量的导数) dE(i)/dθ(i)
                idx = get_band_idx(j, j+1);
                mat_LHS[idx] = Cm_node[i] * dtheta_dpsi[i] + dCs_dpsi[i] * (T[i] - last_T[i]) + dCm_dpsi[i] * 
                              (moist[i] - last_moist[i]) + step_dt * fact[i] * (dEdW_left[i] - dEdW_right[i-1]);
                // 如果土层存在冰
                if ((moist[i] - last_moist[i]) > 1.0e-12) {
                    mat_LHS[idx] -= CONST_LATICE * CONST_RHOFW * dtheta_dpsi[i];
                }
                // 主对角项(本层水量对液态水的导数) dW(i)/dθ(i)
                idx = get_band_idx(j+1, j+1);
                mat_LHS[idx] = fact[i] * (dWdQ_left[lidx] - dWdQ_right[lidx-1]) +
                            dcompress_dpsi[lidx] * dtheta_dpsi[i];
                // 主对角项(本层水量对温度的导数) dW(i)/dT(i)
                idx = get_band_idx(j+1, j);
                mat_LHS[idx] = fact[i] * (dWdE_left[lidx] - dWdE_right[lidx-1]);       
            }

            // 上对角项(上层温度对本层温度的导数) dE(i-1)/dT(i)
            if (lidx > 0) {
                idx = get_band_idx(j-2, j);
                mat_LHS[idx] = fact[i-1] * dEdT_right[i-1];
            }
            // 下对角项(下层温度对本层温度的导数) dE(i+1)/dT(i)
            if (i < Nnode - 1) {
                idx = get_band_idx(j+2, j);
                mat_LHS[idx] = -fact[i+1] * dEdT_left[i];
            }

            // 上对角项(上层能量对本层水量的导数) dE(i-1)/dθ(i)
            if (lidx > 0) {
                idx = get_band_idx(j-2, j+1);
                mat_LHS[idx] = fact[i-1] * dEdW_right[i-1];
            }
            // 下对角项(下层能量对本层水量的导数) dE(i+1)/dθ(i)
            if (lidx < Nsoil - 1) {
                idx = get_band_idx(j+2, j+1);
                mat_LHS[idx] = -fact[i+1] * dEdW_left[i];
            }
            // 上对角项(上层水量对本层水量的导数) dW(i-1)/dθ(i)
            if (lidx > 0) {
                idx = get_band_idx(j-1, j+1);
                mat_LHS[idx] = fact[i-1] * dWdQ_right[lidx-1];
            }
            // 下对角项(下层水量对本层水量的导数) dW(i+1)/dθ(i)
            if (i < Nnode - 1) {
                idx = get_band_idx(j+3, j+1);
                mat_LHS[idx] = -fact[i+1] * dWdQ_left[lidx];
            }
            // 上对角项(上层水量对本层温度的导数) dW(i-1)/dT(i)
            if (lidx > 0) {
                idx = get_band_idx(j-1, j);
                mat_LHS[idx] = fact[i-1] * dWdE_right[lidx-1];
            }
            // 下对角项(下层水量对本层温度的导数) dW(i+1)/dT(i) 
            if (i < Nnode - 1) {
                idx = get_band_idx(j+3, j);
                mat_LHS[idx] = -fact[i+1] * dWdE_left[lidx];
            }
            // 交叉项(地表入渗对于温度的导数)
            idx = get_band_idx(Nsnow*2+1, j);
            mat_LHS[idx] -= fact[Nsnow] * dWdE_surf[lidx];
            // 交叉项(地表入渗对于水量的导数)
            idx = get_band_idx(Nsnow*2+1, j+1);
            mat_LHS[idx] -= fact[Nsnow] * dWdQ_surf[lidx];
        }
    }
    
    return (0);
}