/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate snow and soil layer temperature.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief  
 * Compute bare ground soil and snow surface temperatures.
 *****************************************************************************/
int
SoilTemperature(double   		   step_dt,
                double             pressure,
				cell_data_struct  *cell,
				energy_bal_struct *energy,
				snow_data_struct  *snow,
				soil_con_struct   *soil_con)

{	
    extern parameters_struct param;
    /* Initialize variables */
	size_t i, lidx;
	size_t Nsnow = snow->Nsnow;
    size_t Nsoil = cell->Nsoil;
	size_t Nnode = cell->Nnode;
    double coverage = snow->coverage;
    double frac_h2o = cell->frac_h2o;
    double fact[MAX_NODES];
	double mat_A[MAX_NODES];
	double mat_B[MAX_NODES];
	double mat_C[MAX_NODES];
	double mat_RHS[MAX_NODES];
    double FLOW[MAX_SOILS];
    double EPSLON[MAX_SOILS];
	double *T = energy->T;
	double *soil_T = cell->soil_T;
	double *pack_T = snow->pack_T;
	double *dz_soil = soil_con->dz_soil;
    double *zc_soil = soil_con->zc_soil;
    double *dz_snow = snow->dz_snow;
    double *zc_snow = snow->zc_snow;
    double *Zsum_snow = snow->Zsum_snow;
    double *Zsum_soil = soil_con->Zsum_soil;
	double *kappa_int = energy->kappa_int;
    double *Cs_node = energy->Cs_node;
    double *last_Cs = energy->last_Cs;
    double *last_T = energy->last_T;
    double *matric = cell->matric;
    double *pack_liq = snow->pack_liq;
    double *Wpwp_node = soil_con->Wpwp_node;
    double *last_matric = cell->last_matric;
    double deriv_terms = energy->deriv_terms;
	/* initialization */
    for (i = 0; i < MAX_NODES; i++) {
        fact[i] = 0.0;
        FLOW[i] = 0.0;
        EPSLON[i] = 0.0;
        mat_A[i] = 0.0;
        mat_B[i] = 0.0;
        mat_C[i] = 0.0;
        mat_RHS[i] = 0.0;
    }
    // 计算地表水分通量限制
    size_t IFLAG = 0;
    double Qmax = 0.0;
    double *ice = cell->ice;
    double *liq = cell->liq;
    double *last_ice = cell->last_ice;
    double *last_liq = cell->last_liq;
    double *liquid_flux = cell->liquid_flux;
	double *vapor_flux = cell->vapor_flux;
    double *conv_vapor = cell->conv_vapor;
    double *deric_vapor = cell->deriv_vapor;
    double *conduct_int = cell->conduct_int;
    double *Wsat_node = soil_con->Wsat_node;
    for (i = Nsoil - 2; i > 0; i--) {
        if (ice[i] > 0.0) {
            Qmax = (Wsat_node[i] - liq[i] - ice[i]) * dz_soil[i] / step_dt;
            if (Qmax < 0.0) {
                Qmax = 0.0;
            }
            // 检查净流入水量是否超过最大容量
            if (liquid_flux[i-1] - liquid_flux[i] > Qmax) {
                IFLAG = 1;  
                if (liquid_flux[i-1] > 0.0) {
                    // 上层水向下流 - 限制上层通量
                    conduct_int[i-1] = 0.0;
                    liquid_flux[i-1] = liquid_flux[i] + Qmax;
                    
                    if (liquid_flux[i-1] < 0.0) {
                        liquid_flux[i-1] = 0.0;
                        liquid_flux[i] = liquid_flux[i-1] - Qmax;
                        if (liquid_flux[i] > 0.0) {
                            liquid_flux[i] = 0.0;
                        }
                    }
                } else {
                    // 下层水向上流 - 限制下层通量
                    conduct_int[i] = 0.0;
                    liquid_flux[i] = liquid_flux[i-1] - Qmax;
                    if (liquid_flux[i] > 0.0) {
                        liquid_flux[i] = 0.0;
                    }
                }
            }
        }
        // 设置水通量方向指示符（用于平流热输送）
        if (liquid_flux[i] > 0.0) {
            FLOW[i] = 1.0;   // 水向下流动，热量向下输送
            EPSLON[i] = 0.0;
        } 
        else {
            FLOW[i] = 0.0;   // 水向上流动，热量向上输送
            EPSLON[i] = -1.0;
        }
    }
    // 处理地表节点
    if (liquid_flux[0] > 0.0) {
        FLOW[0] = 1.0;
        EPSLON[0] = 0.0;
    }
    else {
        FLOW[0] = 0.0;
        EPSLON[0] = -1.0;
    }
    // 地表冰层限制
    if (ice[0] > 0.0) {
        Qmax = (Wsat_node[0] - liq[0] - ice[0]) * dz_soil[0] / step_dt;
        if (Qmax < 0.0) {
            Qmax = 0.0;
        }
        if (-liquid_flux[0] > Qmax) {
            IFLAG = 1;
            conduct_int[0] = 0.0;
            liquid_flux[0] = -Qmax;
        }
    }
    // 如果进行了限制，二次调整
    if (IFLAG > 0) {
        for (i = 1; i < Nsoil - 2; i++) {
            if (ice[i] > 0.0) {
                Qmax = (Wsat_node[i] - liq[i] - ice[i]) * dz_soil[i] / step_dt;
                if (Qmax < 0.0) {
                    Qmax = 0.0;
                }
                if (liquid_flux[i-1] - liquid_flux[i] > Qmax) {
                    if (liquid_flux[i] < 0.0) {
                        // 下层水向上流 - 限制下层通量
                        conduct_int[i] = 0.0;
                        liquid_flux[i] = liquid_flux[i-1] - Qmax;
                        if (liquid_flux[i] > 0.0) {
                            liquid_flux[i] = 0.0;
                            liquid_flux[i-1] = liquid_flux[i] + Qmax;
                            if (liquid_flux[i-1] < 0.0) {
								liquid_flux[i-1] = 0.0;
							}
                        }
                    } 
                    else {
                        // 上层水向下流 - 限制上层通量
                        conduct_int[i-1] = 0.0;
                        liquid_flux[i-1] = liquid_flux[i] + Qmax;
                        if (liquid_flux[i-1] < 0.0) {
                            liquid_flux[i-1] = 0.0;
                        }
                    }
                }
            }
        }
    }
    /* Soil and ice thermal properties for the layer */
    prepare_full_energy(pressure, 
                        cell, energy,
                        snow, soil_con);
    
    double capr = 0.34;  // heat capacity ratio                   
    for (i = 0; i < Nsnow; i++) {
		if (i == 0) {
			fact[i] = (0.5 * (zc_snow[i] - Zsum_snow[i] + 
						capr * (zc_snow[i+1] - Zsum_snow[i]))) / step_dt;
		}
		else if (i <= Nsnow - 1) {
			fact[i] = dz_snow[i] / step_dt;
		}
    }
    for (i = 0; i < Nsoil; i++) {
        lidx = Nsnow + i;
        if (Nsnow == 0) {
			fact[lidx] = (0.5 * (zc_soil[i] - Zsum_soil[i] + 
						capr * (zc_soil[i+1] - Zsum_soil[i]))) / step_dt;            
        }
        else {
            fact[lidx] = dz_soil[i] / step_dt;
        }
    }

    double grnd_flux = energy->NetShortGrnd + energy->longwave - 
                        energy->sensible - energy->latent + energy->advection;
    energy->grnd_flux = grnd_flux;
    size_t tmp_Nsnow = Nsnow;
    // ============================================================
    // 第一部分：雪层能量平衡矩阵
    // ============================================================
    for (i = 0; i < Nsnow; i++) {
        if (i == 0) {
            if (pack_liq[i] > 0.0) {
                mat_A[i] = 0.0;
                mat_B[i] = deriv_terms - CONST_RHOFW * CONST_LATICE / step_dt;
                mat_C[i] = 0.0;
            }
            else {
                mat_A[i] = 0.0;
                mat_B[i] = deriv_terms
                        - kappa_int[i]
                        - CONST_LATSUB * deric_vapor[i] * conv_vapor[i]
                        - fact[i] * Cs_node[i];
                mat_C[i] = kappa_int[i]
                        + CONST_LATSUB * deric_vapor[i+1] * conv_vapor[i];
            }
            mat_RHS[i] = grnd_flux
                    - kappa_int[i] * (T[i] - T[i+1])
                    - fact[i] * last_Cs[i] * (T[i] - last_T[i])
                    - CONST_RHOFW * CONST_LATICE * (liq[i] - last_liq[i]) / step_dt
                    - CONST_LATSUB * vapor_flux[i];
        }
        else {
            if (pack_liq[i] > 0.0) {
                mat_A[i] = 0.0;
                mat_B[i] = -CONST_RHOFW * CONST_LATICE / step_dt;
                mat_C[i] = 0.0;
            }
            else {
                mat_A[i] = kappa_int[i-1]
                        + CONST_LATSUB * deric_vapor[i-1] * conv_vapor[i-1];
                if (i < Nsnow - 1) {
                    mat_B[i] = -(kappa_int[i-1] + kappa_int[i])
                            - CONST_LATSUB * deric_vapor[i] *
                            (conv_vapor[i-1] + conv_vapor[i])
                            - fact[i] * Cs_node[i];
                    mat_C[i] = kappa_int[i]
                            + CONST_LATSUB * deric_vapor[i+1] * conv_vapor[i];
                }
                else {
                    /* 最底层雪，与中间层（水/冰川）连接 */
                    mat_B[i] = -kappa_int[i-1]
                            - kappa_int[i]
                            - CONST_LATSUB * deric_vapor[i] *
                            conv_vapor[i-1]
                            - fact[i] * Cs_node[i];
                    mat_C[i] = kappa_int[i];
                }
            }
            if (i < Nsnow - 1) {
                mat_RHS[i] =
                    kappa_int[i-1] * (T[i-1] - T[i])
                - kappa_int[i] * (T[i] - T[i+1])
                - fact[i] * last_Cs[i] * (T[i] - last_T[i])
                - CONST_RHOFW * CONST_LATICE * (liq[i] - last_liq[i]) / step_dt
                - CONST_LATSUB * (vapor_flux[i] - vapor_flux[i-1]);
            }
            else {
                /* 最底层雪，与中间层连接 */
                mat_RHS[i] =
                    kappa_int[i-1] * (T[i-1] - T[i])
                - kappa_int[i] * (T[i] - T[i+1])
                - fact[i] * last_Cs[i] * (T[i] - last_T[i])
                - CONST_RHOFW * CONST_LATICE * (liq[i] - last_liq[i]) / step_dt
                - CONST_LATSUB * (vapor_flux[i] - vapor_flux[i-1]);
            }
        }
    }
    for (i = 0; i < Nsnow; i++) {
        if (i == Nsnow - 1) {
            if (pack_liq[i] > 0.0) {
                mat_A[i] = 0.0;
                mat_B[i] = deriv_terms - CONST_RHOFW * CONST_LATICE / step_dt;
                mat_C[i] = 0.0;
            }
            else {
                mat_A[i] = 0.0;
                mat_B[i] = deriv_terms - kappa_int[i] - CONST_LATSUB * 
                                    deric_vapor[i] * conv_vapor[i] - fact[i] * Cs_node[i];
                mat_C[i] = kappa_int[i] + CONST_LATSUB * deric_vapor[i+1] * conv_vapor[i];
            }
            mat_RHS[i] = grnd_flux - kappa_int[i] * (T[i] - T[i+1]) - fact[i] * last_Cs[i] * 
                         (T[i] - last_T[i]) - CONST_RHOFW * CONST_LATICE * 
                         (liq[i] - last_liq[i]) / step_dt - CONST_LATSUB * vapor_flux[i];
        }
        else if (i < Nsnow - 1) {
            if (pack_liq[i] > 0.0) {
                mat_A[i] = 0.0;
                mat_B[i] = -CONST_RHOFW * CONST_LATICE / step_dt;
                mat_C[i] = 0.0;
            }
            else {
                mat_A[i] = kappa_int[i-1] + CONST_LATSUB * deric_vapor[i-1] * conv_vapor[i-1];
                mat_B[i] = -(kappa_int[i-1] + kappa_int[i]) - CONST_LATSUB * deric_vapor[i] * 
                            (conv_vapor[i-1] + conv_vapor[i]) - fact[i] * Cs_node[i];
                mat_C[i] = kappa_int[i] + CONST_LATSUB * deric_vapor[i+1] * conv_vapor[i];
            }
            mat_RHS[i] = kappa_int[i-1] * (T[i-1] - T[i]) - kappa_int[i] * (T[i] - T[i+1]) - 
                         fact[i] * last_Cs[i] * (T[i] - last_T[i]) - CONST_RHOFW * CONST_LATICE * 
                         (liq[i] - last_liq[i]) / step_dt - CONST_LATSUB * (vapor_flux[i] - vapor_flux[i-1]);
        }
        else {
            if (pack_liq[i] > 0.0) {
                mat_A[i] = 0.0;
                mat_B[i] = -CONST_RHOFW * CONST_LATICE / step_dt;
                mat_C[i] = 0.0;
            } else {
                mat_A[i] = kappa_int[i-1] + CONST_LATSUB * deric_vapor[i-1] * conv_vapor[i-1];
                mat_B[i] = -kappa_int[i-1] - CONST_LATSUB * deric_vapor[i] * conv_vapor[i-1]
                        - fact[i] * Cs_node[i];
                mat_C[i] = 0.0;
            }
            mat_RHS[i] = kappa_int[i-1] * (T[i-1] - T[i]) - fact[i] * last_Cs[i] * (T[i] - last_T[i]) -
                         CONST_RHOFW * CONST_LATICE * (liq[i] - last_liq[i]) / step_dt -
                         CONST_LATSUB * (vapor_flux[i] - vapor_flux[i-1]);
        }
    }
    if (cell->h2osfc > param.TOL_A) {
        lidx = Nsnow - 1;
        if (pack_liq[lidx] = 0.0) {
            mat_A[Nsnow] = kappa_int[lidx];
            mat_B[Nsnow] = -kappa_int[lidx];
            mat_C[Nsnow] = kappa_int[lidx];
            mat_RHS[Nsnow] = kappa_int[lidx] * (T[lidx] - T[Nsnow]);
        }
    }
    // ============================================================
    // 第二部分：土层能量平衡矩阵
    // ============================================================
    for (i = 0; i < Nsoil; i++) {
        lidx = tmp_Nsnow + i;
        double adv_left = 0.0;
        double adv_right = 0.0;
        double latent_term = 0.0;
        double ice_term = 0.0;
        double trans_right = 0.0;   // 右侧界面传输系数
        double trans_left = 0.0;    // 左侧界面传输系数
        // 当前层的水分通量和冰量对热输送的影响
        if (i >= Nsnow && i < Nnode - 2) {
            adv_right = CONST_RHOFW * CONST_CPFWICE * liquid_flux[i] + 
                        CONST_CPWV * vapor_flux[i];
            trans_right = kappa_int[i] + EPSLON[i] * adv_right;
        }
        // 上一层的水分通量和冰量对热输送的影响
        if (i > Nsnow && i < Nnode - 1) {
            adv_left = CONST_RHOFW * CONST_CPFWICE * liquid_flux[i-1] +
                    CONST_CPWV * vapor_flux[i-1];
            latent_term = vapor_flux[i] - vapor_flux[i-1];
            trans_left = kappa_int[i-1] + FLOW[i-1] * adv_left;
        }
        if (i >= Nsnow && i <= Nnode - 2) {
            ice_term = CONST_RHOICE * CONST_LATICE * (ice[i] - last_ice[i]);
        }
        if (i == Nsnow) {
            if (Nsnow == 0 && cell->IS_VEG) {
                mat_A[i] = 0.0;
                mat_B[i] = deriv_terms - trans_right - fact[i] * Cs_node[i];
                mat_C[i] = trans_right;
                mat_RHS[i] = grnd_flux - trans_right * (T[i] - T[i+1]) - CONST_LATVAP *
                        vapor_flux[i] - fact[i] * (last_Cs[i] * (T[i] - last_T[i]) - ice_term);
            }
            else if (Nsnow > 0) {
                trans_left = kappa_int[Nsnow-1] + CONST_LATVAP * conv_vapor[Nsnow-1] * deric_vapor[Nsnow-1];
                mat_A[i] = coverage * trans_left;
                mat_B[i] = deriv_terms - coverage * trans_left - trans_right - fact[i] * Cs_node[i];
                mat_C[i] = trans_right;
                mat_RHS[i] = coverage * trans_left * (T[i-1] - T[i]) - trans_right * (T[i] - T[i+1]) -
                            CONST_LATVAP * vapor_flux[Nsnow-1] - fact[i] * (last_Cs[i] * (T[i] - last_T[i]) - ice_term);
            }
            else if (cell->IS_GLAC) {

            }
            else if (cell->IS_WET) {
                mat_A[i] = frac_h2o * trans_left;
                mat_B[i] = -frac_h2o * trans_left - trans_right - fact[i] * Cs_node[i];
                mat_C[i] = trans_right;
                mat_RHS[i] = frac_h2o * trans_left * (cell->h2osfc_T - T[i]) 
                            - trans_right * (T[i] - T[i+1]) 
                            - CONST_LATVAP * vapor_flux[i]
                            - fact[i] * (last_Cs[i] * (T[i] - last_T[i]) - ice_term);                
            }
        }
        else if (i <= Nnode - 2) {
            mat_A[i] = trans_left;
            mat_B[i] = -trans_left - trans_right - fact[i] * Cs_node[i];
            mat_C[i] = trans_right;
            mat_RHS[i] = trans_left * (T[i-1] - T[i]) - trans_right * (T[i] - T[i+1]) -
                         CONST_LATVAP * latent_term - fact[i] * (last_Cs[i] * (T[i] - last_T[i]) - ice_term);
        }
        else if (i == Nnode - 1) {
            mat_A[i] = trans_left;
            mat_B[i] = -trans_left - fact[i] * Cs_node[i];
            mat_C[i] = 0.0;
            mat_RHS[i] = trans_left * (T[i-1] - T[i]) - fact[i] * (last_Cs[i] * (T[i] - last_T[i]));
        }
    }
    // 修正能量平衡方程的系数
    double tmp_T = 0.0;
    double liq_deriv1 = 0.0;
    double liq_deriv2 = 0.0;
    double tmp_deriv = 0.0;
    double tmp_matric = 0.0;
    for (i = 0; i < Nsoil; i++) {
        if (ice[i] > 0.0) {
            if (last_ice[i] > 0.0) {
                tmp_T = last_T[i];
            }
            else {
                tmp_T = CONST_TKTRIP * CONST_LATICE / CONST_G / 
                                (CONST_LATICE / CONST_G - last_matric[i]);
            }
            liq_deriv1 = water_curve_deriv(i, tmp_T,
                                           last_liq[i], 
                                           last_matric[i],
                                           soil_con);
            if (matric[i] > 0.0) {
                tmp_T = T[i];
                tmp_matric = CONST_LATICE * (tmp_T - CONST_TKFRZ) / tmp_T / CONST_G;
                liq_deriv2 = water_curve_deriv(i, tmp_T, 
                                               liq[i], 
                                               tmp_matric,
                                               soil_con);
                tmp_deriv = SoilWaterRetentionCurve(DERIV_FLAG, i, liq[i],
                                                    tmp_matric, soil_con);
            }
            else {
                liq_deriv2 = water_curve_deriv(i, T[i], liq[i], 
                                               matric[i], 
                                               soil_con);
                tmp_deriv = SoilWaterRetentionCurve(DERIV_FLAG, i, liq[i],
                                                    matric[i], soil_con);
            }
            if (i == 0) {
                mat_B[i] -= fact[i] * 0.5 * CONST_RHOFW * CONST_LATICE * (liq_deriv1 + liq_deriv2) -
                        CONST_RHOFW * CONST_LATICE * conduct_int[i] * liq_deriv2 / liq_deriv1;
            }
            else if (i < Nsoil - 1) {
                mat_B[i] -= fact[i] * 0.5 * CONST_RHOFW * CONST_LATICE * (liq_deriv1 + liq_deriv2) -
                    CONST_RHOFW * CONST_LATICE * (conduct_int[i-1] + conduct_int[i]) * liq_deriv2 / tmp_deriv;
            }
        }
    }

	// 调用LAPACK函数求解线性方程组
	int info = LAPACKE_dgtsv(LAPACK_COL_MAJOR, Nnode,
							 1, mat_A+1, mat_B, mat_C, mat_RHS, Nnode);

	if (info != 0) {
		return ERROR;
	}

	/* 检查收敛 (温度变化) */
	double max_diff = 0.0;
	for (i = 0; i < Nnode; i++) {
		double diff = mat_RHS[i];
        T[i] -= diff;
        if (T[i] < CONST_TKFRZ && i >= Nsnow && i < Nsoil) {
            frozen_soil(i, T[i], liq, ice, matric, soil_con);
        }
        else if (T[i] >= CONST_TKFRZ && i >= Nsnow && i < Nsoil) {
            if (ice[i] > 0.0) {
                liq[i] += ice[i] * CONST_RHOICE / CONST_RHOFW;
            }
            if (liq[i] > Wsat_node[i]) {
                liq[i] = Wsat_node[i];
            }
            ice[i] = 0.0;
            if (liq[i] > Wpwp_node[i]) {
                matric[i] = SoilWaterRetentionCurve(MATRIC_FLAG, i, 
                                                    liq[i], 0.0, soil_con);
            }
        }
		if (fabs(diff) > max_diff) {
			max_diff = fabs(diff);
		}
	}
    // 记录温度变化量
    energy->delt_T = max_diff;

	// 将组合温度T写回各自的温度数组中
	for (i = 0; i < Nsnow; i++) {
		pack_T[i] = T[i];
	}
	// 写回土层温度
	for (i = Nsnow; i < Nnode; i++) {
	    soil_T[i-Nsnow] = T[i];
	}
	
	return(0);
}