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
    /* Initialize variables */
	size_t i;
	size_t Nsnow = snow->Nsnow;
    size_t Nsoil = cell->Nsoil;
	size_t Nnode = cell->Nnode;
    double coverage = snow->coverage;
    double *fact = energy->fact; 
	double mat_A[MAX_NODES];
	double mat_B[MAX_NODES];
	double mat_C[MAX_NODES];
	double mat_RHS[MAX_NODES];
    double FLOW[MAX_SOILS];
    double EPSLON[MAX_SOILS];
	double *T = energy->T;
	double *soil_T = cell->soil_T;
	double *pack_T = snow->pack_T;
	double *dz_node = cell->dz_node;
    double *zc_node = cell->zc_node;
    double *Zsum_node = cell->Zsum_node;
	double *kappa_int = energy->kappa_int;
    double *Cs_node = energy->Cs_node;
    double *last_Cs = energy->last_Cs;
    double *last_T = energy->last_T;
    double *matric = cell->matric;
    double *Wpwp_node = soil_con->Wpwp_node;
    double *last_matric = cell->last_matric;
    double deriv_terms = energy->deriv_terms;
	/* initialization */
    for (i = 0; i < MAX_NODES; i++) {
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
    double *conduct_int = cell->conduct_int;
    double *Wsat_node = soil_con->Wsat_node;
    for (i = Nsoil - 2; i > 0; i--) {
        if (ice[i] > 0.0) {
            Qmax = (Wsat_node[i] - liq[i] - ice[i]) * dz_node[i] / step_dt;
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
        Qmax = (Wsat_node[0] - liq[0] - ice[0]) * dz_node[0] / step_dt;
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
                Qmax = (Wsat_node[i] - liq[i] - ice[i]) * dz_node[i] / step_dt;
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
    for (i = 0; i < Nnode; i++) {
		if (i == 0) {
			fact[i] = (0.5 * (zc_node[i] - Zsum_node[i] + 
						capr * (zc_node[i+1] - Zsum_node[i]))) / step_dt;
		}
		else if (i <= Nnode - 1) {
			fact[i] = dz_node[i] / step_dt;
		}
    }

    double grnd_flux = energy->NetShortGrnd + energy->longwave - 
                        energy->sensible - energy->latent + energy->advection;
    energy->grnd_flux = grnd_flux;
    for (i = 0; i < Nnode; i++) {
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
            ice_term = CONST_RHOICE * CONST_LATICE *
                    (ice[i] - last_ice[i]);
        }
        if (i == Nnode - 1) {
            // 基岩层：不透水，无水分活动
            trans_left = kappa_int[i-1];  // 只有纯热传导，无水流通量贡献
            trans_right = 0.0;
        }
        if (Nsnow == 0 && i == 0) {
            mat_A[i] = 0.0;
            mat_B[i] = deriv_terms - trans_right - fact[i] * Cs_node[i];
            mat_C[i] = trans_right;
            mat_RHS[i] = grnd_flux - trans_right * (T[i] - T[i+1]) - CONST_LATICE *
                    vapor_flux[i] - fact[i] * (last_Cs[i] * (T[i] - last_T[i]) - ice_term);
        }
        else if (i == Nsnow) {
            mat_A[i] = coverage * trans_left;
            mat_B[i] = deriv_terms - coverage * trans_left - trans_right - fact[i] * Cs_node[i];
            mat_C[i] = trans_right;
            mat_RHS[i] = grnd_flux + coverage * trans_left * (T[i-1] - T[i]) - trans_right * (T[i] - T[i+1]) -
                         CONST_LATICE * latent_term - fact[i] * (last_Cs[i] * (T[i] - last_T[i]) - ice_term);
        }
        else if (i <= Nnode - 2) {
            mat_A[i] = trans_left;
            mat_B[i] = -trans_left - trans_right - fact[i] * Cs_node[i];
            mat_C[i] = trans_right;
            mat_RHS[i] = trans_left * (T[i-1] - T[i]) - trans_right * (T[i] - T[i+1]) -
                         CONST_LATICE * latent_term - fact[i] * (last_Cs[i] * (T[i] - last_T[i]) - ice_term);
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