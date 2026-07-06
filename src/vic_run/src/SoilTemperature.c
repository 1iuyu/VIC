/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate snow and soil layer temperature.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief  
 * Compute bare ground soil and snow layer temperatures.
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
	double *kappa_int = energy->kappa_int;
    double *Cs_node = energy->Cs_node;
    double *last_Cs = energy->last_Cs;
    double *last_T = energy->last_T;
    double *matric = cell->matric;
    double *pack_liq = snow->pack_liq;
    double *pack_ice = snow->pack_ice;
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
    double *theta_liq = snow->theta_liq;
    double *last_packliq = snow->last_packliq;
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
                } 
                else {
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
    size_t tmp_Nsnow = Nsnow;
    double capr = 0.34;  // heat capacity ratio
    for (i = 0; i < Nnode; i++) {
        if (i < Nsnow) {
            fact[i] = dz_snow[i] / step_dt;
        }
        else if (i == Nsnow && cell->h2osfc > param.TOL_A) {
            fact[i] = 0.5 * cell->h2osfc / step_dt;
            tmp_Nsnow++;
        }
        else {
            if (i == 0) {
                fact[i] = 0.5 * (zc_soil[i] + capr * zc_soil[i+1]) / step_dt;
            }
            else {
                lidx = i - tmp_Nsnow;
                fact[i] = dz_soil[lidx] / step_dt;
            }
        }
    }
    // 计算雪层和土层地表热通量
    double grnd_flux = energy->shortwave + energy->longwave - 
                            energy->sensible - energy->latent + energy->advection;
    double grnd_snow = energy->NetShortSnow + energy->NetLongSnow -
                    energy->SensibleSnow - energy->LatentSnow + energy->advection;
    double grnd_soil = energy->NetShortSoil + energy->NetLongSoil - 
                    energy->SensibleSoil - energy->LatentSoil + energy->advection;
    energy->grnd_flux = grnd_flux;

    // ============================================================
    //     构建能量平衡三对角矩阵
    // ============================================================
    for (i = 0; i < Nnode; i++) {
        if (i < Nsnow) {
            if (i == 0) {
                if (pack_liq[i] > 0.0) {
                    mat_A[i] = 0.0;
                    mat_B[i] = deriv_terms - CONST_RHOFW * CONST_LATICE / step_dt;
                    mat_C[i] = 0.0;
                    mat_RHS[i] = grnd_snow - CONST_RHOFW * CONST_LATICE * 
                                 (theta_liq[i] - last_packliq[i]) / step_dt;
                }
                else {
                    mat_A[i] = 0.0;
                    mat_B[i] = deriv_terms - kappa_int[i] - CONST_LATSUB * deric_vapor[i] * 
                            conv_vapor[i] - fact[i] * Cs_node[i];
                    mat_C[i] = kappa_int[i] + CONST_LATSUB * deric_vapor[i+1] * conv_vapor[i];
                    mat_RHS[i] = grnd_snow - kappa_int[i] * (T[i] - T[i+1]) - fact[i] * 
                            last_Cs[i] * (T[i] - last_T[i]) - CONST_RHOFW * CONST_LATICE * 
                            (theta_liq[i] - last_packliq[i]) / step_dt - CONST_LATSUB * vapor_flux[i];
                }
            }
            else if (i < Nsnow - 1) {
                if (pack_liq[i] > 0.0) {
                    mat_A[i] = 0.0;
                    mat_B[i] = -CONST_RHOFW * CONST_LATICE / step_dt;
                    mat_C[i] = 0.0;
                    mat_RHS[i] = -CONST_RHOFW * CONST_LATICE * (theta_liq[i] - last_packliq[i]) / step_dt;
                }
                else {
                    mat_A[i] = kappa_int[i-1] + CONST_LATSUB * deric_vapor[i-1] * conv_vapor[i-1];
                    mat_B[i] = -(kappa_int[i-1] + kappa_int[i]) - CONST_LATSUB * deric_vapor[i] *
                                    (conv_vapor[i-1] + conv_vapor[i]) - fact[i] * Cs_node[i];
                    mat_C[i] = kappa_int[i] + CONST_LATSUB * deric_vapor[i+1] * conv_vapor[i];
                    mat_RHS[i] = kappa_int[i-1] * (T[i-1] - T[i]) - kappa_int[i] * (T[i] - T[i+1]) - fact[i] * 
                                 last_Cs[i] * (T[i] - last_T[i]) - CONST_RHOFW * CONST_LATICE * (theta_liq[i] - 
                                 last_packliq[i]) / step_dt - CONST_LATSUB * (vapor_flux[i] - vapor_flux[i-1]);
                }
            }
            else {
                if (pack_liq[i] > 0.0) {
                    mat_A[i] = 0.0;
                    mat_B[i] = -CONST_RHOFW * CONST_LATICE / step_dt;
                    mat_C[i] = 0.0;
                    mat_RHS[i] = -CONST_RHOFW * CONST_LATICE * (theta_liq[i] - last_packliq[i]) / step_dt;
                } 
                else {
                    mat_A[i] = kappa_int[i-1] + CONST_LATSUB * deric_vapor[i-1] * conv_vapor[i-1];
                    mat_B[i] = -(kappa_int[i-1] + kappa_int[i]) - CONST_LATSUB * deric_vapor[i] *
                            conv_vapor[i-1] - fact[i] * Cs_node[i];
                    mat_C[i] = kappa_int[i];
                    mat_RHS[i] = kappa_int[i-1] * (T[i-1] - T[i]) - kappa_int[i] * (T[i] - T[i+1]) - fact[i] * 
                                 last_Cs[i] * (T[i] - last_T[i]) - CONST_RHOFW * CONST_LATICE * (theta_liq[i] - 
                                 last_packliq[i]) / step_dt - CONST_LATSUB * (vapor_flux[i] - vapor_flux[i-1]);
                }
            }
        }
        else if (i == Nsnow && cell->h2osfc > param.TOL_A) {
            if (Nsnow == 0) {
                if (cell->h2osfc_liq == 0.0) {
                    mat_A[i] = 0.0;
                    mat_B[i] = deriv_terms - kappa_int[i] - fact[i] * Cs_node[i];
                    mat_C[i] = kappa_int[i];
                    mat_RHS[i] = grnd_soil - kappa_int[i] * (T[i] - T[i+1]) - 
                                fact[i] * last_Cs[i] * (T[i] - last_T[i]);
                } else {
                    mat_A[i] = 0.0;
                    mat_B[i] = deriv_terms - CONST_RHOFW * CONST_LATICE / step_dt - fact[i] * Cs_node[i];
                    mat_C[i] = 0.0;
                    mat_RHS[i] = grnd_soil - CONST_RHOFW * CONST_LATICE * 
                                (theta_liq[i-1] - last_packliq[i-1]) / step_dt -
                                fact[i] * last_Cs[i] * (T[i] - last_T[i]);
                }
            } 
            else {
                if (cell->h2osfc_liq == 0.0) {
                    mat_A[i] = kappa_int[i-1];
                    mat_B[i] = -(kappa_int[i-1] + kappa_int[i]) - fact[i] * Cs_node[i];
                    mat_C[i] = kappa_int[i];
                    mat_RHS[i] = kappa_int[i-1] * (T[i-1] - T[i]) - kappa_int[i] * 
                            (T[i] - T[i+1]) - fact[i] * last_Cs[i] * (T[i] - last_T[i]);
                } else {
                    mat_A[i] = 0.0;
                    mat_B[i] = -CONST_RHOFW * CONST_LATICE / step_dt - fact[i] * Cs_node[i];
                    mat_C[i] = 0.0;
                    mat_RHS[i] = -CONST_RHOFW * CONST_LATICE * (theta_liq[i-1] - last_packliq[i-1]) / 
                                 step_dt - fact[i] * last_Cs[i] * (T[i] - last_T[i]);
                }
            }
        }
        else {
            lidx = i - tmp_Nsnow;
            double adv_left = 0.0;
            double adv_right = 0.0;
            double latent_term = 0.0;
            double ice_term = 0.0;
            double trans_right = 0.0;   // 右侧界面传输系数
            double trans_left = 0.0;    // 左侧界面传输系数
            // 当前层的水分通量和冰量对热输送的影响
            if (i < Nnode - 2) {
                adv_right = CONST_RHOFW * CONST_CPFWICE * liquid_flux[lidx] + 
                            CONST_CPWV * vapor_flux[i];
                trans_right = kappa_int[i] + EPSLON[lidx] * adv_right;
            }
            else if (i == Nnode - 2) {
                trans_right = kappa_int[i];
            }
            // 上一层的水分通量和冰量对热输送的影响
            if (lidx == 0) {
                if (cell->h2osfc > param.TOL_A) {
                    trans_left = kappa_int[i-1];
                }
                else if (Nsnow > 0) {
                    trans_left = coverage * (kappa_int[Nsnow-1] + CONST_LATVAP * 
                                    conv_vapor[Nsnow-1] * deric_vapor[Nsnow-1]);
                }
            }
            else if (i < Nnode - 1) {
                adv_left = CONST_RHOFW * CONST_CPFWICE * liquid_flux[lidx-1] +
                                                    CONST_CPWV * vapor_flux[i-1];
                trans_left = kappa_int[i-1] + FLOW[lidx-1] * adv_left;
                latent_term = vapor_flux[i] - vapor_flux[i-1];
            }
            if (i <= Nnode - 2) {
                ice_term = CONST_RHOICE * CONST_LATICE * (ice[lidx] - last_ice[lidx]);
            }
            if (i == 0) {
                if (cell->h2osfc > param.TOL_A || Nsnow > 0) {
                    mat_A[i] = trans_left;
                    mat_B[i] = -trans_left - trans_right - fact[i] * Cs_node[i];
                    mat_C[i] = trans_right;
                    mat_RHS[i] = trans_left * (T[i-1] - T[i]) - trans_right * 
                            (T[i] - T[i+1]) - CONST_LATVAP * vapor_flux[i] - 
                            fact[i] * (last_Cs[i] * (T[i] - last_T[i]) - ice_term);
                }
                else {
                    mat_A[i] = 0.0; // 表土层裸露于大气
                    mat_B[i] = deriv_terms - trans_right - fact[i] * Cs_node[i];
                    mat_C[i] = trans_right;
                    mat_RHS[i] = grnd_soil - trans_right * (T[i] - T[i+1]) - 
                                 CONST_LATVAP * vapor_flux[i] - fact[i] * 
                                 (last_Cs[i] * (T[i] - last_T[i]) - ice_term);
                }
            }
            else if (i <= Nnode - 2) {
                mat_A[i] = trans_left;
                mat_B[i] = -trans_left - trans_right - fact[i] * Cs_node[i];
                mat_C[i] = trans_right;
                mat_RHS[i] = trans_left * (T[i-1] - T[i]) - trans_right * 
                            (T[i] - T[i+1]) - CONST_LATVAP * latent_term - 
                            fact[i] * (last_Cs[i] * (T[i] - last_T[i]) - ice_term);
            }
            else {
                mat_A[i] = kappa_int[i-1];
                mat_B[i] = -kappa_int[i-1] - fact[i] * Cs_node[i];
                mat_C[i] = 0.0;
                mat_RHS[i] = kappa_int[i-1] * (T[i-1] - T[i]) - fact[i] * 
                                                (last_Cs[i] * (T[i] - last_T[i]));
            }
        }
    }
    // 修正能量平衡方程的系数
    double liq_deriv1 = 0.0;
    double liq_deriv2 = 0.0;
    double mat_deriv = 0.0;
    double tmp_matric = 0.0;
    for (i = 0; i < Nsoil; i++) {
        lidx = tmp_Nsnow + i;
        double tmp_tkfrz = CONST_TKTRIP;
        if (ice[i] > 0.0) {
            double total_liq = liq[i] + ice[i] * CONST_RHOICE / CONST_RHOFW;
            if (total_liq > Wsat_node[i]) {
                total_liq = Wsat_node[i];
            }
            double tmp_mat = SoilWaterRetentionCurve(MATRIC_FLAG, i,
                                                     total_liq, 0.0, soil_con);

            if (tmp_mat < 0.0) {
                tmp_tkfrz = CONST_TKTRIP * tmp_mat /
                            (CONST_LATICE / CONST_G - tmp_mat) + CONST_TKTRIP;
            }
        }
        if (matric[i] >= 0.0) {
            tmp_matric = CONST_LATICE * (T[lidx] - tmp_tkfrz) / T[lidx] / CONST_G;
            liq_deriv2 = water_curve_deriv(i, T[lidx], 
                                            liq[i], 
                                            tmp_matric,
                                            soil_con);
            mat_deriv = SoilWaterRetentionCurve(DERIV_FLAG, i, liq[i],
                                                tmp_matric, soil_con);
        }
        else {
            liq_deriv2 = water_curve_deriv(i, T[lidx], liq[i], 
                                            matric[i], 
                                            soil_con);
            mat_deriv = SoilWaterRetentionCurve(DERIV_FLAG, i, liq[i],
                                                matric[i], soil_con);
        }
        if (last_matric[i] < 0.0 || last_ice[i] > 0.0) {
            liq_deriv1 = water_curve_deriv(i, last_T[lidx],
                                            last_liq[i], 
                                            last_matric[i],
                                            soil_con);
        }
        else if (ice[i] > 0.0) {
            liq_deriv1 = liq_deriv2;
        }
        if ((liq_deriv1 > 0.0 || liq_deriv2 > 0.0) && (mat_deriv > 0.0)) {
            if (i == 0) {
                mat_B[lidx] -= fact[lidx] * 0.5 * CONST_RHOFW * CONST_LATICE * (liq_deriv1 + liq_deriv2) +
                        CONST_RHOFW * CONST_LATICE * conduct_int[i] * liq_deriv2 / mat_deriv;
            }
            else if (i < Nsoil) {
                mat_B[lidx] -= fact[lidx] * 0.5 * CONST_RHOFW * CONST_LATICE * (liq_deriv1 + liq_deriv2) +
                    CONST_RHOFW * CONST_LATICE * (conduct_int[i-1] + conduct_int[i]) * liq_deriv2 / mat_deriv;
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
        if (i < Nsnow) {
            if (pack_liq[i] == 0.0) {
                T[i] -= diff;
                if (T[i] > CONST_TKFRZ) {
                    pack_liq[i] = CONST_RHOICE * CONST_CPICE * dz_snow[i] * 
                                (T[i] - CONST_TKFRZ) / (CONST_RHOFW * CONST_LATICE);
                    pack_ice[i] -= pack_liq[i];
                    T[i] = CONST_TKFRZ;
                }
            }
            else {
                double tmp_liq = pack_liq[i];
                pack_liq[i] -= diff;
                if (pack_liq[i] < 0.0) {
                    T[i] = CONST_TKFRZ + pack_liq[i] * CONST_RHOFW * CONST_LATICE / 
                                    (CONST_RHOICE * CONST_CPICE * dz_snow[i]);
                    pack_ice[i] += tmp_liq;
                    pack_liq[i] = 0.0;
                }
                else {
                    pack_ice[i] += diff;
                    T[i] = CONST_TKFRZ;
                }
            }
        }
        else if (i == Nsnow && cell->h2osfc > param.TOL_A) {
            if (cell->h2osfc > param.TOL_A) {
                if (cell->h2osfc_liq == 0.0) {
                    T[i] -= diff;
                    if (T[i] > CONST_TKFRZ) {
                        cell->h2osfc_liq = CONST_RHOICE * CONST_CPICE * cell->h2osfc * 
                                        (T[i] - CONST_TKFRZ) / (CONST_RHOFW * CONST_LATICE);
                        cell->h2osfc_ice -= cell->h2osfc_liq;
                        T[i] = CONST_TKFRZ;
                    }
                }
                else {
                    double tmp_liq = cell->h2osfc_liq;
                    cell->h2osfc_liq -= diff;
                    if (cell->h2osfc_liq < 0.0) {
                        T[i] = CONST_TKFRZ + cell->h2osfc_liq * CONST_RHOFW * CONST_LATICE / 
                                        (CONST_RHOICE * CONST_CPICE * cell->h2osfc);
                        cell->h2osfc_ice += tmp_liq;
                        cell->h2osfc_liq = 0.0;
                    }
                    else {
                        cell->h2osfc_ice += diff;
                        T[i] = CONST_TKFRZ;
                    }
                }
            }
        }
        else if (i < Nnode - 1) {
            T[i] -= diff;
            lidx = i - tmp_Nsnow;
            // 判断是否需要处理相变
            CalcPhaseChange(lidx, energy,
                            cell, soil_con);
        }
        if (fabs(diff) > max_diff) {
            max_diff = fabs(diff);
        }
	}
    // 判断能量收敛标志
    if (max_diff < 0.01) {
        energy->energy_flag = true;
    }
    else {
        energy->energy_flag = false;
    }

	// 将组合温度T写回各自的温度数组中
	for (i = 0; i < Nnode; i++) {
        if (i < Nsnow) {
            pack_T[i] = T[i];
        }
        else if (i == Nsnow && cell->h2osfc > param.TOL_A) {
            cell->h2osfc_T = T[Nsnow];
        }
        else {
            lidx = i - tmp_Nsnow;
            soil_T[lidx] = T[i];
        }	
	}
	
	return(0);
}