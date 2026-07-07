/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate glacier accumulation and melt
 *****************************************************************************/

#include "vic_run.h"

/*****************************************************************************
  * @brief   Calculate glacier accumulation and melt using an energy balance
  *          approach for a two layer model
 *****************************************************************************/
int
calc_water_bal(double             step_dt,
               double             pressure,
               energy_bal_struct *energy,
               cell_data_struct  *cell,
               snow_data_struct  *snow,
               soil_con_struct   *soil_con)
{
    extern parameters_struct param;
    size_t i, j, k;
    size_t idx = 0;
    size_t lidx, nidx;
    size_t IFLAG = 0;
    size_t Nsoil = cell->Nsoil;
    double fact[MAX_SOILS];
    double mat_A[MAX_SOILS];
    double mat_B[MAX_SOILS];
    double mat_C[MAX_SOILS];
    double mat_RHS[MAX_SOILS];
    // 指针赋值
    double *dz_soil = soil_con->dz_soil;
    double *ice = cell->ice;
    double *liq = cell->liq;
    double *last_liq = cell->last_liq;
    double *last_ice = cell->last_ice;
	double *matric = cell->matric;
    double *zc_soil = soil_con->zc_soil;
    double *liquid_flux = cell->liquid_flux;
    double *transp_sink = cell->transp_sink;
	double *vapor_flux = cell->vapor_flux;
    double *lateral_flow = cell->lateral_flow;
    double *deriv_vapor = cell->deriv_vapor;
    double *ksat_node = soil_con->Ksat_node;
    double *Wpwp_node = soil_con->Wpwp_node;
    double *Wsat_node = soil_con->Wsat_node;
    double *conduct_int = cell->conduct_int;
    // 计算地表水分通量限制
    double Qmax = 0.0;
    int ISAT[MAX_SOILS];
    for (i = 0; i < MAX_SOILS; i++) {
        ISAT[i] = 0;
        fact[i] = 0.0;
        mat_A[i] = 0.0;
        mat_B[i] = 0.0;
        mat_C[i] = 0.0;
        mat_RHS[i] = 0.0;
    }

    // Compute the vapor flux between nodes
    calc_vapor_flux(pressure, cell, snow, soil_con);

    // Compute the hydraulic conductivity between nodes
    soil_hydraulic_conductivity(cell, soil_con);

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
    size_t tmp_Nsnow = snow->Nsnow;
    if (cell->h2osfc > param.TOL_A) {
        tmp_Nsnow++;
    }
    // 计算固定的时间步长系数
    for (i = 0; i < Nsoil; i++) {
		fact[i] = dz_soil[i] / step_dt;
    }
    double lateral_max = 0.0;
    double seepage = 0.0;
    double deric_matric = 0.0;
    double tmp_matric = 0.0;
    double tmp_liq = 0.0;
    double esoil = cell->esoil;
    double deriv_evap = energy->deriv_evap;
    // 确定表层的系数
    for (i = 0; i < Nsoil; i++) {
        lidx = tmp_Nsnow + i;
        if (i == 0) {
            if (matric[i] >= 0.0) {
                // 允许侧向流出剖面
                lateral_flow[i] = ksat_node[i] * sin(soil_con->slope * 
                                (CONST_PI / 180.0)) * dz_soil[i];
                
                // 限制侧向流，使得该层不会仅因侧向流而脱饱和
                lateral_max = mat_RHS[i] - liquid_flux[i] - vapor_flux[lidx] / 
                                CONST_RHOFW - transp_sink[i] / CONST_RHOFW;
                if (lateral_max < 0.0) {
                    lateral_max = 0.0;
                }
                if (lateral_flow[i] > lateral_max) {
                    lateral_flow[i] = lateral_max;
                }
            } 
            else {
                lateral_flow[i] = 0.0;
            }
            mat_RHS[i] = esoil / CONST_RHOFW - liquid_flux[i] - vapor_flux[lidx] / CONST_RHOFW -
                         transp_sink[i] / CONST_RHOFW - lateral_flow[i] - fact[i] * (liq[i] - 
                         last_liq[i] + CONST_RHOICE / CONST_RHOFW * (ice[i] - last_ice[i]));

            // 存在冰 - 计算冰含量的系数
            if (ice[i] > 0.0) {
                mat_A[i] = 0.0;
                mat_B[i] = -fact[i] * CONST_RHOICE / CONST_RHOFW;
                mat_C[i] = 0.0;
                seepage = 0.0;
            } 
            else {
                // 无冰存在 - 计算液态水平衡系数
                if (matric[i] < 0.0 || mat_RHS[i] < 0.0) {
                    // 无渗漏
                    deric_matric = 
                            SoilWaterRetentionCurve(DERIV_FLAG, i,
                                                    liq[i], matric[i], 
                                                    soil_con);

                    if (matric[i] >= 0.0) {
                        // 节点饱和 - 检查边界是否定义不良
                        if (fabs(mat_B[i]) < conduct_int[i] * param.TOL_A || conduct_int[i] <= 0.0) {
                            // 边界定义不良 - 将节点视为饱和
                            ISAT[idx] = i;
                            idx++;
                        }
                    }
                    mat_A[i] = 0.0;
                    mat_B[i] = deriv_evap - (conduct_int[i] + deriv_vapor[lidx] / 
                                            CONST_RHOFW) - fact[i] * deric_matric;
                    mat_C[i] = conduct_int[i] + deriv_vapor[lidx] / CONST_RHOFW;
                    seepage = 0.0;
                } 
                else {
                    // 允许表层渗漏
                    mat_A[i] = conduct_int[i];
                    mat_B[i] = 1.0;
                    mat_C[i] = conduct_int[i] + deriv_vapor[lidx] / CONST_RHOFW;
                    // 流入节点的净流量等于表层渗漏量
                    seepage = mat_RHS[i];
                    mat_RHS[i] = 0.0;
                }
            }
        }
        else if (i <= Nsoil - 2) {
            if (ice[i] > 0.0) {
                // 存在冰 - 计算冰含量的系数
                mat_A[i] = 0.0;
                mat_B[i] = -fact[i] * CONST_RHOICE / CONST_RHOFW;
                mat_C[i] = 0.0;
                lateral_flow[i] = 0.0;
            } 
            else {
                // 无冰存在 - 计算液态水平衡系数
                deric_matric = 
                        SoilWaterRetentionCurve(DERIV_FLAG, i,
                                                liq[i], matric[i], 
                                                soil_con);
                
                if (matric[i] >= 0.0) {
                    // 节点饱和
                    ISAT[idx] = i;
                    idx++;
                    
                    // 计算侧向流出剖面
                    lateral_flow[i] = ksat_node[i] * sin(soil_con->slope * 
                                        (CONST_PI / 180.0)) * dz_soil[i];
                    
                    // 限制侧向流，使该层不会脱饱和
                    lateral_max = liquid_flux[i-1] - liquid_flux[i] +
                            (vapor_flux[lidx-1] - vapor_flux[lidx]) / CONST_RHOFW - transp_sink[i] / CONST_RHOFW -
                            fact[i] * (liq[i] - last_liq[i] + CONST_RHOICE / CONST_RHOFW * (ice[i] - last_ice[i]));
                    
                    if (lateral_max < 0.0) {
                        lateral_max = 0.0;
                    }
                    if (lateral_flow[i] > lateral_max) {
                        lateral_flow[i] = lateral_max;
                    }
                } 
                else {
                    lateral_flow[i] = 0.0;
                }
                mat_A[i] = conduct_int[i-1] + deriv_vapor[lidx-1] / CONST_RHOFW;
                mat_B[i] = -(conduct_int[i-1] + conduct_int[i] + (deriv_vapor[lidx-1] + 
                             deriv_vapor[lidx]) / CONST_RHOFW) - fact[i] * deric_matric;
                mat_C[i] = conduct_int[i] + deriv_vapor[lidx] / CONST_RHOFW;
            }
            
            mat_RHS[i] = liquid_flux[i-1] - liquid_flux[i] + (vapor_flux[lidx-1] - vapor_flux[lidx]) / 
                         CONST_RHOFW - transp_sink[i] / CONST_RHOFW - lateral_flow[i] -
                         fact[i] * (liq[i] - last_liq[i] + CONST_RHOICE / CONST_RHOFW * (ice[i] - last_ice[i]));
        }
        else if (i == Nsoil - 1) {
            if (ice[i] > 0.0) {
                // 存在冰 - 计算冰含量的系数
                mat_A[i] = 0.0;
                mat_B[i] = -fact[i] * CONST_RHOICE / CONST_RHOFW;
                mat_C[i] = 0.0;
                lateral_flow[i] = 0.0;
            } 
            else {
                // 无冰存在 - 计算液态水平衡系数
                deric_matric = 
                        SoilWaterRetentionCurve(DERIV_FLAG, i,
                                                liq[i], matric[i], 
                                                soil_con);
                
                if (matric[i] >= 0.0) {
                    // 节点饱和
                    ISAT[idx] = i;
                    idx++;
                    
                    // 计算侧向流出剖面
                    lateral_flow[i] = ksat_node[i] * sin(soil_con->slope * 
                                        (CONST_PI / 180.0)) * dz_soil[i];
                    
                    // 限制侧向流，使该层不会脱饱和
                    lateral_max = liquid_flux[i-1] - liquid_flux[i] +
                            (vapor_flux[lidx-1] - vapor_flux[lidx]) / CONST_RHOFW - transp_sink[i] / CONST_RHOFW -
                            fact[i] * (liq[i] - last_liq[i] + CONST_RHOICE / CONST_RHOFW * (ice[i] - last_ice[i]));
                    
                    if (lateral_max < 0.0) {
                        lateral_max = 0.0;
                    }
                    if (lateral_flow[i] > lateral_max) {
                        lateral_flow[i] = lateral_max;
                    }
                } 
                else {
                    lateral_flow[i] = 0.0;
                }
                mat_A[i] = conduct_int[i-1] + deriv_vapor[lidx-1] / CONST_RHOFW;
                mat_B[i] = -(conduct_int[i-1] + conduct_int[i] + (deriv_vapor[lidx-1] + 
                            deriv_vapor[lidx]) / CONST_RHOFW) - fact[i] * deric_matric;
                mat_C[i] = 0.0;
            }
            
            mat_RHS[i] = liquid_flux[i-1] + vapor_flux[lidx-1] / CONST_RHOFW - transp_sink[i] / 
                         CONST_RHOFW - lateral_flow[i] - fact[i] * (liq[i] - last_liq[i] + 
                         CONST_RHOICE / CONST_RHOFW * (ice[i] - last_ice[i]));
        }
    }

    // 如果发生渗漏，调整矩阵
    if (seepage > 0.0) {
        mat_A[0] = 0.0;
        mat_C[0] = 0.0;
    }
    
    if (idx > 0) {
        // 存在饱和节点 - 检查矩阵奇异性
        int jtop = -1;
        ISAT[idx] = Nsoil;
        
        for (j = 0; j < idx; j++) {

            k = ISAT[j];
            // 检查矩阵不连续性并保存顶部节点
            if (ISAT[j] == 0) {
                jtop = ISAT[j];
            } 
            else {
                // 检查上方节点的导水率是否足够小，可视为零
                if (conduct_int[ISAT[j]-1] <= fabs(mat_B[k]) * 1.0e-7) {
                    jtop = ISAT[j];
                }
                // 检查饱和区上方非饱和节点的比储水是否可忽略
                if (-conduct_int[ISAT[j]-1] - mat_B[k-1] <= fabs(mat_B[k]) * 1.0e-7 && jtop == -1) {
                    jtop = ISAT[j];
                }
            }
            // 如果找到了饱和区顶部，检查是否需要修正
            if (jtop != -1) {
                if (conduct_int[ISAT[j]] <= fabs(mat_B[k]) * 1.0e-7) {
                    // 饱和节点被零导水率包围；解未定义
                    // 重置顶部饱和层的mat_B，如同它刚好在进气值以下
                    nidx = jtop;
                    tmp_liq = Wsat_node[nidx] - 0.001;
                    tmp_matric = 
                            SoilWaterRetentionCurve(MATRIC_FLAG, nidx,
                                                    tmp_liq, 0.0, 
                                                    soil_con);
                    deric_matric = 
                            SoilWaterRetentionCurve(DERIV_FLAG, nidx,
                                                    tmp_liq,
                                                    tmp_matric, 
                                                    soil_con);
                    
                    mat_B[k] = mat_B[k] - fact[k] * deric_matric;
                    // 重置顶部饱和层
                    jtop = -1;
                }
            }
            // 如果下一个节点不饱和，重置nidx
            if (ISAT[j] + 1 != ISAT[j+1]) {
                jtop = -1;
            }
        }
    }

    // 调用LAPACK函数求解线性方程组
    int info = LAPACKE_dgtsv(LAPACK_COL_MAJOR, Nsoil, 
                                1, mat_A+1, mat_B, mat_C, mat_RHS, Nsoil);

    if (info != 0) {
        return ERROR;
    }

    /* 检查收敛 */
    double max_diff = 0.0;
	for (i = 0; i < Nsoil; i++) {
		double diff = mat_RHS[i];
        if (ice[i] > 0.0) {
            if (diff > 0.2) {
                double energy_flux = fabs(CONST_RHOICE * CONST_LATICE * fact[i] * diff);
                if (energy_flux > 2000.0) {
                    mat_RHS[i] = diff / fabs(diff) * 2000.0 * step_dt / 
                    (CONST_RHOICE * CONST_LATICE * dz_soil[i]);
                }
            }
            ice[i] -= mat_RHS[i];
            if (ice[i] < 0.0) {
                double total_liq = liq[i] + ice[i] * CONST_RHOICE / CONST_RHOFW;
                double tol_liq = (liq[i] + Wpwp_node[i]) / 2.0;
                if (total_liq > tol_liq) {
                    liq[i] = total_liq;
                }
                else {
                    liq[i] = tol_liq;
                }
                if (liq[i] > Wpwp_node[i]) {
                    matric[i] = SoilWaterRetentionCurve(MATRIC_FLAG, i, 
                                                        liq[i], 0.0, soil_con);
                }
                ice[i] = 0.0;
            }
        }
        else {
            double tmp_mat = matric[i];
            if (fabs(tmp_mat) < 1.0) {
                tmp_mat = 1.0;
            }
            // 计算基质势的相对变化量
            diff /= tmp_mat;
            matric[i] -= diff;
            liq[i] = SoilWaterRetentionCurve(MOIST_FLAG, i,
                                             0.0, matric[i], soil_con);
        }
        if (fabs(diff) > max_diff) {
            max_diff = fabs(diff);
        }
	}
    if (matric[0] > 0.0) {
        double excess = matric[0];  // 超出的水头
        for (j = 1; j < Nsoil; j++) {
            if (matric[j] > 0.0) {
                matric[j] -= excess;
                liq[j] = SoilWaterRetentionCurve(MOIST_FLAG, j,
                                                 0.0, matric[j], soil_con);
            }
            else {
                break;
            }
        }
        // 将地表基质势强制设为 0（自由水面）
        matric[0] = 0.0;
        liq[0] = SoilWaterRetentionCurve(MOIST_FLAG, 0,
                                         0.0, matric[0], soil_con);
    }
    // 判断水量收敛标志
    if (max_diff < 0.01) {
        energy->moist_flag = true;
    }
    else {
        energy->moist_flag = false;
    }
    
    return (0);
}
