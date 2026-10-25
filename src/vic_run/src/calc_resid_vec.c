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
calc_resid_vec(double             step_dt,
               double            *mat_RHS,
               energy_bal_struct *energy,
               cell_data_struct  *cell,
               snow_data_struct  *snow,
               soil_con_struct   *soil_con)
{
    extern parameters_struct param;
    size_t i, j, k;
    size_t idx = 0;
    size_t lidx, nidx;
    size_t Nsoil = cell->Nsoil;
    size_t Nsnow = snow->Nsnow;
    size_t Nnode = cell->Nnode;
    double fact[MAX_SOILS] = {0};
    double FLOW[MAX_SOILS] = {0};
    double EPSLON[MAX_SOILS] = {0};
    double zc_node[MAX_SNOWS+1] = {0};
    // 指针赋值
    double *dz_soil = soil_con->dz_soil;
    double *dz_snow = snow->dz_snow;
    double *zc_soil = soil_con->zc_soil;
    double *ice = energy->ice;
    double *liq = energy->liq;
    double *T = energy->T;
    double *last_liq = energy->last_liq;
    double *last_ice = energy->last_ice;
	double *matric = cell->matric;
    double *moist = energy->moist;
    double *soil_T = cell->soil_T;
    double *last_moist = energy->last_moist;
    double *liquid_flux = cell->liquid_flux;
    double *transp_sink = cell->transp_sink;
	double *vapor_flux = cell->vapor_flux;
    double *lateral_flow = cell->lateral_flow;
    double *dQvdSMP = cell->dQvdSMP;
    double *conv_vapor = cell->conv_vapor;
    double *drhodT = cell->drhodT;
    double *compress = cell->compress;
    double *ksat_node = soil_con->Ksat_node;
    double *Wsat_node = soil_con->Wsat_node;
    double *AbsSnowLyr = energy->AbsSnowLyr;
    double *kappa_int = energy->kappa_int;
    double *enthalpy = energy->enthalpy;
    double *last_enthalpy = energy->last_enthalpy;
    double *conduct_int = cell->conduct_int;
    int ISAT[MAX_SOILS] = {0};
    size_t IFLAG = 0;
    double Qmax = 0.0;
    // 计算固定的时间步长系数
    size_t tmp_Nsnow = Nsnow;
    double capr = 0.34;  // heat capacity ratio
    for (i = 0; i < Nnode; i++) {
        if (Nsnow > 0 && i < Nsnow) {
            if (i == 0) {
                fact[i] = step_dt / 0.5 / (zc_node[i] + capr * zc_node[i+1]);
            }
            else {
                fact[i] = step_dt / dz_snow[i];
            }
        }
        else if (i == Nsnow && cell->h2osfc > param.TOL_A) {
            fact[i] = step_dt / (0.5 * cell->h2osfc);
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
    double lateral_max = 0.0;
    double coverage = snow->coverage;
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
                lateral_max = cell->infiltration + cell->soil_dew + cell->soil_frost - cell->soil_evap - 
                              cell->soil_sublim - liquid_flux[lidx] - vapor_flux[lidx] / CONST_RHOFW - 
                              transp_sink[i] / CONST_RHOFW;
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
        }
        else if (i <= Nsoil - 1) {

            if (matric[i] >= 0.0) {

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
        }
    }

    // 计算雪层和土层地表热通量
    double grnd_snow = 0.0;
    double grnd_soil = 0.0;
    double grnd_flux = 0.0;
    if (Nsnow > 0) {
        grnd_snow = energy->AbsSnowLyr[0] + energy->NetLongSnow -
                        energy->SensibleSnow - energy->LatentSnow + energy->advection;
        grnd_soil = energy->NetShortSoil + energy->NetLongSoil - 
                        energy->SensibleSoil - energy->LatentSoil + energy->advection;
        grnd_flux = grnd_snow * coverage + (1.0 - coverage) * grnd_soil;
    }
    else {
        grnd_flux = energy->shortwave + energy->longwave - 
                                energy->sensible - energy->latent + energy->advection;
        grnd_snow = grnd_flux;
        grnd_soil = grnd_flux;
    }
    energy->grnd_flux = grnd_flux;
    energy->grnd_snow = grnd_snow;
    energy->grnd_soil = grnd_soil;

    // 计算雪层和土层表面水通量
    double flux_snow = 0.0;
    double flux_soil = 0.0;
    double flux_grnd = 0.0;
    if (Nsnow > 0) {
        flux_snow = (snow->snow_dew + snow->snow_frost - 
                        snow->snow_evap - snow->snow_sublim) / CONST_RHOFW;
        flux_soil = (cell->infiltration + cell->soil_dew + cell->soil_frost - 
                    cell->soil_evap + cell->soil_sublim) / CONST_RHOFW;
        flux_grnd = flux_soil + flux_snow;
    }
    else {
        flux_grnd = (cell->infiltration + cell->soil_dew + cell->soil_frost - 
                    cell->soil_evap + cell->soil_sublim) / CONST_RHOFW;
        flux_snow = 0.0;
        flux_soil = 0.0;
    }
    // 限制最大出流
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

    double adv_left = 0.0;
    double adv_right = 0.0;
    double latent_term = 0.0;
    double trans_right = 0.0;   // 右侧界面传输系数
    double trans_left = 0.0;    // 左侧界面传输系数
    // ============================================================
    //     构建能量和水量状态向量
    // ============================================================
    for (i = 0, j = 0; i < Nnode; i++, j+=2) {
        // 单位体积内因冰量变化而释放或吸收潜热通量
        double ice_term = CONST_RHOICE * CONST_LATICE * (ice[i] - last_ice[i]);
        if (i < Nsnow) {
            // 构建右端残差项
            if (Nsnow == 1) {
                // 能量残差项
                mat_RHS[j] = enthalpy[i] - last_enthalpy[i] - coverage * (kappa_int[i] * (T[i] - T[i+1]) + 
                             CONST_LATSUB * vapor_flux[i] + grnd_snow) * fact[i] - ice_term; // [j/m3]
                // 水量残差项
                mat_RHS[j+1] = (moist[i] - last_moist[i]) - (coverage * liquid_flux[i] - liquid_flux[i+1] + 
                            (coverage * vapor_flux[i] - vapor_flux[i+1]) / CONST_RHOFW + flux_snow) * fact[i];
            }
            else if (i == 0) {
                // 能量残差项
                mat_RHS[j] = enthalpy[i] - last_enthalpy[i] - (kappa_int[i] * (T[i] - T[i+1]) -
                             CONST_LATVAP * vapor_flux[i] + grnd_snow) * fact[i] - ice_term;
                // 水量残差项
                mat_RHS[j+1] = (moist[i] - last_moist[i]) - (liquid_flux[i] - liquid_flux[i+1] + 
                               (vapor_flux[i] - vapor_flux[i+1]) / CONST_RHOFW + flux_snow) * fact[i];
            }
            else if (i < Nsnow - 1) {
                // 能量残差项
                mat_RHS[j] = enthalpy[i] - last_enthalpy[i] - (kappa_int[i-1] * (T[i-1] - T[i]) - 
                             kappa_int[i] * (T[i] - T[i+1]) + CONST_LATSUB * (vapor_flux[i] - 
                             vapor_flux[i-1]) + AbsSnowLyr[i]) * fact[i] - ice_term;
                // 水量残差项
                mat_RHS[j+1] = (moist[i] - last_moist[i]) - (liquid_flux[i] - liquid_flux[i+1] + 
                               (vapor_flux[i] - vapor_flux[i+1]) / CONST_RHOFW) * fact[i];           
            }
            else {
                // 能量残差项
                mat_RHS[j] = enthalpy[i] - last_enthalpy[i] - (kappa_int[i-1] * (T[i-1] - T[i]) - 
                             coverage * kappa_int[i] * (T[i] - T[i+1]) + CONST_LATSUB * (vapor_flux[i] * 
                             coverage - vapor_flux[i-1]) + AbsSnowLyr[i]) * fact[i] - ice_term;
                // 水量残差项
                mat_RHS[j+1] = (moist[i] - last_moist[i]) - (liquid_flux[i] - liquid_flux[i+1] + 
                               (vapor_flux[i] - vapor_flux[i+1]) / CONST_RHOFW + flux_snow) * fact[i];        
            }
        }
        else {
            lidx = i - Nsnow;
            if (i < Nnode - 1) {
                adv_right = CONST_RHOFW * CONST_CPFWICE * liquid_flux[lidx] + 
                            CONST_CPWV * vapor_flux[i];
                trans_right = kappa_int[i] + EPSLON[lidx] * adv_right;
            } 
            else {
                trans_right = kappa_int[i]; // bottom boundary
            }
            // 上一层的水分通量和冰量对热输送的影响
            if (lidx == 0) {
                trans_left = 0.0;
            }
            else if (i == Nsnow) {
                if (cell->h2osfc > param.TOL_A) {
                    trans_left = kappa_int[i-1];
                }
                else if (Nsnow > 0) {
                    trans_left = coverage * (kappa_int[Nsnow-1] + CONST_LATSUB * 
                                    conv_vapor[Nsnow-1] * drhodT[Nsnow-1]);
                }
            }
            else {
                adv_left = CONST_RHOFW * CONST_CPFWICE * liquid_flux[lidx-1] +
                                                    CONST_CPWV * vapor_flux[i-1];
                trans_left = kappa_int[i-1] + FLOW[lidx-1] * adv_left;
                latent_term = vapor_flux[i] - vapor_flux[i-1];
            }
            // 构建右端残差项
            if (i == Nsnow) {
                if (cell->h2osfc > param.TOL_A || Nsnow > 0) {
                    // 能量残差项
                    mat_RHS[j] = enthalpy[i] - last_enthalpy[i] + ((1.0 - coverage) * grnd_soil + coverage * 
                                 AbsSnowLyr[i] + coverage * trans_left * (T[i-1] - T[i]) - trans_right *
                                 (T[i] - T[i+1]) - CONST_LATVAP * vapor_flux[i]) * fact[i] - ice_term;
                    // 水量残差项
                    mat_RHS[j+1] = (moist[i] - last_moist[i]) - (liquid_flux[i] - liquid_flux[i+1] + 
                                   (vapor_flux[i] - vapor_flux[i+1]) / CONST_RHOFW + transp_sink[i] / 
                                   CONST_RHOFW + lateral_flow[lidx] + flux_snow) * fact[i] - compress[i];
                }
                else {
                    // 能量残差项
                    mat_RHS[j] = enthalpy[i] - last_enthalpy[i] + (grnd_soil - trans_right * 
                                (T[i] - T[i+1]) - CONST_LATVAP * vapor_flux[i]) * fact[i] - ice_term;
                    // 水量残差项
                    mat_RHS[j+1] = (moist[i] - last_moist[i]) - (liquid_flux[i] - liquid_flux[i+1] + 
                                   (vapor_flux[i] - vapor_flux[i+1]) / CONST_RHOFW + transp_sink[i] / 
                                   CONST_RHOFW + lateral_flow[lidx] + flux_grnd) * fact[i] - compress[i];
                }
            }
            else if (i < Nnode - 1) {
                // 能量残差项
                mat_RHS[j] = enthalpy[i] - last_enthalpy[i] - (trans_left * (T[i-1] - T[i]) - trans_right * 
                            (T[i] - T[i+1]) - CONST_LATVAP * latent_term) * fact[i] - ice_term;
                // 水量残差项
                mat_RHS[j+1] = (moist[i] - last_moist[i]) - (liquid_flux[i] - liquid_flux[i+1] + 
                                (vapor_flux[i] - vapor_flux[i+1]) / CONST_RHOFW + transp_sink[i] / 
                                CONST_RHOFW + lateral_flow[lidx]) * fact[i] - compress[i];
            }
            else {
                // 能量残差项
                mat_RHS[j] = enthalpy[i] - last_enthalpy[i] - (trans_right * (T[i] - T[i+1]) - 
                            CONST_LATVAP * vapor_flux[i]) * fact[i];
            }
        }
    }

    return (0);
}