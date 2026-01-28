/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate the phase change (melting/freezing) of snow water and soil water.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    Calculate the aerodynamic resistance for each vegetation layer,
 *           based on Monin-Obukhov (M-O) Similarity Theory (MOST).
 *****************************************************************************/
int
CalcPhaseChange(double             step_dt,
                energy_bal_struct *energy,
                cell_data_struct  *cell,
                snow_data_struct  *snow,
                soil_con_struct   *soil_con)
{
    extern option_struct     options;

    size_t  i, j, nidx;
    double  RestTerm;
    double  old_swq;
    double  ratio;
    double *T = energy->T;
    double *liq = cell->liq;
    double *moist = cell->moist;
    double *pack_T = snow->pack_T;
    double *soil_T = cell->soil_T;
    double *pack_ice = snow->pack_ice;
    double *pack_liq = snow->pack_liq;
    double *dz_node = cell->dz_node;          // layer and thermal node thickness (m)
    int *PhaseChange = cell->PhaseChange;         // 0 = none, 1 = melt, 2 = freeze
    double *fusion_fact = energy->fusion_fact;
    double *fusion_flux = energy->fusion_flux;
    
    double  init_ice[MAX_NODES + MAX_SNOWS];
    double  init_liq[MAX_NODES + MAX_SNOWS];
    double  mass_ice[MAX_NODES + MAX_SNOWS];
    double  mass_liq[MAX_NODES + MAX_SNOWS];
    double  EnergyRes[MAX_NODES + MAX_SNOWS];
    double  supercool[MAX_NODES + MAX_SNOWS];
    double  init_moist[MAX_NODES + MAX_SNOWS];

    // 数组大小计算
    size_t Nnode = options.Nnode;
    size_t Nsnow = snow->Nsnow;
    size_t Total_Layer = MAX_NODES + MAX_SNOWS;
    double melt_grnd = 0.0;
    double latent_fusion = 0.0;
    
    // 初始化数组
    for (i = 0; i < Total_Layer; i++) {
        mass_ice[i] = 0.0;
        mass_liq[i] = 0.0;
        init_ice[i] = 0.0;
        init_liq[i] = 0.0;
        supercool[i] = 0.0;
        EnergyRes[i] = 0.0;
        init_moist[i] = 0.0;
    }

    // 雪层水质量
    if (Nsnow > 0) {
        for (i = 0; i < Nsnow; i++) {
            mass_ice[i] = pack_ice[i];
            mass_liq[i] = pack_liq[i];
        }
    }
    // 土壤层水质量
    for (i = 0; i < Nnode; i++) {
        mass_ice[Nsnow + i] = (moist[i] - liq[i]) * dz_node[i] * MM_PER_M;
        mass_liq[Nsnow + i] = liq[i] * dz_node[i] * MM_PER_M; 
    }
    // 初始化状态
    for (i = 0; i < Nnode + Nsnow; i++) {
        init_ice[i] = mass_ice[i];
        init_liq[i] = mass_liq[i];
        init_moist[i] = mass_ice[i] + mass_liq[i];
    }

    // 计算土壤过冷水含量
    for (i = Nsnow; i < Nnode + Nsnow; i++) {
        j = i - Nsnow;
        if (T[i] < CONST_TKFRZ) {
            SuperCooling(j, T[i], liq[j], moist[j],
                         &supercool[i], soil_con);

            supercool[i] *= dz_node[j];
        }
    }

    // 确定融化或冻结状态
    for (i = 0; i < Nnode + Nsnow; i++) {
        if (mass_ice[i] > 0.0 && T[i] >= CONST_TKFRZ) {
            PhaseChange[i] = 1;
        }
        else if (mass_liq[i] > supercool[i] && T[i] < CONST_TKFRZ) {
            PhaseChange[i] = 2;
        }
        else if (Nsnow == 0 && snow->swq > 0.0 && i == Nsnow) {
            if (T[i] >= CONST_TKFRZ) {
                PhaseChange[i] = 1;
            }
        }
    }

    // 计算融化和冻结的能量盈余和损失
    for (i = 0; i < Nnode + Nsnow; i++) {
        if (PhaseChange[i] > 0) {
            EnergyRes[i] = (T[i] - CONST_TKFRZ) / fusion_fact[i];
            T[i] = CONST_TKFRZ;
        }
        if (PhaseChange[i] == 1 && EnergyRes[i] < 0.) {
            EnergyRes[i] = 0.0;
            PhaseChange[i] = 0;
        }
        if (PhaseChange[i] == 2 && EnergyRes[i] > 0.) {
            EnergyRes[i] = 0.0;
            PhaseChange[i] = 0;
        }
        fusion_flux[i] = EnergyRes[i] * step_dt / CONST_LATICE;
    }

    // 无分层雪的融化速率
    if (Nsnow == 0 && snow->swq > 0. && fusion_flux[0] > 0.) {
        old_swq = snow->swq;
        snow->swq = max(0.0, old_swq - fusion_flux[0]);
        ratio = snow->swq / old_swq;
        // 限制雪深在合理范围内
        snow->snow_depth = max(0.0, ratio * snow->snow_depth);
        snow->snow_depth = min(max(snow->snow_depth, snow->swq / 500.0), snow->swq / 50.0);
        RestTerm = EnergyRes[0] - CONST_LATICE * (old_swq - snow->swq) / step_dt;
        if (RestTerm > 0.0) {
            fusion_flux[0] = RestTerm * step_dt / CONST_LATICE;
            EnergyRes[0] = RestTerm;
        }
        else {
            fusion_flux[0] = 0.0;
            EnergyRes[0] = 0.0;
        }
        melt_grnd = max(0.0, old_swq - snow->swq) / step_dt;
        snow->pack_melt = old_swq - snow->swq;
        latent_fusion = CONST_LATICE * melt_grnd;
    }

    // 多层雪和土壤的融化和冻结速率
    for (i = 0; i < Nnode + Nsnow; i++) {
        bool is_snow = (i < Nsnow);
        if (PhaseChange[i] > 0 && fabs(EnergyRes[i]) > 0.0) {
            RestTerm = 0.0;
            if (fusion_flux[i] > 0.0) {
                mass_ice[i] = max(0.0, init_ice[i] - fusion_flux[i]);
                RestTerm = EnergyRes[i] - CONST_LATICE * (init_ice[i] - mass_ice[i]) / step_dt;
            }
            else if (fusion_flux[i] < 0.0) {
                if (is_snow) {      // snow layer
                    mass_ice[i] = min(init_moist[i], init_ice[i] - fusion_flux[i]);
                }
                else {              // soil layer
                    if (init_moist[i] < supercool[i]) {
                        mass_ice[i] = 0.0;
                    }
                    else {
                        mass_ice[i] = min(init_moist[i] - supercool[i], 
                                                init_ice[i] - fusion_flux[i]);
                        mass_ice[i] = max(mass_ice[i], 0.0);
                    }
                }
                RestTerm = EnergyRes[i] - CONST_LATICE * (init_ice[i] - mass_ice[i]) / step_dt;
            }
            mass_liq[i] = max(0.0, init_moist[i] - mass_ice[i]);

            if (fabs(RestTerm) > 0.) {
                T[i] += fusion_fact[i] * RestTerm;
                if (is_snow) {
                    if (mass_liq[i] * mass_ice[i] > 0.) {
                        T[i] = CONST_TKFRZ; // 相变界面温度为0
                    }
                    if (mass_ice[i] == 0.0) {
                        T[i] = CONST_TKFRZ;
                        EnergyRes[i + 1] += RestTerm;
                        fusion_flux[i + 1] *= step_dt / CONST_LATICE;
                    }
                }
            }
            latent_fusion += CONST_LATICE * (init_ice[i] - mass_ice[i]) / step_dt;
            // 雪融化速率
            if (is_snow) {
                melt_grnd += max(0.0, (init_ice[i] - mass_ice[i])) / step_dt;
            }
        }
    }

    // 更新雪和土壤的冰和液态水含量
    if (Nsnow > 0) {
        for (i = 0; i < Nsnow; i++) {
            pack_ice[i] = mass_ice[i];
            pack_liq[i] = mass_liq[i];
        }
    }
    for (i = 0; i < Nnode; i++) {
        liq[i] = mass_liq[Nsnow + i] / (dz_node[i] * MM_PER_M);
        moist[i] = (mass_liq[Nsnow + i] + mass_ice[Nsnow + i]) / (dz_node[i] * MM_PER_M);
    }
	// 将组合温度T写回各自的温度数组中
	if (Nsnow > 0) {
	    // 写回雪层温度（需要反转）
	    for (i = 0; i < Nsnow; i++) {
	        pack_T[i] = T[i];
	    }
	}

	// 写回土层温度（直接映射）
	for (i = 0; i < Nnode; i++) {
	    soil_T[i] = T[Nsnow + i];
	}
    return (0);
}
