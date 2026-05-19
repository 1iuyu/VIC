/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate the phase change (melting/freezing) of snow water and soil water.
 *****************************************************************************/

#include "vic_run.h"

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

    size_t  i, j, lidx;
    double  RestTerm;
    double  old_swq;
    double  ratio;
    double *liq = cell->liq;
    double *ice = cell->ice;
    double *moist = cell->moist;
    double *pack_T = snow->pack_T;
    double *soil_T = cell->soil_T;
    double *fact = energy->fact;
    double *pack_ice = snow->pack_ice;
    double *pack_liq = snow->pack_liq;
    double *dz_node = cell->dz_node;          // layer and thermal node thickness (m)
    int *MELTING = energy->MELTING;         // 0 = none, 1 = melt, 2 = freeze
    double *fusion_flux = energy->fusion_flux;
    double *bexp_node = soil_con->bexp_node;
    double *Wsat_node = soil_con->Wsat_node;
    double *psisat_node = soil_con->psisat_node;
    double init_ice[MAX_THERM];    // initial ice content of each thermal node (kg/m^2)
    double init_liq[MAX_THERM];
    double mass_ice[MAX_THERM];
    double mass_liq[MAX_THERM];
    double EnergyRes[MAX_THERM];
    double supercool[MAX_THERM];
    double init_moist[MAX_THERM];

    // 数组大小计算
    size_t Nsoil = cell->Nsoil;
    size_t Nsnow = snow->Nsnow;
    size_t Ntherm = Nsoil + Nsnow;
    double latent_fusion = 0.0;
    double pack_melt = 0.0; // 雪层融化的水量[mm]
    double pack_frze = 0.0; // 雪层冻结的水量[mm]
    
    // 初始化数组
    for (i = 0; i < MAX_THERM; i++) {
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
    for (i = 0; i < Nsoil; i++) {
        lidx = Nsnow + i;
        mass_ice[lidx] = (moist[i] - liq[i]) * dz_node[lidx] * MM_PER_M;
        mass_liq[lidx] = liq[i] * dz_node[lidx] * MM_PER_M; 
    }
    // 初始化状态
    for (i = 0; i < Ntherm; i++) {
        init_ice[i] = mass_ice[i];
        init_liq[i] = mass_liq[i];
        init_moist[i] = mass_ice[i] + mass_liq[i];
    }

    // 计算土壤过冷水含量
/*    for (i = 0; i < Nsoil; i++) {
        j = Nsnow + i;
        if (soil_T[i] < CONST_TKFRZ) {
            supercool[j] = SuperCooling(soil_T[i], 
                                        liq[i], moist[i],
                                        bexp_node[i],
                                        Wsat_node[i],
                                        psisat_node[i]);

            supercool[j] *= dz_node[j] * MM_PER_M;
        }
    } */

    // 确定融化或冻结状态
    for (i = 0; i < Nsnow; i++) {
        if (mass_ice[i] > 0.0 && pack_T[i] >= CONST_TKFRZ) {
            MELTING[i] = 1;
        }
        else if (mass_liq[i] > 0.0 && pack_T[i] < CONST_TKFRZ) {
            MELTING[i] = 2;
        }
    }
    for (i = Nsnow; i < Ntherm; i++) {
        j = i - Nsnow;
        if (mass_ice[i] > 0.0 && soil_T[j] >= CONST_TKFRZ) {
            MELTING[i] = 1;
        }
        else if (mass_liq[i] > supercool[i] && soil_T[j] < CONST_TKFRZ) {
            MELTING[i] = 2;
        }
        else if (Nsnow == 0 && snow->swq > 0.0 && i == Nsnow) {
            if (soil_T[i] >= CONST_TKFRZ) {
                MELTING[i] = 1;
            }
        }
    }
    // 计算融化和冻结的能量盈余和损失
    for (i = 0; i < Ntherm; i++) {
        lidx = i - Nsnow;
        if (MELTING[i] > 0) {
            if (i < Nsnow) {
                EnergyRes[i] = (pack_T[i] - CONST_TKFRZ) / fact[i];
                pack_T[i] = CONST_TKFRZ;
            }
            else {
                EnergyRes[i] = (soil_T[i] - CONST_TKFRZ) / fact[i];
                soil_T[lidx] = CONST_TKFRZ;
            }
        }
        // 特殊情况
        if (MELTING[i] == 1 && EnergyRes[i] < 0.) {
            EnergyRes[i] = 0.0;
            MELTING[i] = 0;
        }
        if (MELTING[i] == 2 && EnergyRes[i] > 0.) {
            EnergyRes[i] = 0.0;
            MELTING[i] = 0;
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
        pack_melt = max(0.0, old_swq - snow->swq);
        latent_fusion = CONST_LATICE * pack_melt / step_dt;
    }

    // 多层雪和土壤的融化和冻结速率
    for (i = 0; i < Ntherm; i++) {
        bool is_snow = (i < Nsnow);
        lidx = i - Nsnow;
        if (MELTING[i] > 0 && fabs(EnergyRes[i]) > 0.0) {
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
                if (is_snow) {
                    pack_T[i] += fact[i] * RestTerm;
                }
                else {
                    soil_T[lidx] += fact[i] * RestTerm;
                }
                if (is_snow) {
                    if (mass_liq[i] * mass_ice[i] > 0.) {
                        pack_T[i] = CONST_TKFRZ; // 相变界面温度为0
                    }
                    if (mass_ice[i] == 0.) {
                        pack_T[i] = CONST_TKFRZ; // 相变界面温度为0
                        EnergyRes[i+1] += RestTerm; 
                        fusion_flux[i+1] = EnergyRes[i+1] * step_dt / CONST_LATICE;
                    }
                }
            }
            latent_fusion += CONST_LATICE * (init_ice[i] - mass_ice[i]) / step_dt;
            // 雪融化速率
            if (MELTING[i] == 1 && is_snow) {
                pack_melt += max(0.0, (init_ice[i] - mass_ice[i]));
            }
            if (MELTING[i] == 2 && is_snow) {
                pack_frze += max(0.0, (mass_ice[i] - init_ice[i]));
            }
        }
    }
    // 写回结构体
    snow->pack_melt = pack_melt;
    snow->pack_frze = pack_frze;

    // 更新雪和土壤的冰和液态水含量
    if (Nsnow > 0) {
        for (i = 0; i < Nsnow; i++) {
            pack_ice[i] = mass_ice[i];
            pack_liq[i] = mass_liq[i];
        }
    }
    for (i = 0; i < Nsoil; i++) {
        lidx = Nsnow + i;
        liq[i] = mass_liq[lidx] / (dz_node[lidx] * MM_PER_M);
        ice[i] = mass_ice[lidx] / (dz_node[lidx] * MM_PER_M);
        moist[i] = (mass_liq[lidx] + mass_ice[lidx]) / (dz_node[lidx] * MM_PER_M);
    }

    return (0);
}