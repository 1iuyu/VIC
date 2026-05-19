/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate the phase change (melting/freezing) of glacier water.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    Calculate the phase change (melting/freezing) of glacier water.
 *****************************************************************************/
int
PhaseChangeGlac(double             step_dt,
                energy_bal_struct *energy,
                cell_data_struct  *cell,
                snow_data_struct  *snow,
                soil_con_struct   *soil_con)
{
    extern option_struct     options;

    size_t  i, j, lidx;
    double  old_swq;
    double *liq = cell->liq;
    double *ice = cell->ice;
    double *moist = cell->moist;
    double *pack_T = snow->pack_T;
    double *soil_T = cell->soil_T;
    double *pack_ice = snow->pack_ice;
    double *pack_liq = snow->pack_liq;
    double *dz_soil = soil_con->dz_soil;
    int *MELTING = energy->MELTING;
    double *fusion_flux = energy->fusion_flux;
    double *fact = energy->fact;
    double init_ice[MAX_THERM];
    double init_liq[MAX_THERM];
    double mass_ice[MAX_THERM];
    double mass_liq[MAX_THERM];
    double RestTerm[MAX_THERM];
    double EnergyRes[MAX_THERM];
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
        RestTerm[i] = 0.0;
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
        lidx = i + Nsnow;
        mass_ice[lidx] = (moist[i] - liq[i]) * dz_soil[i] * MM_PER_M;
        mass_liq[lidx] = liq[i] * dz_soil[i] * MM_PER_M;
    }

    for (i = 0; i < Ntherm; i++) {
        init_ice[i] = mass_ice[i];
        init_liq[i] = mass_liq[i];
        init_moist[i] = mass_ice[i] + mass_liq[i];
    }

    // determine melting or freezing state
    for (i = 0; i < Nsnow; i++) {
        if (mass_ice[i] > 0. && pack_T[i] >= CONST_TKFRZ) {
            MELTING[i] = 1; // melting
        }
        if (mass_liq[i] > 0. && pack_T[i] < CONST_TKFRZ) {
            MELTING[i] = 2; // freezing
        }
    }
    for (i = Nsnow; i < Ntherm; i++) {
        j = i - Nsnow;
        if (mass_ice[i] > 0.0 && soil_T[j] >= CONST_TKFRZ) {
            MELTING[i] = 1;
        }
        else if (mass_liq[i] > 0.0 && soil_T[j] < CONST_TKFRZ) {
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

    // The rate of melting and freezing for multi-layer snow
    for (i = 0; i < Nsnow; i++) {
        if (MELTING[i] > 0 && fabs(EnergyRes[i]) > 0.0) {
            RestTerm[i] = 0.0;

            if (fusion_flux[i] > 0.0) {
                mass_ice[i] = max(0.0, init_ice[i] - fusion_flux[i]);
                RestTerm[i] = EnergyRes[i] - CONST_LATICE * (init_ice[i] - mass_ice[i]) / step_dt;
            }
            else if (fusion_flux[i] < 0.0) {
                mass_ice[i] = min(init_moist[i], init_ice[i] - fusion_flux[i]);
                RestTerm[i] = EnergyRes[i] - CONST_LATICE * (init_ice[i] - mass_ice[i]) / step_dt;
            }
            mass_liq[i] = max(0.0, init_moist[i] - mass_ice[i]);
            if (fabs(RestTerm[i]) > 0.) {
                pack_T[i] += fact[i] * RestTerm[i];
                
                if (mass_liq[i] * mass_ice[i] > 0.) {
                    pack_T[i] = CONST_TKFRZ;
                }  
            }
            
            latent_fusion += CONST_LATICE * (init_ice[i] - mass_ice[i]) / step_dt;
            // snow melting rate
            // 雪融化速率
            if (MELTING[i] == 1) {
                pack_melt += max(0.0, (init_ice[i] - mass_ice[i]));
            }
            if (MELTING[i] == 2) {
                pack_frze += max(0.0, (mass_ice[i] - init_ice[i]));
            }
        }
    }
    // ice layer water mass
    if (Nsnow == 0 && snow->swq > 0. && fusion_flux[0] > 0.) {
        old_swq = snow->swq;
        snow->swq = max(0.0, old_swq - fusion_flux[0]);
        double ratio = snow->swq / old_swq;
        snow->snow_depth = max(0.0, ratio * snow->snow_depth);
        snow->snow_depth = min(max(snow->snow_depth, snow->swq / 500.0), snow->swq / 50.0);
        RestTerm[0] = EnergyRes[0] - CONST_LATICE * (old_swq - snow->swq) / step_dt;
        if (RestTerm[0] > 0.0) {
            fusion_flux[0] = RestTerm[0] * step_dt / CONST_LATICE;
            EnergyRes[0] = RestTerm[0];
            MELTING[0] = 1;
        }
        else {
            fusion_flux[0] = 0.0;
            EnergyRes[1] = 0.0;
            MELTING[0] = 0;
        }
        pack_melt = max(0.0, old_swq - snow->swq);
        latent_fusion = CONST_LATICE * pack_melt / step_dt;
    }
    // The rate of melting and freezing for glacier ice
    for (i = Nsnow; i < Ntherm; i++) {
        lidx = i - Nsnow;
        if (MELTING[i] > 0 && fabs(EnergyRes[i]) > 0.0) {
            RestTerm[i] = 0.0;

            if (fusion_flux[i] > 0.0) {
                mass_ice[i] = max(0.0, init_ice[i] - fusion_flux[i]);
                RestTerm[i] = EnergyRes[i] - CONST_LATICE * (init_ice[i] - mass_ice[i]) / step_dt;
            }
            else if (fusion_flux[i] < 0.0) {
                mass_ice[i] = min(init_moist[i], init_ice[i] - fusion_flux[i]);
                RestTerm[i] = EnergyRes[i] - CONST_LATICE * (init_ice[i] - mass_ice[i]) / step_dt;
            }

            mass_liq[i] = max(0.0, init_moist[i] - mass_ice[i]);
            if (fabs(RestTerm[i]) > 0.) {
                soil_T[lidx] += fact[i] * RestTerm[i];
            }
            latent_fusion += CONST_LATICE * (init_ice[i] - mass_ice[i]) / step_dt;
        }
    }

    // Deal with residuals in ice/soil
    if (any(soil_T, Nsoil, > CONST_TKFRZ) && any(soil_T, Nsoil, < CONST_TKFRZ)) {
        for (i = 0; i < Nsoil; i++) {
            if (soil_T[i] > CONST_TKFRZ) {
                RestTerm[i] = (soil_T[i] - CONST_TKFRZ) / fact[i];
                for (j = 0; j < Nsoil; j++) {
                    if (i != j && soil_T[j] < CONST_TKFRZ && RestTerm[i] > 0.1) {
                        RestTerm[j] = (soil_T[j] - CONST_TKFRZ) / fact[j];
                        if (fabs(RestTerm[j]) > RestTerm[i]) {
                            RestTerm[j] += RestTerm[i];
                            soil_T[j] = CONST_TKFRZ + RestTerm[j] * fusion_flux[j];
                            RestTerm[i] = 0.0;
                        }
                        else {
                            RestTerm[i] += RestTerm[j];
                            RestTerm[j] = 0.0;
                            soil_T[j] = CONST_TKFRZ;
                        }
                    }
                }
                soil_T[i] = CONST_TKFRZ + RestTerm[i] * fact[i];
            }
        }
    }

    // now remove excess cold by increasing temperture
    if (any(soil_T, Nsoil, > CONST_TKFRZ) && any(soil_T, Nsoil, < CONST_TKFRZ)) {
        for (i = 0; i < Nsoil; i++) {
            if (soil_T[i] < CONST_TKFRZ) {
                RestTerm[i] = (soil_T[i] - CONST_TKFRZ) / fact[i];
                for (j = 0; j < Nsoil; j++) {
                    if (i != j && soil_T[j] > CONST_TKFRZ && RestTerm[i] < 0.1) {
                        RestTerm[j] = (soil_T[j] - CONST_TKFRZ) / fact[j];
                        if (RestTerm[j] > fabs(RestTerm[i])) {
                            RestTerm[j] += RestTerm[i];
                            soil_T[j] = CONST_TKFRZ + RestTerm[j] * fact[j];
                            RestTerm[i] = 0.0;
                        }
                        else {
                            RestTerm[i] += RestTerm[j];
                            RestTerm[j] = 0.0;
                            soil_T[j] = CONST_TKFRZ;
                        }
                    }
                }
                soil_T[i] = CONST_TKFRZ + RestTerm[i] * fact[i];
            }
        }
    }

    // now remove excess heat by melting ice
    if (any(soil_T, Nsoil, > CONST_TKFRZ) && any(mass_ice, Nsoil, > 0.)) {
        for (i = 0; i < Nsoil; i++) {
            if (soil_T[i] < CONST_TKFRZ) {
                RestTerm[i] = (soil_T[i] - CONST_TKFRZ) / fact[i];
                fusion_flux[i] = RestTerm[i] * step_dt / CONST_LATICE;
                for (j = 0; j < Nsoil; j++) {
                    if (i != j && mass_ice[j] > 0. && fusion_flux[i] > 0.1) {
                        if (mass_ice[j] > fusion_flux[i]) {
                            mass_ice[j] -= fusion_flux[i];
                            latent_fusion += CONST_LATICE * fusion_flux[i] / step_dt;
                            soil_T[j] = CONST_TKFRZ;
                            fusion_flux[i] = 0.0;
                        }
                        else {
                            fusion_flux[i] -= mass_ice[j];
                            latent_fusion += CONST_LATICE * mass_ice[j] / step_dt;
                            mass_ice[j] = 0.0;
                            soil_T[j] = CONST_TKFRZ;
                        }
                        mass_liq[j] = max(0.0, init_moist[j] - mass_ice[j]);
                    }
                }
                RestTerm[i] = fusion_flux[i] * CONST_LATICE / step_dt;
                soil_T[i] = CONST_TKFRZ + RestTerm[i] * fact[i];
            }
        }
    }

    // snow remove excess cold by refreezing liquid (may not be necessary with above loop)
    if (any(soil_T, Nsoil, < CONST_TKFRZ) && any(mass_liq, Nsoil, > 0.)) {
        for (i = 0; i < Nsoil; i++) {
            lidx = i + Nsnow;
            if (soil_T[i] < CONST_TKFRZ) {
                // 计算该层需要的冻结能量（负值表示冷量亏缺）
                RestTerm[i] = (soil_T[i] - CONST_TKFRZ) / fact[lidx];
                fusion_flux[lidx] = RestTerm[i] * step_dt / CONST_LATICE;
                
                // 遍历所有层寻找可冻结的液态水
                for (j = 0; j < Nsoil; j++) {
                    size_t jdx = j + Nsnow;
                    if (i != j && mass_liq[jdx] > 0. && fusion_flux[lidx] < -1e-10) {
                        if (mass_liq[jdx] >= fabs(fusion_flux[lidx])) {
                            // j层有足够液态水满足i层的冻结需求
                            mass_liq[jdx] += fusion_flux[lidx];  // fusion_flux是负的
                            mass_ice[jdx] -= fusion_flux[lidx];
                            latent_fusion += CONST_LATICE * fusion_flux[lidx] / step_dt;
                            fusion_flux[lidx] = 0.0;
                        }
                        else {
                            // j层液态水不够，全部冻结
                            fusion_flux[lidx] += mass_liq[jdx];  // 减少剩余冷量需求
                            latent_fusion -= CONST_LATICE * mass_liq[jdx] / step_dt;
                            mass_ice[jdx] += mass_liq[jdx];
                            mass_liq[jdx] = 0.0;
                        }
                        // 更新j层的液态水含量
                        mass_liq[jdx] = max(0.0, init_moist[jdx] - mass_ice[jdx]);
                    }
                }
                
                // 根据剩余的冷量需求更新i层温度
                RestTerm[i] = fusion_flux[lidx] * CONST_LATICE / step_dt;
                soil_T[i] = CONST_TKFRZ + RestTerm[i] * fact[lidx];
                
                // 如果冷量被完全满足，温度应该达到冻结点
                if (fabs(fusion_flux[lidx]) < 1e-10) {
                    soil_T[i] = CONST_TKFRZ;
                }
            }
        }
    }
    // 写回结构体
    snow->pack_melt = pack_melt;
    snow->pack_frze = pack_frze;
    
    // update snow and soil ice and liquid content
    if (Nsnow > 0) {
        for (i = 0; i < Nsnow; i++) {
            pack_ice[i] = mass_ice[i];
            pack_liq[i] = mass_liq[i];
        }
    }
    for (i = 0; i < Nsoil; i++) {    // glacier
        lidx = Nsnow + i;
        liq[i] = mass_liq[lidx] / (dz_soil[i] * MM_PER_M);
        liq[i] = max(0.0, min(1.0, liq[i]));
        ice[i] = mass_ice[lidx] / (dz_soil[i] * MM_PER_M);
        moist[i] = 1.0;
    }

    return (0);
}
        


    
