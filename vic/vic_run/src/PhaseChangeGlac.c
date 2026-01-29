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
PhaseChangeGlac(double             step_dt,
                energy_bal_struct *energy,
                cell_data_struct  *cell,
                snow_data_struct  *snow,
                soil_con_struct   *soil_con)
{
    extern option_struct     options;

    size_t  i, j, nidx;
    double  pond_melt;
    double  old_swq;
    double *T = energy->T;
    double *liq = cell->liq;
    double *moist = cell->moist;
    double *pack_T = snow->pack_T;
    double *soil_T = cell->soil_T;
    double *pack_ice = snow->pack_ice;
    double *pack_liq = snow->pack_liq;
    double *dz_soil = soil_con->dz_soil;
    int *PhaseChange = cell->PhaseChange;
    double *fusion_fact = energy->fusion_fact;
    double *fusion_flux = energy->fusion_flux;

    double  init_ice[MAX_NODES + MAX_SNOWS];
    double  init_liq[MAX_NODES + MAX_SNOWS];
    double  mass_ice[MAX_NODES + MAX_SNOWS];
    double  mass_liq[MAX_NODES + MAX_SNOWS];
    double  RestTerm[MAX_NODES + MAX_SNOWS];
    double  EnergyRes[MAX_NODES + MAX_SNOWS];
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
    for (i = 0; i < Nnode; i++) {
        mass_ice[Nsnow + i] = (moist[i] - liq[i]) * dz_soil[i] * MM_PER_M;
        mass_liq[Nsnow + i] = liq[i] * dz_soil[i] * MM_PER_M;
    }

    for (i = 0; i < Nnode + Nsnow; i++) {
        init_ice[i] = mass_ice[i];
        init_liq[i] = mass_liq[i];
        init_moist[i] = mass_ice[i] + mass_liq[i];
    }

    // determine melting or freezing state
    for (i = 0; i < Nsnow; i++) {
        if (mass_ice[i] > 0. && T[i] >= CONST_TKFRZ) {
            PhaseChange[i] = 1; // melting
        }
        if (mass_liq[i] > 0. && T[i] < CONST_TKFRZ) {
            PhaseChange[i] = 2; // freezing
        }
    }
    // Calculate the energy surplus and loss for melting and freezing
    for (i = 0; i < Nsnow; i++) {
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

    // The rate of melting and freezing for multi-layer snow
    for (i = 0; i < Nsnow; i++) {
        if (PhaseChange[i] > 0 && fabs(EnergyRes[i]) > 0.0) {
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
                T[i] += fusion_fact[i] * RestTerm[i];
                
                if (mass_liq[i] * mass_ice[i] > 0.) {
                    T[i] = CONST_TKFRZ;
                }  
            }
            
            latent_fusion += CONST_LATICE * (init_ice[i] - mass_ice[i]) / step_dt;
            // snow melting rate
            melt_grnd += max(0.0, (init_ice[i] - mass_ice[i])) / step_dt;
        }
    }

    // 确定融化或冻结状态
    for (i = Nsnow; i < Nnode + Nsnow; i++) {
        if (mass_ice[i] > 0.0 && T[i] >= CONST_TKFRZ) {
            PhaseChange[i] = 1;
        }
        else if (mass_liq[i] > 0. && T[i] < CONST_TKFRZ) {
            PhaseChange[i] = 2;
        }
        else if (Nsnow == 0 && snow->swq > 0.0 && i == Nsnow) {
            if (T[i] >= CONST_TKFRZ) {
                PhaseChange[i] = 1;
            }
        }
    }

    // 计算融化和冻结的能量盈余和损失
    for (i = Nsnow; i < Nnode + Nsnow; i++) {
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
    // ice layer water mass
    if (Nsnow == 0 && snow->swq > 0. && fusion_flux[1] > 0.) {
        old_swq = snow->swq;
        snow->swq = max(0.0, old_swq - fusion_flux[1]);
        double ratio = snow->swq / old_swq;
        snow->snow_depth = max(0.0, ratio * snow->snow_depth);
        snow->snow_depth = min(max(snow->snow_depth, snow->swq / 500.0), snow->swq / 50.0);
        RestTerm[0] = EnergyRes[0] - CONST_LATICE * (old_swq - snow->swq) / step_dt;
        if (RestTerm[0] > 0.0) {
            fusion_flux[0] = RestTerm[0] * step_dt / CONST_LATICE;
            EnergyRes[0] = RestTerm[0];
            PhaseChange[0] = 1;
        }
        else {
            fusion_flux[1] = 0.0;
            EnergyRes[1] = 0.0;
            PhaseChange[0] = 0;
        }
        snow->pack_melt = old_swq - snow->swq;
    }
    // The rate of melting and freezing for glacier ice
    for (i = Nsnow; i < Nnode + Nsnow; i++) {
        if (PhaseChange[i] > 0 && fabs(EnergyRes[i]) > 0.0) {
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
                T[i] += fusion_fact[i] * RestTerm[i];
            }
            latent_fusion += CONST_LATICE * (init_ice[i] - mass_ice[i]) / step_dt;
        }
    }

    // Deal with residuals in ice/soil
    if (any(T, Nnode, > CONST_TKFRZ) && any(T, Nnode, < CONST_TKFRZ)) {
        for (i = Nsnow; i < Nnode + Nsnow; i++) {
            if (T[i] > CONST_TKFRZ) {
                RestTerm[i] = (T[i] - CONST_TKFRZ) / fusion_fact[i];
                for (j = Nsnow; j < Nnode + Nsnow; j++) {
                    if (i != j && T[j] < CONST_TKFRZ && RestTerm[i] > 0.1) {
                        RestTerm[j] = (T[j] - CONST_TKFRZ) / fusion_fact[j];
                        if (fabs(RestTerm[j]) > RestTerm[i]) {
                            RestTerm[j] += RestTerm[i];
                            T[j] = CONST_TKFRZ + RestTerm[j] * fusion_flux[j];
                            RestTerm[i] = 0.0;
                        }
                        else {
                            RestTerm[i] += RestTerm[j];
                            RestTerm[j] = 0.0;
                            T[j] = CONST_TKFRZ;
                        }
                    }
                }
                T[i] = CONST_TKFRZ + RestTerm[i] * fusion_fact[i];
            }
        }
    }

    // now remove excess cold by increasing temperture
    if (any(T, Nnode, > CONST_TKFRZ) && any(T, Nnode, < CONST_TKFRZ)) {
        for (i = Nsnow; i < Nnode + Nsnow; i++) {
            if (T[i] < CONST_TKFRZ) {
                RestTerm[i] = (T[i] - CONST_TKFRZ) / fusion_fact[i];
                for (j = Nsnow; j < Nnode + Nsnow; j++) {
                    if (i != j && T[j] > CONST_TKFRZ && RestTerm[i] < 0.1) {
                        RestTerm[j] = (T[j] - CONST_TKFRZ) / fusion_fact[j];
                        if (RestTerm[j] > fabs(RestTerm[i])) {
                            RestTerm[j] += RestTerm[i];
                            T[j] = CONST_TKFRZ + RestTerm[j] * fusion_fact[j];
                            RestTerm[i] = 0.0;
                        }
                        else {
                            RestTerm[i] += RestTerm[j];
                            RestTerm[j] = 0.0;
                            T[j] = CONST_TKFRZ;
                        }
                    }
                }
                T[i] = CONST_TKFRZ + RestTerm[i] * fusion_fact[i];
            }
        }
    }

    // now remove excess heat by melting ice
    if (any(T, Nnode, > CONST_TKFRZ) && any(mass_ice, Nnode, > 0.)) {
        for (i = Nsnow; i < Nnode + Nsnow; i++) {
            if (T[i] < CONST_TKFRZ) {
                RestTerm[i] = (T[i] - CONST_TKFRZ) / fusion_fact[i];
                fusion_flux[i] = RestTerm[i] * step_dt / CONST_LATICE;
                for (j = Nsnow; j < Nnode + Nsnow; j++) {
                    if (i != j && mass_ice[j] > 0. && fusion_flux[i] > 0.1) {
                        if (mass_ice[j] > fusion_flux[i]) {
                            mass_ice[j] -= fusion_flux[i];
                            latent_fusion += CONST_LATICE * fusion_flux[i] / step_dt;
                            T[j] = CONST_TKFRZ;
                            fusion_flux[i] = 0.0;
                        }
                        else {
                            fusion_flux[i] -= mass_ice[j];
                            latent_fusion += CONST_LATICE * mass_ice[j] / step_dt;
                            mass_ice[j] = 0.0;
                            T[j] = CONST_TKFRZ;
                        }
                        mass_liq[j] = max(0.0, init_moist[j] - mass_ice[j]);
                    }
                }
                RestTerm[i] = fusion_flux[i] * CONST_LATICE / step_dt;
                T[i] = CONST_TKFRZ + RestTerm[i] * fusion_fact[i];
            }
        }
    }

    // snow remove excess cold by refreezing liquid (may not be necessary with above loop)
    if (any(T, Nnode, < CONST_TKFRZ) && any(mass_liq, Nnode, > 0.)) {
        for (i = Nsnow; i < Nnode + Nsnow; i++) {
            if (T[i] < CONST_TKFRZ) {
                RestTerm[i] = (T[i] - CONST_TKFRZ) / fusion_fact[i];
                fusion_flux[i] = RestTerm[i] * step_dt / CONST_LATICE;
                for (j = Nsnow; j < Nnode + Nsnow; j++) {
                    if (i != j && mass_liq[j] > 0. && fusion_flux[i] < -0.1) {
                        if (mass_liq[j] > fabs(fusion_flux[i])) {
                            mass_ice[j] -= fusion_flux[i];
                            latent_fusion += CONST_LATICE * fusion_flux[i] / step_dt;
                            T[j] = CONST_TKFRZ;
                            fusion_flux[i] = 0.0;
                        }
                        else {
                            fusion_flux[i] += mass_ice[j];
                            latent_fusion -= CONST_LATICE * mass_liq[j] / step_dt;
                            mass_ice[j] = init_moist[j];
                            T[j] = CONST_TKFRZ;
                        }
                        mass_liq[j] = max(0.0, init_moist[j] - mass_ice[j]);
                    }
                }
                RestTerm[i] = fusion_flux[i] * CONST_LATICE / step_dt;
                T[i] = CONST_TKFRZ + RestTerm[i] * fusion_fact[i];
            }
        }
    }

    // update snow and soil ice and liquid content
    if (Nsnow > 0) {
        for (i = 0; i < Nsnow; i++) {
            pack_ice[i] = mass_ice[i];
            pack_liq[i] = mass_liq[i];
        }
    }
    for (i = 0; i < Nnode; i++) {    // glacier
        liq[i] = mass_liq[Nsnow + i] / (dz_soil[i] * MM_PER_M);
        liq[i] = max(0.0, min(1.0, liq[i]));
        moist[i] = 1.0;
    }
	// 将组合温度T写回各自的温度数组中
	if (Nsnow > 0) {
	    // 写回雪层温度
	    for (i = 0; i < Nsnow; i++) {
	        pack_T[i] = T[i];
	    }
	}

	// 写回土层温度
	for (i = 0; i < Nnode; i++) {
	    soil_T[i] = T[Nsnow + i];
	}
    return (0);
}
        


    
