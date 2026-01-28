/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate snow and soil layer temperature.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief  
 * Solves the glacier energy balance.
 *****************************************************************************/
int
SoilHeatFlux(double   		    step_dt,
			 double             fcanopy,
			 double             longwave,
			 cell_data_struct  *cell,
			 energy_bal_struct *energy,
			 snow_data_struct  *snow,
			 soil_con_struct   *soil_con)

{	
	extern option_struct     options;
	extern parameters_struct param;
	size_t i, j;
	size_t Nnode = options.Nnode;
	size_t Nsnow = snow->Nsnow;
	size_t Total_Layer = Nnode + Nsnow;
	size_t Total_Nodes = MAX_NODES + MAX_SNOWS;
	double *T = energy->T;
	double *pack_T = snow->pack_T;
	double *soil_T = cell->soil_T;
	double *dz_soil = soil_con->dz_soil;
	double *dz_node = cell->dz_node;
	double *Zsum_node = cell->Zsum_node;
	double *zc_node = cell->zc_node;
	double *kappa_int = energy->kappa_int;
	double *kappa_node = energy->kappa_node;
	double *Cs_node = energy->Cs_node;
	double *dz_snow = snow->dz_snow;
	double fact[MAX_NODES + MAX_SNOWS];
	double depthInv[MAX_NODES + MAX_SNOWS];
	double grad_temp[MAX_NODES + MAX_SNOWS];
	double heat_capacity[MAX_NODES + MAX_SNOWS];
	double EnergyExcess[MAX_NODES + MAX_SNOWS];
	/* initialization */
	for (i = 0; i < Total_Nodes; i++) {
		depthInv[i] = 0.;
		grad_temp[i] = 0.;
		heat_capacity[i] = 0.;
		EnergyExcess[i] = 0.;
	}
    for (i = 0; i < Nsnow; i++) {
        T[i] = pack_T[i];
    }
    for (i = 0; i < Nnode; i++) {
        T[Nsnow + i] = soil_T[i];
    }

	double emit_grnd = param.EMISS_GRND * CONST_STEBOL * pow(energy->Tgrnd, 4);
	double emit_slope = 4.0 * param.EMISS_GRND * CONST_STEBOL * 
						pow(energy->Tgrnd, 3);
	double emit_snow = param.EMISS_SNOW * CONST_STEBOL * pow(T[0], 4);
	double emit_soil = param.EMISS_GRND * CONST_STEBOL * pow(T[Nsnow], 4);
	double emit_water = param.EMISS_H2O * CONST_STEBOL * pow(T[Nsnow], 4);

	energy->grnd_flux = energy->NetShortGrnd + (1 - fcanopy) * param.EMISS_GRND * longwave - emit_grnd - (energy->LatentGrnd);
	energy->snow_flux = energy->NetShortSnow + (1 - fcanopy) * param.EMISS_SNOW * longwave - emit_snow - (energy->LatentGrnd);
	energy->soil_flux = energy->NetShortSoil + (1 - fcanopy) * param.EMISS_GRND * longwave - emit_soil - (energy->LatentGrnd);

	double capr = 0.5;  // heat capacity ratio
	double fn[MAX_NODES + MAX_SNOWS];
    for (i = 0; i < Total_Layer; i++) {
		if (i == 0) {
			fact[i] = step_dt / heat_capacity[i] * dz_node[i] / 
							(0.5 * zc_node[i] - Zsum_node[i - 1] + 
									capr * (dz_node[i + 1] - Zsum_node[i - 1]));
			fn[i] = kappa_int[i] * (T[i + 1] - T[i]) / (zc_node[i + 1] - zc_node[i]);
		}
		else if (i < Total_Layer - 1) {
			fact[i] = step_dt / heat_capacity[i];
			fn[i] = kappa_int[i] * (T[i + 1] - T[i]) / (zc_node[i + 1] - zc_node[i]);
		}
		else if (i == Total_Layer - 1) {
			fact[i] = step_dt / heat_capacity[i];
			fn[i] = 0.;
		}

    }



	for (i = 0; i < Total_Layer; i++) {
		if (Nsnow == 0) {
			// 没有雪层，全部为土壤层
	        if (i == 0) {
	            // 第一个土壤层（表层）
	            heat_capacity[i] = dz_node[i] * Cs_node[i];
	            if (Total_Layer > 1) {
	                depthInv[i] = 2.0 / (dz_node[i] + dz_node[i + 1]);
	                grad_temp[i] = 2.0 * (T[i] - T[i + 1]) / 
	                               (dz_node[i] + dz_node[i + 1]);
	            }
	            EnergyExcess[i] = kappa_node[i] * grad_temp[i] - energy->grnd_flux;
	        }
	        else if (i < Total_Layer - 1) {
	            // 中间土壤层
	            heat_capacity[i] = dz_node[i] * Cs_node[i];
	            depthInv[i] = 2.0 / (dz_node[i] + dz_node[i + 1]);
	            grad_temp[i] = 2.0 * (T[i] - T[i + 1]) / 
	                           (dz_node[i] + dz_node[i + 1]);
	            EnergyExcess[i] = kappa_node[i] * grad_temp[i] - 
	                              kappa_node[i - 1] * grad_temp[i - 1];
	        }
	        else {
	            // 最后一个土壤层
	            heat_capacity[i] = dz_node[i] * Cs_node[i];
	            EnergyExcess[i] = -kappa_node[i - 1] * grad_temp[i - 1];
	        }
	    }
	    else {
	    	// 包含雪层和土壤层
			if (i == 0) {	// 顶层雪
				heat_capacity[i] = dz_snow[i] * Cs_node[i];

				if (Nsnow > 1) {  // 下面是雪层
					depthInv[i] = 2.0 / (dz_snow[i] + dz_snow[i - 1]);
					grad_temp[i] = 2.0 * (T[i] - T[i + 1]) /
												(dz_snow[i] + dz_snow[i - 1]);
				}
				else {  // // Nsnow = 1，下面是土层
					depthInv[i] = 2.0 / (dz_snow[i] + dz_node[0]);;
					grad_temp[i] = 2.0 * (T[i] - T[i + 1]) / 
															(dz_snow[i] + dz_node[0]);
				}
				EnergyExcess[i] = kappa_node[i] * grad_temp[i] - energy->grnd_flux;
			}
			else if (i < Nsnow - 1) {	// 中间层雪
				heat_capacity[i] = dz_snow[i] * Cs_node[i];
				depthInv[i] = 2.0 / (dz_snow[i] + dz_snow[i - 1]);
				grad_temp[i] = 2.0 * (T[i] - T[i + 1]) /
														(dz_snow[i] + dz_snow[i - 1]);
				EnergyExcess[i] = kappa_node[i] * grad_temp[i] - 
										  kappa_node[i - 1] * grad_temp[i - 1];
			}
			else if (i == Nsnow - 1) {	// 最底层雪,对应depth[0]
				heat_capacity[i] = dz_snow[i] * Cs_node[i];
				depthInv[i] = 2.0 / (dz_snow[i] + dz_node[0]);
				grad_temp[i] = 2.0 * (T[i] - T[i + 1]) / 
														(dz_snow[i] + dz_node[0]);
				EnergyExcess[i] = kappa_node[i] * grad_temp[i] - 
										  kappa_node[i - 1] * grad_temp[i - 1];				
			}
			else if (i >= Nsnow && i < Total_Layer - 1) {  // 上层土壤
				j = i - Nsnow;
				heat_capacity[i] = dz_node[j] * Cs_node[i];
				depthInv[i] = 2.0 / (dz_node[j] + dz_node[j + 1]);;
				grad_temp[i] = 2.0 * (T[i] - T[i + 1]) / 
														(dz_node[j] + dz_node[j + 1]);
				EnergyExcess[i] = kappa_node[i] * grad_temp[i] - 
										  kappa_node[i - 1] * grad_temp[i - 1];
			}
			else {	// 最底层土壤
				heat_capacity[i] = dz_node[i - Nsnow] * Cs_node[i];
				EnergyExcess[i] = -kappa_node[i - 1] * grad_temp[i - 1];
			}
		}
	}

}
