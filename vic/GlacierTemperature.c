/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate snow pack energy balance
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief  
 * Solves the glacier energy balance. See SnowPackEnergyBalance for comparison.
 *****************************************************************************/
int
GlacierTemperature(double   		  step_dt,
				   cell_data_struct  *cell,
				   energy_bal_struct *energy,
				   snow_data_struct  *snow,
				   soil_con_struct   *soil_con)

{	
	extern option_struct       options;

	size_t i, j, lidx;
	double mat_A[MAX_NODES + MAX_SNOWS];
	double mat_B[MAX_NODES + MAX_SNOWS];
	double mat_C[MAX_NODES + MAX_SNOWS];
	double mat_RHS[MAX_NODES + MAX_SNOWS];
	double depthInv[MAX_NODES + MAX_SNOWS];
	double grad_temp[MAX_NODES + MAX_SNOWS];
	double heat_capacity[MAX_NODES + MAX_SNOWS];
	double EnergyExcess[MAX_NODES + MAX_SNOWS];
	double grnd_flux = energy->grnd_flux;
	double *depth = snow->dz_snow;
	double *dz_node = soil_con->dz_soil;
	double *Cs_node = energy->Cs_node;
	double *T = energy->T;
	double *pack_T = snow->pack_T;
	double *soil_T = cell->soil_T;
	double *kappa_node = energy->kappa_node;
	size_t Total_Nodes = MAX_NODES + MAX_SNOWS;
	size_t Nlayer = options.Nnode + snow->Nsnow;
	size_t Nsnow = snow->Nsnow;
	size_t Nnode = options.Nnode;
	/* initialization */
	for (i = 0; i < Total_Nodes; i++) {
		mat_A[i] = 0.;
		mat_B[i] = 0.;
		mat_C[i] = 0.;
		mat_RHS[i] = 0.;
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

	for (i = 0; i < Nlayer; i++) {
		if (Nsnow == 0) {
			// 没有雪层，全部为土壤层
	        if (i == 0) {
	            // 第一个土壤层（表层）
	            heat_capacity[i] = dz_node[i] * Cs_node[i];
	            if (Nlayer > 1) {
	                depthInv[i] = 2.0 / (dz_node[i] + dz_node[i + 1]);
	                grad_temp[i] = 2.0 * (T[i] - T[i + 1]) / 
	                               (dz_node[i] + dz_node[i + 1]);
	            }
	            EnergyExcess[i] = kappa_node[i] * grad_temp[i] - grnd_flux;
	        }
	        else if (i < Nlayer - 1) {
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
				lidx = Nsnow - 1;
				heat_capacity[i] = depth[lidx] * Cs_node[i];

				if (Nsnow > 1) {  // 下面是雪层
					depthInv[i] = 2.0 / (depth[lidx] + depth[lidx - 1]);
					grad_temp[i] = 2.0 * (T[i] - T[i + 1]) /
												(depth[lidx] + depth[lidx - 1]);
				}
				else {  // // Nsnow = 1，下面是土层
					depthInv[i] = 2.0 / (depth[lidx] + dz_node[0]);;
					grad_temp[i] = 2.0 * (T[i] - T[i + 1]) / 
															(depth[lidx] + dz_node[0]);
				}
				EnergyExcess[i] = kappa_node[i] * grad_temp[i] - grnd_flux;
			}
			else if (i < Nsnow - 1) {	// 中间层雪
				lidx = Nsnow - 1 - i;
				heat_capacity[i] = depth[lidx] * Cs_node[i];
				depthInv[i] = 2.0 / (depth[lidx] + depth[lidx - 1]);;
				grad_temp[i] = 2.0 * (T[i] - T[i + 1]) / 
														(depth[lidx] + depth[lidx - 1]);
				EnergyExcess[i] = kappa_node[i] * grad_temp[i] - 
										  kappa_node[i - 1] * grad_temp[i - 1];
			}
			else if (i == Nsnow - 1) {	// 最底层雪,对应depth[0]
				heat_capacity[i] = depth[0] * Cs_node[i];
				depthInv[i] = 2.0 / (depth[0] + dz_node[0]);
				grad_temp[i] = 2.0 * (T[i] - T[i + 1]) / 
														(depth[0] + dz_node[0]);
				EnergyExcess[i] = kappa_node[i] * grad_temp[i] - 
										  kappa_node[i - 1] * grad_temp[i - 1];				
			}
			else if (i >= Nsnow && i < Nlayer - 1) {  // 上层土壤
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

	/* prepare the matrix coefficients for the tri-diagonal matrix */
	for (i = 0; i < Nlayer; i++) {
		if (i == 0) {
			mat_A[i] = 0.;
			mat_C[i] = -kappa_node[i] * depthInv[i] / heat_capacity[i];
			mat_B[i] = -mat_C[i];
		}
		else if (i < Nlayer - 1) {
			mat_A[i] = -kappa_node[i - 1] * depthInv[i - 1] / heat_capacity[i];
			mat_C[i] = -kappa_node[i] * depthInv[i] / heat_capacity[i];
			mat_B[i] = -(mat_A[i] + mat_C[i]);
		}
		else {
			mat_A[i] = -kappa_node[i - 1] * depthInv[i - 1] / heat_capacity[i];
			mat_C[i] = 0.;
			mat_B[i] = -(mat_A[i] + mat_C[i]);			
		}
		mat_RHS[i] = -EnergyExcess[i] / heat_capacity[i];
	}
	/* compute soil temperatures */
	solve_temperature(mat_A, mat_B, mat_C,
					  mat_RHS, step_dt,
					  Nlayer, cell, energy, snow);
	return(0);
}
