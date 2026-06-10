/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate snow pack energy balance
 *****************************************************************************/

#include "vic_run.h"

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
	/* Initialize variables */
	size_t i, idx;
	size_t Nsnow = snow->Nsnow;
	size_t Nnode = cell->Nnode;
	double mat_A[Nnode];
	double mat_B[Nnode];
	double mat_C[Nnode];
	double mat_RHS[Nnode];
	double *fact = energy->fact;
	double *T = energy->T;
	double *pack_T = snow->pack_T;
	double *soil_T = cell->soil_T;
	double *Cs_node = energy->Cs_node;
	double *dz_node = cell->dz_node;
	double *zc_node = cell->zc_node;
	double *Zsum_node = cell->Zsum_node;
	double *kappa_int = energy->kappa_int;
	double coverage = snow->coverage;
	double deriv_grnd = energy->deriv_grnd;
	
	/* initialization */
	for (i = 0; i < Nnode; i++) {
		mat_A[i] = 0.;
		mat_B[i] = 0.;
		mat_C[i] = 0.;
		mat_RHS[i] = 0.;
	}
	// 填充温度数组T
    for (i = 0; i < Nsnow; i++) {
        T[i] = pack_T[i];
    }
    for (i = Nsnow; i < Nnode; i++) {
		idx = i - Nsnow;
        T[i] = soil_T[idx];
    }

    double capr = 0.34;  // heat capacity ratio                   
    for (i = 0; i < Nnode; i++) {
		if (i == 0) {
			fact[i] = step_dt / Cs_node[i] /
							(0.5 * (zc_node[i] - Zsum_node[i] + 
									capr * (zc_node[i+1] - Zsum_node[i])));
		}
		else if (i <= Nnode - 1) {
			fact[i] = step_dt / (Cs_node[i] * dz_node[i]);
		}
    }
    double grnd_flux = energy->grnd_flux;
    for (i = 0; i < Nnode; i++) {
        if (Nsnow == 0 && i == 0) {
            mat_A[i] = 0.0;
            mat_B[i] = 1.0 + fact[i] * kappa_int[i] - fact[i] * deriv_grnd;
            mat_C[i] = -fact[i] * kappa_int[i];
            mat_RHS[i] = T[i] + fact[i] * (grnd_flux - deriv_grnd * T[i]);
        }
        else if (i == Nsnow) {
            mat_A[i] = -coverage * fact[i] * kappa_int[i-1];
            mat_B[i] = 1.0 + fact[i] * (kappa_int[i] +
                   coverage * kappa_int[i-1]) - (1.0 - coverage) * fact[i] * deriv_grnd;
            mat_C[i] = -fact[i] * kappa_int[i];
            mat_RHS[i] = T[i] + fact[i] * ((1 - coverage) * (grnd_flux - deriv_grnd * T[i]));
        }
        else if (i <= Nnode - 2) {
            mat_A[i] = -fact[i] * kappa_int[i-1];
            mat_B[i] = 1.0 + fact[i] * (kappa_int[i-1] + kappa_int[i]);
            mat_C[i] = -fact[i] * kappa_int[i];
            mat_RHS[i] = T[i];
        }
        else if (i == Nnode - 1) {
            mat_A[i] = -fact[i] * kappa_int[i-1];
            mat_B[i] = 1.0 + fact[i] * kappa_int[i-1];
            mat_C[i] = 0.0;
            mat_RHS[i] = T[i];
        }
    }

	// 调用LAPACK函数求解线性方程组
	int info = LAPACKE_dgtsv(LAPACK_COL_MAJOR, Nnode, 
							 1, mat_A+1, mat_B, mat_C, mat_RHS, Nnode);

	if (info != 0) {
		return ERROR;
	}

	/* 检查收敛 */
	double max_diff = 0.0;
	for (i = 0; i < Nnode; i++) {
		double diff = fabs(mat_RHS[i] - T[i]);
		if (diff > max_diff) {
			max_diff = diff;
		}
	}
	energy->delt_T = max_diff;

    for (i = 0; i < Nnode; i++) {
        T[i] = mat_RHS[i];
    }
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
