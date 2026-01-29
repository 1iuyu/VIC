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
SoilTemperature(double   		   step_dt,
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
	double grnd_flux = energy->grnd_flux;
	double snow_flux = energy->snow_flux;
	double flux_slope = energy->flux_slope;
	double *fact = energy->fact;
	double *fn = energy->fn;
	double *T = energy->T;
	double *kappa_int = energy->kappa_int;
	double *zc_node = cell->zc_node;
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
	}
	double dzp, dzm;
	/* prepare the matrix coefficients for the tri-diagonal matrix */
	for (i = 0; i < Nsnow; i++) {
		if (i == 0) {
			dzp = zc_node[i + 1] - zc_node[i];
			mat_RHS[i + 1] = T[i] + fact[i] * (snow_flux - flux_slope * T[i] + fn[i]);
		}
		else {
			dzm = zc_node[i] - zc_node[i - 1];
			dzp = zc_node[i + 1] - zc_node[i];
		}

	}
	for (i = 0; i < Nnode; i++) {
		lidx = Nsnow + i;
		if (i == 0) {
			dzp = zc_node[lidx + 1] - zc_node[lidx];
			dzm = zc_node[lidx] - zc_node[lidx - 1];
			mat_RHS[lidx] = T[lidx] + fact[lidx] * (grnd_flux - flux_slope * T[lidx] + fn[lidx]);
		}
		else if (i == Nnode - 1) {
			mat_RHS[lidx] = T[lidx] - fn[lidx] + fact[lidx] * fn[lidx];
		}
		else {
			mat_RHS[lidx] = T[lidx] + fact[lidx] * (fn[lidx] - fn[lidx - 1]);
		}

	}
	/* compute soil temperatures */
	solve_temperature(mat_A, mat_B, mat_C,
					  mat_RHS, step_dt,
					  Nlayer, cell, energy, snow);
	return(0);
}

/******************************************************************************
* @brief    Calculate the saturated area and runoff
******************************************************************************/
void
solve_temperature(double   			*mat_A,
                  double   			*mat_B,
                  double   			*mat_C,
                  double  		    *mat_RHS,
                  double    		 step_dt,
                  size_t			 Nlayer,
				  cell_data_struct  *cell,
				  energy_bal_struct *energy,
				  snow_data_struct  *snow)
{
	extern option_struct       options;

    size_t i;
    size_t Nsnow = snow->Nsnow;
    size_t Nnode = options.Nnode;
    double *T = energy->T;
    double *pack_T = snow->pack_T;
    double *soil_T = cell->soil_T;

    for (i = 0; i < Nlayer; i++) {
        mat_A[i] = mat_A[i] * step_dt;
        mat_B[i] = 1.0 + mat_B[i] * step_dt;
        mat_C[i] = mat_C[i] * step_dt;
        mat_RHS[i] = mat_RHS[i] * step_dt;
    }

    tridiag(mat_A, mat_B, mat_C, mat_RHS, Nlayer);

    for (i = 0; i < Nlayer; i++) {
        T[i] += mat_RHS[i];
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
}