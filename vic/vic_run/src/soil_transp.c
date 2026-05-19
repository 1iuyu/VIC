/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate the soil layer transpiration factor.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief  
 * Compute the soil layer transpiration factor.
 *****************************************************************************/
int
soil_transp(cell_data_struct *cell,
			soil_con_struct  *soil_con)
{
    extern parameters_struct param;
    size_t i;
    size_t Nroot = cell->Nroot;
    double psi_wp = param.SRESP_PSIWILT;
    double total_transp = 0.;
    double wet_fact = 0.0;
    double psi_soil = 0.0;
    double *transp_fact = cell->transp_fact;
    double *liq = cell->liq;
    double *Zsum_soil = soil_con->Zsum_soil;
    double *dz_soil = soil_con->dz_soil;
    double *Wsat_node = soil_con->Wsat_node;
    double *bexp_node = soil_con->bexp_node;
    double *psisat_node = soil_con->psisat_node;

    for (i = 0; i < Nroot; i++) {
        psi_soil = max(psi_wp, -psisat_node[i] * 
                            pow(max(0.01, liq[i]) /
                                    Wsat_node[i], -bexp_node[i]));
        wet_fact = (1.0 - psi_soil / psi_wp) / 
                        (1.0 + psisat_node[i] / psi_wp);

        wet_fact = min(1.0, max(0.0, wet_fact));

        transp_fact[i] = max(param.TOL_A, dz_soil[i] / Zsum_soil[Nroot] * wet_fact);

        total_transp += transp_fact[i];
    }

    total_transp = min(max(total_transp, param.TOL_A), 1.0);
    cell->total_transp = total_transp;
    for (i = 0; i < Nroot; i++) {
        transp_fact[i] /= total_transp;
    }
    
    return(0);
}