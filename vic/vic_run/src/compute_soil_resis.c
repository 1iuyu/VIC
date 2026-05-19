/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate soil resistance (Ra_evap and soilbeta).
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief   Compute the resistance factor for soil evaporation calculation.
 *****************************************************************************/
void
compute_soil_resis(cell_data_struct *cell,
                   soil_con_struct  *soil_con)
{
    extern parameters_struct param;
    /* Allocate temp arrays */
    double *liq = cell->liq;
    double *porosity = cell->porosity;  // 有效孔隙度
    double *bexp_node = soil_con->bexp_node;
    double *Wsat_node = soil_con->Wsat_node;
    double *psisat_node = soil_con->psisat_node;

    double Ra_evap = 0.0;
    double aird = Wsat_node[0] * pow(psisat_node[0] / 1.e4, 1.0 / bexp_node[0]);
    double d0 = 2.12e-5 * pow(cell->soil_T[0] / CONST_TKFRZ, 1.75);
    double eps = Wsat_node[0] - aird;
    double dg = eps * d0 * pow(eps / Wsat_node[0], 3.0 / max(3.0, bexp_node[0]));
            
    double dsl = 15.0 * max(0.001, (0.8 * porosity[0] - liq[0])) /
                 max(0.001, (0.8 * Wsat_node[0] - aird));
            
    dsl = min(max(dsl,0.0), 200.0);
    Ra_evap = dsl / (dg * eps * 1.e3) + 20.0;
    Ra_evap = min(param.MAX_LIMIT, Ra_evap);
      
    cell->Ra_evap = Ra_evap;
}
