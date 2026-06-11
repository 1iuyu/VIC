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
compute_soil_resis(double            coverage,
                   cell_data_struct *cell,
                   soil_con_struct  *soil_con)
{
    extern parameters_struct param;
    /* Allocate temp arrays */
    double Ra_evap = 0.0;
    if (cell->IS_VEG) {
        double *liq = cell->liq;
        double *porosity = cell->porosity;  // 有效孔隙度
        double *Wsat_node = soil_con->Wsat_node;
        double aird = SoilWaterRetentionCurve(MOIST_FLAG, 0, 0.0, -1.0e4, soil_con);
        double d0 = CONST_VAPDIFF * pow(cell->soil_T[0] / CONST_TKFRZ, 1.75);
        double eps = Wsat_node[0] - aird;
        double eps100 = Wsat_node[0] - SoilWaterRetentionCurve(MOIST_FLAG, 0, 0.0, -1.0, soil_con);
        // 3POE 方案的曲折度计算
        double tao = Wsat_node[0] * Wsat_node[0] * 
                    pow(eps / Wsat_node[0], 2.0 + log(pow(eps100, 0.25)) / log(eps100 / Wsat_node[0]));
        // 计算气体扩散度
        double dg = d0 * tao;
                
        double dsl = 0.015 * max(param.TOL_A, (0.8 * porosity[0] - liq[0])) /
                    max(param.TOL_A, (0.8 * Wsat_node[0] - aird));
                
        dsl = min(max(dsl, 0.0), 0.2);
        Ra_evap = dsl / (dg * eps * 1.e3) + 20.0;
        // 添加积雪覆盖修正
        if (1.0 - coverage + coverage * Ra_evap > 0.0) {
            Ra_evap = Ra_evap / (1.0 - coverage + coverage * Ra_evap);
        }
        else {
            Ra_evap = 0.0;
        }
        Ra_evap = min(param.MAX_LIMIT, Ra_evap);
    }
    else if (cell->IS_URBAN) {
        Ra_evap = param.MAX_LIMIT;
    }
    else {
        Ra_evap = 0.0;
    }
    cell->Ra_evap = Ra_evap;
}
