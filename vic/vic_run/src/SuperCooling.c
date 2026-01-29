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
SuperCooling(size_t           i,
             double           T,
             double           liq,
             double           moist,
             double          *supercool,
             soil_con_struct *soil_con)
{
    double      ice;
    double      delta_ice;

    double bexp_node = soil_con->bexp_node[i];
    double *Wsat_node = soil_con->Wsat_node;
    double *psi_sat_node = soil_con->psi_sat_node;
    if (bexp_node > 5.5) {
        bexp_node = 5.5;
    }
    double CK = 8.0;
    size_t iter = 0;
    size_t count = 0;
    if (T > CONST_TKFRZ - 1.0e-3) {
        *supercool = moist;
    }
    else {
        // initial guess for Ice (frozen content)
        if (CK != 0.0) {
            ice = moist - liq;
            if (ice > moist - 0.02) {
                ice = moist - 0.02;
            }
            // start the iterations
            if (ice < 0.0) {
                ice = 0.0;
            }
            while (iter < 10 && count == 0) {
                iter++;

                // 计算DF
                double DF = log(psi_sat_node[i] * CONST_G / CONST_LATICE * pow(1.0 + CK * ice, 2.0) * 
                    pow(Wsat_node[i] / (moist - ice), bexp_node)) - log(-(T - CONST_TKFRZ) / T);

                double Denom = 2.0 * CK / (1.0 + CK * ice) + bexp_node / (moist - ice);
                double tmp_ice = ice - DF / Denom;
                if (tmp_ice > moist - 0.02) {
                    tmp_ice = moist - 0.02;
                }
                if (tmp_ice < 0.0) {
                    tmp_ice = 0.0;
                }
                delta_ice = fabs(tmp_ice - ice);
                ice = tmp_ice;
                if (delta_ice < 0.005) {
                    count++;
                    break;
                }
            }
            // 如果超过10次迭代，使用显式方法
            if (count == 0) {
                log_warn("Flerchinger used in NEW version.");
                double FlerFac = pow((CONST_LATICE / (CONST_G * -psi_sat_node[i])) * 
                                    ((T - CONST_TKFRZ) / T), -1.0 / bexp_node) * Wsat_node[i];
                if (FlerFac < 0.02) {
                    FlerFac = 0.02;
                }
                (*supercool) = min(FlerFac, moist);
            }
            (*supercool) = moist - ice;
        }
    }
    
    return (0);
}
