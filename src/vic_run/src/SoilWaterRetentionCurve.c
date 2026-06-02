/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine defines the relationship between the volume liquid content
 * and the matric potential. 
 *  [1] Calculates the matric potential based on the moisture content.
 *  [2] Calculates the moisture content based on the matric potential.
 *  [3] Calculates the derivative of the matric potential with respect to
 *          the moisture content.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief  
 * This subroutine defines the relationship between the volume liquid content
 * and the matric potential.
 *****************************************************************************/
double
SoilWaterRetentionCurve(int              flag,
                        size_t           idx,
                        double           liq,
                        double           matric,
                        soil_con_struct *soil_con)
{
    extern option_struct options;
    double Se, tmp;
    double *Wpwp_node = soil_con->Wpwp_node;
    double *Wsat_node = soil_con->Wsat_node;
    double *bexp_node = soil_con->bexp_node;
    double *Ksat_node = soil_con->Ksat_node;
    double *psisat_node = soil_con->psisat_node;
    /*================================
        CAMPBELL 模型
    =================================*/
    if (options.SWRC == SWRC_CAMPBELL) {

        if (flag == MATRIC_FLAG) {
            /* θ → ψ */
            if (liq < Wsat_node[idx]) {
                return -psisat_node[idx] *
                            pow((liq / Wsat_node[idx]), -bexp_node[idx]);
            }
            else {
                return psisat_node[idx];
            }
        }
        else if (flag == MOIST_FLAG) {
            /* ψ → θ */
            if (matric < 0.0) {
                return Wsat_node[idx] *
                            pow((-matric / psisat_node[idx]), -1.0 / bexp_node[idx]);
            }
            else {
                return Wsat_node[idx];
            }
        }
        else if (flag == DERIV_FLAG) {
            /* dψ/dθ */
            if (matric < psisat_node[idx]) {
                return -liq / (bexp_node[idx] * matric);
            }
            else {
                return 0.0;
            }
        }
        else if (flag == CONDUCT_FLAG) {
            if (liq < Wsat_node[idx]) {
                return Ksat_node[idx] * pow(liq / Wsat_node[idx], 
                                        2.0 * bexp_node[idx] + 3.0);
            }
            else {
                return Ksat_node[idx];
            }
        }
        else {
            log_err("Unknown FLAG option");
        }
    }

    /* ================================
        VAN GENUCHTEN 模型
    ================================= */
    else if (options.SWRC == SWRC_VAN_GENUCHTEN) {

        double *alpha_node = soil_con->alpha_node;
        double m = 1.0 - 1.0 / bexp_node[idx];

        if (flag == MATRIC_FLAG) {
            /* θ → ψ */
            Se = (liq - Wpwp_node[idx]) /
                    (Wsat_node[idx] - Wpwp_node[idx]);

            if (Se >= 1.0) {
                return 0.0;
            }
            else if (Se <= 0.0) {
                return -1e6; // 防止数值爆炸
            }
            else {
                return -1.0 / fabs(alpha_node[idx]) *
                    pow(pow(Se, -1.0 / m) - 1.0,
                        1.0 / bexp_node[idx]);
            }
        }
        else if (flag == MOIST_FLAG) {
            /* ψ → θ */
            if (matric < psisat_node[idx]) {
                tmp = pow(1.0 + pow(fabs(alpha_node[idx] * matric),
                                    bexp_node[idx]), m);
                return Wpwp_node[idx] +
                            (Wsat_node[idx] - Wpwp_node[idx]) / tmp;
            }
            else {
                return Wsat_node[idx];
            }
        }
        else if (flag == DERIV_FLAG) {
            /* dψ/dθ */
            if (matric < psisat_node[idx]) {
                tmp = pow(fabs(alpha_node[idx] * 
                            matric), bexp_node[idx]);
                return -(liq - Wpwp_node[idx]) * 
                        m * bexp_node[idx] * tmp / (matric * (1.0 + tmp));
            }
            else {
                return 0.0;
            }
        }
        else if (flag == CONDUCT_FLAG) {
            double *lpar_node = soil_con->lpar_node;
            Se = (liq - Wpwp_node[idx]) / (Wsat_node[idx] - Wpwp_node[idx]);
            Se = max(0.0, min(1.0, Se));
            if (Se >= 1.0) {
                return Ksat_node[idx];
            }
            else {
                // Van Genuchten-Mualem 公式
                // K = Ks * Se^L * [1 - (1 - Se^(1/m))^m]^2
                double Se_pow = pow(Se, 1.0 / m);  // Se^(1/m)
                double term = 1.0 - pow(1.0 - Se_pow, m);  // [1 - (1-Se^(1/m))^m]
                return Ksat_node[idx] * pow(Se, lpar_node[idx]) * term * term;
            }
        }
        else {
            log_err("Unknown FLAG option");
        }
    }
    else {
        log_err("Unknown Soil Water Retention Curve FLAG option");
    }

    return 0.0;
}