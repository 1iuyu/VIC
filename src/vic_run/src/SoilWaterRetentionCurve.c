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
    double Se, tmp;
    double *Wpwp_node = soil_con->Wpwp_node;
    double *Wsat_node = soil_con->Wsat_node;
    double *expt_node = soil_con->expt_node;
    double *Ksat_node = soil_con->Ksat_node;

    /* ================================
        VAN GENUCHTEN 模型
    ================================= */
    double *alpha_node = soil_con->alpha_node;
    double *mpar_node = soil_con->mpar_node;

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
                pow(pow(Se, -1.0 / mpar_node[idx]) - 1.0,
                    1.0 / expt_node[idx]);
        }
    }
    else if (flag == MOIST_FLAG) {
        /* ψ → θ */
        if (matric < 0.0) {
            tmp = pow(1.0 + pow(fabs(alpha_node[idx] * matric),
                                expt_node[idx]), mpar_node[idx]);
            return Wpwp_node[idx] +
                        (Wsat_node[idx] - Wpwp_node[idx]) / tmp;
        }
        else {
            return Wsat_node[idx];
        }
    }
    else if (flag == DERIV_FLAG) {
        /* dθ/dψ */
        if (matric < 0.0) {
            tmp = pow(fabs(alpha_node[idx] * 
                        matric), expt_node[idx]);
            return -(liq - Wpwp_node[idx]) * 
                    mpar_node[idx] * expt_node[idx] * tmp / (matric * (1.0 + tmp));
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
            double Se_pow = pow(Se, 1.0 / mpar_node[idx]);
            double term = 1.0 - pow(1.0 - Se_pow, mpar_node[idx]);
            return Ksat_node[idx] * pow(Se, lpar_node[idx]) * term * term;
        }
    }
    else {
        log_err("Unknown FLAG option");
    }

    return 0.0;
}