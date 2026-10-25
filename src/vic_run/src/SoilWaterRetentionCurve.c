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
                        double           soil_T,
                        soil_con_struct *soil_con)
{
    extern parameters_struct param;
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
            return -1e4; // 防止数值爆炸
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
    else if (flag == DERIVW_FLAG) {
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
    else if (flag == DERIVM_FLAG) {
        /* dψ/dθ */
        if (liq < Wsat_node[idx]) {
            Se = (liq - Wpwp_node[idx]) / (Wsat_node[idx] - Wpwp_node[idx]);
            if (Se < MIN_SOILMOIST) {
                Se = MIN_SOILMOIST;
            }
            double y1 = pow(Se, -1.0 / mpar_node[idx]) - 1.0;
            double d1 = (-1.0 / mpar_node[idx]) * pow(Se, -1.0 / mpar_node[idx] - 1.0) / 
                        (Wsat_node[idx] - Wpwp_node[idx]);
            if (y1 < param.TOL_A) {
                y1 = param.TOL_A;
            }
            // double y2 = pow(y1, 1.0 / expt_node[idx]);
            double d2 = (1.0 / expt_node[idx]) * pow(y1, 1.0 / expt_node[idx] - 1.0);
            return d1 * d2 / alpha_node[idx];
        } 
        else {
            return (0.0);
        }
    }
    else if (flag == DERIVK_FLAG) {
        /* dK/dψ */
        if (matric >= 0.0) {
            return 0.0;
        }
        double f_x1 = pow(matric * alpha_node[idx], expt_node[idx] - 1.0);
        double f_x2 = pow(1.0 + pow(matric * alpha_node[idx], expt_node[idx]), -mpar_node[idx]);
        double d_x1 = alpha_node[idx] * (expt_node[idx] - 1.0) * pow(matric * alpha_node[idx], expt_node[idx] - 2.0);
        double d_x2 = alpha_node[idx] * expt_node[idx] * pow(matric * alpha_node[idx], expt_node[idx] - 1.0) *
              (-mpar_node[idx]) * pow(1.0 + pow(matric * alpha_node[idx], expt_node[idx]), -mpar_node[idx] - 1.0);

        double f_nm = (1.0 - f_x1 * f_x2) * (1.0 - f_x1 * f_x2);
        double d_nm = (-d_x1 * f_x2 - f_x1 * d_x2) * 2.0 * (1.0 - f_x1 * f_x2);

        double f_dm = pow(1.0 + pow(matric * alpha_node[idx], expt_node[idx]), mpar_node[idx] / 2.0);
        double d_dm = alpha_node[idx] * expt_node[idx] * pow(matric * alpha_node[idx], expt_node[idx] - 1.0) *
                        (mpar_node[idx] / 2.0) * pow(1.0 + pow(matric * alpha_node[idx], expt_node[idx]), 
                        mpar_node[idx] / 2.0 - 1.0);

        double orig_val = Ksat_node[idx] * (f_nm / f_dm);
        double orig_d = Ksat_node[idx] * (d_nm * f_dm - d_dm * f_nm) / (f_dm * f_dm);

        if (matric < -param.TOL_A) {
            return orig_d;
        } 
        else {
            double t = (matric + param.TOL_A) / param.TOL_A;
            if (t < 0.0) {
                t = 0.0;
            }
            if (t > 1.0) {
                t = 1.0;
            }
            double s = 3.0 * t * t - 2.0 * t * t * t;
            double ds_dpsi = (6.0 * t * (1.0 - t)) / param.TOL_A;
            return orig_d * (1.0 - s) + (Ksat_node[idx] - orig_val) * ds_dpsi;
        }
    }
    else if (flag == DERIVT_FLAG) {
        /* dθ/dT */
        if (soil_T > CONST_TKFRZ) {
            return 0.0;
        }
        else {
            double kappa = CONST_LATICE / (CONST_G * CONST_TKFRZ);
            double xtemp = alpha_node[idx] * kappa * (soil_T - CONST_TKFRZ);
            /* 链式法则求导 */
            return (alpha_node[idx] * kappa) * expt_node[idx] * pow(xtemp, expt_node[idx] - 1.0) *
                (-mpar_node[idx]) * pow(1.0 + pow(xtemp, expt_node[idx]), -mpar_node[idx] - 1.0) *
                (Wsat_node[idx] - Wpwp_node[idx]);
        }
    }
    else if (flag == CONDUCT_FLAG) {
        double *lpar_node = soil_con->lpar_node;
        if (matric >= 0.0) {
            return Ksat_node[idx];
        }
        // Van Genuchten-Mualem 公式
        double alpha_h = fabs(alpha_node[idx] * matric);
        double alpha_pow = pow(alpha_h, expt_node[idx]);
        double denom = pow(1.0 + alpha_pow, lpar_node[idx] * mpar_node[idx]);
        double term = 1.0 - 1.0 / (1.0 + alpha_pow);
        double term_pow = pow(term, mpar_node[idx]);
        double numera = 1.0 - term_pow;
        double numera_pow = numera * numera;
        return Ksat_node[idx] * numera_pow / denom;
    }
    else {
        log_err("Unknown FLAG option");
    }

    return 0.0;
}