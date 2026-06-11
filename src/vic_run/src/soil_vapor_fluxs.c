/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine calculate the vapor fluxes between soil nodes due to changes
 * in potential gradient and temperature gradient. (The fluxes resulting from 
 * the temperature gradient will be multiplied by an enhancement factor.)
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    Calculate the vapor fluxes between soil nodes.
 *****************************************************************************/
int
calc_vapor_flux(double             pressure,
                cell_data_struct  *cell,
                energy_bal_struct *energy,
                snow_data_struct  *snow,
                soil_con_struct   *soil_con)
{
    /* Initialize variables */
    size_t i, lidx;
    double air = 0.0;
    double vapor_diff = 0.0;
    double coef_vapor = 0.6;
    double vapor_exp  = 1.0;
    double enhance = 0.0;
    double ice_corr = 0.0;
    double qsdT = 0.0;
    double qsat_T = 0.0;
    double esat_T = 0.0;
    double air_density = 0.0;
    double enhanc_fact1 = 9.5;
    double enhanc_fact2 = 3.0;
    double enhanc_fact3 = 0.0;
    double enhanc_fact4 = 1.0;
    double enhanc_fact5 = 4.0;
    double rel_humid[MAX_NODES];
    double diff_therm[MAX_NODES];
    double diff_vapor[MAX_NODES];
    double conv_temp[MAX_NODES];
    double conv_vapor[MAX_NODES];
    double *ice = cell->ice;
    double *liq = cell->liq;
    double *T = energy->T;
    double *Wsat_node = soil_con->Wsat_node;
    double *zc_soil = soil_con->zc_soil;
    double *Zsum_soil = soil_con->Zsum_soil;
    double *matric = cell->matric;
    double *vapor_flux = cell->vapor_flux;
    double *clay_node = soil_con->clay_node;
    double *deric_vapor = cell->deriv_vapor;
    // 初始化变量
    for (i = 0; i < MAX_NODES; i++) {
        rel_humid[i] = 0.0;
        diff_therm[i] = 0.0;
        diff_vapor[i] = 0.0;
        conv_temp[i] = 0.0;
        conv_vapor[i] = 0.0;
    }
    size_t Nsnow = snow->Nsnow;
    for (i = 0; i < Nsnow; i++) {
        vapor_diff = 0.00009 * (CONST_PSTD / pressure) *
                         pow(1.0 + T[i] / CONST_TKFRZ, 14.0);

        svp_flags(T[i], pressure, 
                  NULL, &qsat_T, 
                  NULL, &qsdT, 
                  QSAT | QSDT);
        air_density = pressure / (CONST_RDAIR * T[i]);
        ice_corr = exp(CONST_MWWV * CONST_LATICE * (T[i] - CONST_TKFRZ)
                                       / (CONST_RGAS * T[i] * T[i]));
        diff_vapor[i] = qsat_T * air_density * ice_corr;
        diff_therm[i] = qsdT * air_density * ice_corr;
    }

    size_t Nsoil = cell->Nsoil;
    for (i = 0; i < Nsoil; i++) {
        lidx = i + Nsnow;
        air = Wsat_node[i] - ice[i] - liq[i];
        if (air > 0.0 && matric[i] < 0.0) {
            // Calculate vapor fluxes due to temperature gradient
            vapor_diff = CONST_VAPDIFF *
                pow(T[lidx] / CONST_TKFRZ, 2.0) *
                (CONST_PSTD / pressure);
            vapor_diff *= coef_vapor * pow(air, vapor_exp);

            svp_flags(T[lidx], pressure, 
                      &esat_T, &qsat_T, 
                      NULL, &qsdT, 
                      ESAT | QSAT | QSDT);

            // Calculate relative humidity
            rel_humid[lidx] = exp(CONST_MWWV * CONST_G / CONST_RGAS / T[lidx] * matric[i]);

            // 计算实际水汽压
            double e_actual = esat_T * rel_humid[lidx];
            
            // 计算土层空气密度
            air_density = (pressure - 0.378 * e_actual) / (CONST_RDAIR * T[lidx]);

            // Calculate vapor fluxes due to potential gradient
            if (clay_node[i] <= 0.02) {
                enhanc_fact3 = Wsat_node[i] * (1.0 + 2.6 / sqrt(0.02));
            } else {
                enhanc_fact3 = Wsat_node[i] * (1.0 + 2.6 / sqrt(clay_node[i]));
            }
            // 计算指数项
            double expon = -pow(enhanc_fact3 * liq[i] / Wsat_node[i], enhanc_fact5);
            if (expon <= -50.0) {
                expon = 0.0;
            } else {
                expon = exp(expon);
            }
            enhance = enhanc_fact1 + enhanc_fact2 * liq[i] / Wsat_node[i] 
                        - (enhanc_fact1 - enhanc_fact4) * expon;

            diff_therm[i] = vapor_diff * enhance * rel_humid[lidx] * qsdT * air_density;
            diff_vapor[i] = vapor_diff * qsat_T * air_density;
        }
        else {
            diff_therm[i] = 0.0;
            diff_vapor[i] = 0.0;
            rel_humid[lidx] = 1.0;
        }
    }
    
    // 计算层间调和平均传导系数
    for (i = 0; i < Nsoil - 1; i++) {
        double dzp = zc_soil[i+1] - zc_soil[i];

        conv_temp[i] = diff_therm[i] * diff_therm[i+1] * dzp /
                        (diff_therm[i] * (zc_soil[i+1] - Zsum_soil[i]) +
                         diff_therm[i+1] * (Zsum_soil[i] - zc_soil[i]));
        conv_vapor[i] = diff_vapor[i] * diff_vapor[i+1] * dzp /
                        (diff_vapor[i] * (zc_soil[i+1] - Zsum_soil[i]) +
                         diff_vapor[i+1] * (Zsum_soil[i] - zc_soil[i]));
        
        // 水汽通量 = 热梯度项 + 湿度梯度项
        vapor_flux[i] = conv_temp[i] * (T[i] - T[i+1]) / dzp + 
                            conv_vapor[i] * (rel_humid[i] - rel_humid[i+1]) / dzp;
        deric_vapor[i] = CONST_MWWV * CONST_G / CONST_RGAS / T[i] * conv_vapor[i] / dzp;
    }

    return (0);
}