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
                snow_data_struct  *snow,
                soil_con_struct   *soil_con)
{
    extern parameters_struct param;
    /* Initialize variables */
    size_t i, lidx;
    double air = 0.0;
    double dzp = 0.0;
    double coef_vapor = 0.6;
    double vapor_exp  = 1.0;
    double enhance = 0.0;
    double dair_dT = 0.0;
    double qsdT = 0.0;
    double qsaT = 0.0;
    double esaT = 0.0;
    double conv_temp = 0.0;
    double air_density = 0.0;
    double enhanc_fact = 0.0;
    double rel_humid[MAX_SOILS] = {0};
    double vapor_diff[MAX_NODES] = {0};
    double diff_therm[MAX_SOILS] = {0};
    double diff_vapor[MAX_NODES] = {0};
    double *ice = cell->ice;
    double *liq = cell->liq;
    double *pack_T = snow->pack_T;
    double *soil_T = cell->soil_T;
    double *zc_snow = snow->zc_snow;
    double *Zsum_snow = snow->Zsum_snow;
    double *Wsat_node = soil_con->Wsat_node;
    double *zc_soil = soil_con->zc_soil;
    double *Zsum_soil = soil_con->Zsum_soil;
    double *matric = cell->matric;
    double *vapor_flux = cell->vapor_flux;
    double *clay_node = soil_con->clay_node;
    double *conv_vapor = cell->conv_vapor;
    double *drhodT = cell->drhodT;
    double *dQvdSMP = cell->dQvdSMP;
    // 初始化变量
    size_t Nsnow = snow->Nsnow;
    size_t tmp_Nsnow = Nsnow;
    // 雪层水汽扩散
    if (Nsnow > 0) {
        for (i = 0; i < Nsnow; i++) {
            vapor_diff[i] = 0.00009 * (CONST_PSTD / pressure) *
                            pow(pack_T[i] / CONST_TKFRZ, 14.0);
            // 计算冰面饱和比湿
            svp_flags(pack_T[i], pressure, 
                      NULL, &qsaT, 
                      NULL, &qsdT, QSAT | QSDT);
            air_density = pressure / (CONST_RDAIR * pack_T[i]);
            dair_dT = -air_density / pack_T[i];
            diff_vapor[i] = qsaT * air_density;
            drhodT[i] = qsdT * air_density + qsaT * dair_dT;
        }
    }

    // 表层水汽扩散
    if (cell->h2osfc > param.TOL_A) {
        double h2osfc_T = cell->h2osfc_T;
        svp_flags(h2osfc_T, pressure, 
                  NULL, &qsaT, 
                  NULL, &qsdT, QSAT | QSDT);
        air_density = pressure / (CONST_RDAIR * h2osfc_T);
        if (cell->IS_GLAC) {
            diff_vapor[Nsnow] = qsaT * air_density;
            dair_dT = -air_density / h2osfc_T;
            drhodT[Nsnow] = qsdT * air_density + qsaT * dair_dT; 
            vapor_diff[Nsnow] = 0.00009 * (CONST_PSTD / pressure) *
                                pow(h2osfc_T / CONST_TKFRZ, 14.0);
        }
        else if (cell->IS_WET) {
            diff_vapor[Nsnow] = qsaT * air_density;
            drhodT[Nsnow] = qsdT * air_density + qsaT * dair_dT;
            vapor_diff[Nsnow] = CONST_VAPDIFF * pow(h2osfc_T / CONST_TKFRZ, 2.0) *
                                (CONST_PSTD / pressure); // 水层：使用自由水面扩散系数
        }
        tmp_Nsnow++;
    }
    else {
        svp_flags(soil_T[0], pressure, 
                  NULL, &qsaT, 
                  NULL, &qsdT, QSAT | QSDT);
        air_density = pressure / (CONST_RDAIR * soil_T[0]);
        double relh = exp(CONST_MWWV * CONST_G * matric[0] / 
                         (CONST_RGAS * soil_T[0]));
        drhodT[Nsnow] = qsdT * air_density + qsaT * dair_dT * relh;
    }
    // 土层水汽扩散
    size_t Nsoil = cell->Nsoil;
    for (i = 0; i < Nsoil; i++) {
        lidx = i + tmp_Nsnow;
        air = Wsat_node[i] - ice[i] - liq[i];
        if (air > 0.0 && matric[i] < 0.0) {
            // Calculate vapor fluxes due to temperature gradient
            vapor_diff[lidx] = CONST_VAPDIFF *
                pow(soil_T[i] / CONST_TKFRZ, 2.0) *
                (CONST_PSTD / pressure);
            vapor_diff[lidx] *= coef_vapor * pow(air, vapor_exp);

            svp_flags(soil_T[i], pressure,
                      &esaT, &qsaT, 
                      NULL, &qsdT, 
                      ESAT | QSAT | QSDT);

            // Calculate relative humidity
            rel_humid[i] = exp(CONST_MWWV * CONST_G / CONST_RGAS / soil_T[i] * matric[i]);
            
            // 计算土层空气密度
            air_density = pressure / (CONST_RDAIR * soil_T[i]);

            // 饱和水汽密度（kg/m³）
            double sat_vap_dens = qsaT * air_density;
            
            // 饱和水汽密度的温度导数（kg/m³/K）
            double dair_dT = -air_density / soil_T[i];
            double sat_vap_dens_dT = qsdT * air_density + qsaT * dair_dT;

            // Calculate vapor fluxes due to potential gradient
            if (clay_node[i] <= 0.02) {
                enhanc_fact = Wsat_node[i] * (1.0 + 2.6 / sqrt(0.02));
            } else {
                enhanc_fact = Wsat_node[i] * (1.0 + 2.6 / sqrt(clay_node[i]));
            }
            // 计算指数项
            double expon = -pow(enhanc_fact * liq[i] / Wsat_node[i], 4.0);
            if (expon <= -50.0) {
                expon = 0.0;
            } else {
                expon = exp(expon);
            }
            enhance = 9.5 + 3.0 * liq[i] / Wsat_node[i] - 8.5 * expon;

            diff_therm[i] = vapor_diff[lidx] * enhance * rel_humid[i] * sat_vap_dens_dT;
            diff_vapor[lidx] = vapor_diff[lidx] * sat_vap_dens;
        }
        else {
            diff_therm[i] = 0.0;
            diff_vapor[lidx] = 0.0;
            rel_humid[i] = 1.0;
        }
    }
    
    // 计算层间调和平均传导系数
    for (i = 0; i < Nsnow; i++) {
        if (i == Nsnow - 1) {
            if (cell->h2osfc > param.TOL_A) {
                dzp = zc_snow[i] + 0.5 * cell->h2osfc;
                conv_vapor[i] = vapor_diff[i] * vapor_diff[i+1] * dzp /
                    (vapor_diff[i] * 0.5 * cell->h2osfc + vapor_diff[i+1] * zc_snow[i]);
                conv_vapor[i] /= dzp;
            } 
            else {
                dzp = zc_soil[0] + zc_snow[i];
                conv_vapor[i] = vapor_diff[i] * vapor_diff[i+1] * dzp /
                    (vapor_diff[i] * zc_soil[0] +
                    vapor_diff[i+1] * zc_snow[i]);
                conv_vapor[i] /= dzp;
            }
        } 
        else {
            dzp = zc_snow[i+1] - zc_snow[i];
            conv_vapor[i] = vapor_diff[i] * vapor_diff[i+1] * dzp /
                (vapor_diff[i] * (zc_snow[i+1] - Zsum_snow[i]) +
                vapor_diff[i+1] * (Zsum_snow[i] - zc_snow[i]));
            conv_vapor[i] /= dzp;
        }
        vapor_flux[i] = conv_vapor[i] * (diff_vapor[i] - diff_vapor[i+1]);
    }
    if (cell->h2osfc > param.TOL_A) {
        dzp = zc_soil[0] + 0.5 * cell->h2osfc;
        conv_vapor[Nsnow] = vapor_diff[Nsnow] * vapor_diff[Nsnow+1] * dzp /
            (vapor_diff[Nsnow] * zc_soil[0] + vapor_diff[Nsnow+1] * 0.5 * cell->h2osfc);
        conv_vapor[Nsnow] /= dzp;
        vapor_flux[Nsnow] = conv_vapor[Nsnow] * (diff_vapor[Nsnow] - diff_vapor[Nsnow+1]);
    }
    for (i = 0; i < Nsoil - 1; i++) {
        lidx = tmp_Nsnow + i;
        dzp = zc_soil[i+1] - zc_soil[i];
        conv_temp = diff_therm[i] * diff_therm[i+1] * dzp /
                        (diff_therm[i] * (zc_soil[i+1] - Zsum_soil[i]) +
                        diff_therm[i+1] * (Zsum_soil[i] - zc_soil[i]));
        conv_temp /= dzp;
        conv_vapor[lidx] = diff_vapor[lidx] * diff_vapor[lidx+1] * dzp /
                        (diff_vapor[lidx] * (zc_soil[i+1] - Zsum_soil[i]) +
                        diff_vapor[lidx+1] * (Zsum_soil[i] - zc_soil[i]));
        conv_vapor[lidx] /= dzp;
        // 水汽通量 = 热梯度项 + 湿度梯度项
        vapor_flux[lidx] = conv_temp * (soil_T[i] - soil_T[i+1]) + 
                            conv_vapor[lidx] * (rel_humid[i] - rel_humid[i+1]);
        dQvdSMP[i] = CONST_MWWV * CONST_G / CONST_RGAS / soil_T[i] * conv_vapor[lidx];
    }
    dQvdSMP[Nsoil-1] = CONST_MWWV * CONST_G / CONST_RGAS / soil_T[Nsoil-1] *
                                     conv_vapor[tmp_Nsnow + Nsoil - 1];

    return (0);
}