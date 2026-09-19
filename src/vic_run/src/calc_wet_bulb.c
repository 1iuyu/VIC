/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine calculates the wet-bulb temperature based on the air 
 * temperature and relative humidity.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief  This subroutine calculates the wet-bulb temperature.
 *****************************************************************************/
double
calc_wet_bulb(double  air_temp,
              double  pressure,
              double  rel_humid)
{
    extern parameters_struct param;
    const double kappad = 0.2854;  // Heat Capacity
    const double WB_A = 2675.0; // Constant used for extreme cold temparatures (K)
    const double LAMBDA = 3.504;
    double foftk, fdT;
    double wb_temp, wb_temp_new;
    double rs_teq, dlnes_teq;
    double hot, cold;
    double wet_bulb = air_temp;
    bool converged = false;

    double pmb = pressure * 0.01;
    double pi = pow(pmb / CONST_PSREF, kappad);
    double temp_diff = air_temp - CONST_TKFRZ + 243.5;
    double esat_mb = param.SVP_FRZ * exp(17.67 * (air_temp - CONST_TKFRZ) / temp_diff);
    double adjust_pres = CONST_PSREF * pow(pi, LAMBDA);
    double humidity = param.SVP_RDAIR * esat_mb / (adjust_pres - esat_mb);
    double vapemb = esat_mb * rel_humid / 100.0;
    double mixr = humidity * rel_humid / 100.0 * GRAMS_PER_KG;

    double tl = 1.0 / (1.0 / (air_temp - 55.0) - log(rel_humid / 100.0) / 2840.0) + 55.0;
    double theta_dl = air_temp * pow(CONST_PSREF / (pmb - vapemb), kappad) * pow(air_temp / tl, mixr * 0.00028);
    double epott = theta_dl * exp((3.036 / tl - 0.00178) * mixr * (1.0 + 0.000448 * mixr));
    double Teq = epott * pi;
    double X = pow(CONST_TKFRZ / Teq, LAMBDA);
    double D = pow(0.1859 * pmb / CONST_PSREF + 0.6512, -1.0);
    double k1 = -38.5 * pi * pi + 137.81 * pi - 53.737;
    double k2 = -4.392 * pi * pi + 56.831 * pi - 0.384;
    // Calculates the regime requirements of wet bulb equations.
    if (Teq > 355.15) {
        hot = 1.0;
    } else {
        hot = 0.0;
    }
    if (X >= 1.0 && X <= D) {
        cold = 0.0;
    } else {
        cold = 1.0;
    }

    if (X > D) {
        wet_bulb_Qsat(Teq, pressure,
                     &dlnes_teq, &rs_teq, &foftk, &fdT);

        wb_temp = Teq - CONST_TKFRZ - WB_A * rs_teq / (1.0 + WB_A * rs_teq * dlnes_teq);
    }
    else {
        wb_temp = k1 - 1.21 * cold - 1.45 * hot - (k2 - 1.21 * cold) * X + (0.58 / X) * hot;
    }
    /* iterate wet-bulb temperature */
    size_t iter = 0;
    do {
        // 限制温度
        if (wb_temp > 50.0) {
            break;
        }
        wet_bulb_Qsat(wb_temp + CONST_TKFRZ, pressure, 
                      &dlnes_teq, &rs_teq, &foftk, &fdT);

        wb_temp_new = wb_temp - (foftk - X) / fdT;

        if (fabs(wb_temp - wb_temp_new) < param.TOL_WETBULB) {
            converged = true;
            wb_temp = wb_temp_new;
            break;
        }

        wb_temp = wb_temp_new;

        iter++;
 
    } while (iter < param.MAX_ITER_WETBULB);

    if (!converged) {
        wet_bulb = air_temp;
        log_warn("WARNING: Wet_Bulb algorithm failed to converge. Setting to air_temp");
    }
    else {
        wet_bulb = wb_temp + CONST_TKFRZ;
    }
    
    return (wet_bulb);
}

/******************************************************************************
* @brief    Calculate saturation mixing ratio and the change in saturation 
*           mixing ratio with respect to temperature.
******************************************************************************/
void 
wet_bulb_Qsat(double    air_temp, 
              double    pressure,
              double   *dlnes_mbdT,
              double   *humidity,
              double   *foftk,
              double   *fdT)
{
    extern parameters_struct param;
    /* 常量定义 */
    const double lambd_a = 3.504;      /* 热容量的倒数 */
    const double alpha   = 17.67;      /* 计算水汽压的常数 */
    const double beta    = 243.5;      /* 计算水汽压的常数 */
    const double vkp     = 0.2854;     /* 热容量 */
    const double y0      = 3036.0;     /* 常数 */
    const double y1      = 1.78;       /* 常数 */
    const double y2      = 0.448;      /* 常数 */

    /* Pa -> mb */
    double p_tmb = pressure * 0.01;
    double tcfbdiff = air_temp - CONST_TKFRZ + beta;
    double esat_mb = param.SVP_FRZ * exp(alpha * (air_temp - CONST_TKFRZ) / tcfbdiff);
    *dlnes_mbdT = alpha * beta / (tcfbdiff * tcfbdiff);
    double pminuse = p_tmb - esat_mb;
    double de_mbdT = esat_mb * (*dlnes_mbdT);
    double ndimpress = pow(p_tmb / CONST_PSREF, vkp);
    double p0ndplam = CONST_PSREF * pow(ndimpress, lambd_a);
    (*humidity) = param.SVP_RDAIR * esat_mb / (p0ndplam - esat_mb);
    double prersdt = param.SVP_RDAIR * p_tmb / (pminuse * pminuse);
    double rsdT = prersdt * de_mbdT;

    double rsy2rs2 = (*humidity) + y2 * (*humidity) * (*humidity);
    double oty2rs = 1.0 + 2.0 * y2 * (*humidity);
    double y0tky1 = y0 / air_temp - y1;
    double goftk = y0tky1 * ((*humidity) + y2 * (*humidity) * (*humidity));
    double gdT = -y0 * rsy2rs2 / (air_temp * air_temp) + y0tky1 * oty2rs * rsdT;

    (*foftk) = pow(CONST_TKFRZ / air_temp, lambd_a) * pow(1.0 - esat_mb / 
             p0ndplam, vkp * lambd_a) * exp(-lambd_a * goftk);

    (*fdT) = -lambd_a * (1.0 / air_temp + vkp * de_mbdT / pminuse + gdT) * (*foftk);
}