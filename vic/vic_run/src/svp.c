/******************************************************************************
* @section DESCRIPTION
*
* Calculate values related to the saturated vapor pressure.
******************************************************************************/

#include "vic_run.h"

/******************************************************************************
* @brief        This routine computes the saturated vapor pressure
*
* @note         Handbook of Hydrology eqn 4.2.2.
******************************************************************************/
void
svp_flags(double  TK,
          double  pressure,
          double *esat_T,
          double *qsat_T,
          double *esdT,
          double *qsdT,
          int     flags)
{
    extern parameters_struct param;

    double es_local = 0.0;
    double esdT_local = 0.0;
    double temp = K_TO_C(TK);
    temp = min(100.0, max(-75.0, temp));
    
    // 计算饱和水汽压（如果需要）
    if (flags & (ESAT | QSAT | QSDT)) {
        if (temp > 0.0) {
            es_local = param.SVP_A0 + temp * (param.SVP_A1 + temp * 
                          (param.SVP_A2 + temp * (param.SVP_A3 + temp * 
                          (param.SVP_A4 + temp * (param.SVP_A5 + temp * 
                          (param.SVP_A6 + temp * (param.SVP_A7 + temp * 
                          param.SVP_A8)))))));
        }
        else {
            es_local = param.SVP_B0 + temp * (param.SVP_B1 + temp * 
                          (param.SVP_B2 + temp * (param.SVP_B3 + temp * 
                          (param.SVP_B4 + temp * (param.SVP_B5 + temp * 
                          (param.SVP_B6 + temp * (param.SVP_B7 + temp * 
                          param.SVP_B8)))))));
        }
    }
    
    // 计算esat_T
    if ((flags & (ESAT | QSAT | QSDT)) && esat_T != NULL) {
        *esat_T = es_local * 100.0;
    }
    
    // 计算esdT（如果需要）
    if (flags & (ESDT | QSDT)) {
        if (temp > 0.0) {
            esdT_local = param.SVP_C0 + temp * (param.SVP_C1 + temp * 
                        (param.SVP_C2 + temp * (param.SVP_C3 + temp * 
                        (param.SVP_C4 + temp * (param.SVP_C5 + temp * 
                        (param.SVP_C6 + temp * (param.SVP_C7 + temp * 
                        param.SVP_C8)))))));
        }
        else {
            esdT_local = param.SVP_D0 + temp * (param.SVP_D1 + temp * 
                        (param.SVP_D2 + temp * (param.SVP_D3 + temp * 
                        (param.SVP_D4 + temp * (param.SVP_D5 + temp * 
                        (param.SVP_D6 + temp * (param.SVP_D7 + temp * 
                        param.SVP_D8)))))));
        }
    }
    if ((flags & (ESDT | QSDT)) && esdT != NULL) {
        *esdT = esdT_local * 100.0;
    } 
    // 计算qsat_T和qsdT（如果需要）
    if (flags & (QSAT | QSDT)) {
        double vp = 1.0 / (pressure - 0.378 * (es_local * 100.0));
        double vp1 = param.SVP_RDAIR * vp;
        
        if (flags & QSAT) {
            *qsat_T = es_local * 100.0 * vp1;
        }
        
        if (flags & QSDT) {
            double vp2 = vp1 * vp;
            *qsdT = esdT_local * 100.0 * vp2 * pressure;
        }
    }
}