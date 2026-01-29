/******************************************************************************
* @section DESCRIPTION
*
* Calculate values related to the saturated vapor pressure.
******************************************************************************/

#include <vic_run.h>

/******************************************************************************
* @brief        This routine computes the saturated vapor pressure
*
* @note         Handbook of Hydrology eqn 4.2.2.
******************************************************************************/
void
svp(double  temp,
    double *SVP_liq,
    double *SVP_ice)
{
    extern parameters_struct param;
    
    temp = K_TO_C(temp);
    *SVP_liq = 0.;
    *SVP_ice = 0.;
    
    (*SVP_liq) = 100 * (param.SVP_A0 + temp * (param.SVP_A1 + temp * 
                      (param.SVP_A2 + temp * (param.SVP_A3 + temp * 
                      (param.SVP_A4 + temp * (param.SVP_A5 + temp * 
                                              param.SVP_A6))))));
    (*SVP_ice) = 100 * (param.SVP_B0 + temp * (param.SVP_B1 + temp * 
                      (param.SVP_B2 + temp * (param.SVP_B3 + temp * 
                      (param.SVP_B4 + temp * (param.SVP_B5 + temp * 
                                              param.SVP_B6))))));  

}

/******************************************************************************
* @brief        This routine computes the gradient of d(svp)/dT
*
* @note         Handbook of Hydrology eqn 4.2.3
******************************************************************************/
void
svp_slope(double  temp,
          double *liq_slope,
          double *ice_slope)
{
    extern parameters_struct param;

    temp = K_TO_C(temp);
    *liq_slope = 0.;
    *ice_slope = 0.;

    (*liq_slope) = 100 * (param.SVP_C0 + temp * (param.SVP_C1 + temp * 
                      (param.SVP_C2 + temp * (param.SVP_C3 + temp * 
                      (param.SVP_C4 + temp * (param.SVP_C5 + temp * 
                                              param.SVP_C6))))));
    (*ice_slope) = 100 * (param.SVP_D0 + temp * (param.SVP_D1 + temp * 
                      (param.SVP_D2 + temp * (param.SVP_D3 + temp * 
                      (param.SVP_D4 + temp * (param.SVP_D5 + temp * 
                                              param.SVP_D6))))));
}
/******************************************************************************
* @brief  This routine computes the saturated surface specific humidity 
*         and changing rate to temperature.
******************************************************************************/
void
s_humid(double  temp,
        double  pressure,
        double *ratio,
        double *ratio_slope)
{
    extern parameters_struct param;

    double tmp_esat = param.SVP_FRZ * exp(CONST_LATVAP / CONST_RWV * 
                                (1.0 / CONST_TKFRZ - 1.0 / temp));
    (*ratio) = param.SVP_RDAIR * tmp_esat / (pressure / PA_PER_KPA - tmp_esat);

    /* convert from  g/g to g/kg */
    (*ratio) *= PA_PER_KPA;
    (*ratio_slope) = (*ratio / (1 + *ratio)) * param.SVP3 / pow(temp - param.SVP2, 2.0);

    /* convert from  g/kg to g/g */
    (*ratio) /= PA_PER_KPA;
}