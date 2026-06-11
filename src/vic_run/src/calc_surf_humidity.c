/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes the lake variables for each new grid cell.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    This routine Calculate surface humidities, as well as a few 
 *           intermediate variables that are needed in the humidity calculations
 *****************************************************************************/
int
calc_surf_humidity(double             Tgrnd,
                   double             pressure,
                   double             Qair,
                   energy_bal_struct *energy,
                   snow_data_struct  *snow,
                   cell_data_struct  *cell,
                   soil_con_struct   *soil_con)
{
    double rh_grnd = 1.0;
    double qsat_Tgrnd = 0.0;
    double qsdT = 0.0;
    double alpha_soil = 0.0;
    double ice_factor = 0.0;
    double coverage = snow->coverage;
    double *soil_T = cell->soil_T;
    double *matric = cell->matric;
    double *pack_T = snow->pack_T;
    if (cell->IS_VEG) {
        // 计算降低地面饱和比湿度的因子（-）
        alpha_soil = exp(matric[0] / (CONST_RWV / CONST_G) / soil_T[0]);
        rh_grnd = (1.0 - coverage) * alpha_soil + coverage;
    }
    else if (cell->IS_URBAN) {

    }
    else if (cell->IS_GLAC || cell->IS_WET) {
        rh_grnd = 1.0;
    }
    // compute humidities individually for snow, soil for vegetated
    if (cell->IS_VEG) {
        svp_flags(Tgrnd, pressure, 
                NULL, &qsat_Tgrnd, 
                NULL, &qsdT, QSAT | QSDT);
        // soil humidity
        cell->Qair_soil = qsat_Tgrnd * alpha_soil;
        /* 2) snow 顶层饱和比湿 */
        if (snow->Nsnow > 0) {
            ice_factor = exp(CONST_LATICE * (pack_T[0] - CONST_TKFRZ) / 
                                        (CONST_RWV * pow(pack_T[0], 2.0)));
            cell->Qair_snow = qsat_Tgrnd * ice_factor;
        } else {
            cell->Qair_snow = cell->Qair_soil;
        }
        
        cell->Qair_grnd = coverage * cell->Qair_snow
                + (1.0 - coverage) * cell->Qair_soil;
        energy->qsdT = coverage * qsdT * ice_factor +
                        (1.0 - coverage) * alpha_soil * qsdT;        
    }
    else {
        
    }
    
    return (0);

}