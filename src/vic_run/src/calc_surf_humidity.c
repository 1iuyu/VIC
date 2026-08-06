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
calc_surf_humidity(double            pressure,
                   double            Qair,
                   double            Tgrnd,
                   snow_data_struct *snow,
                   cell_data_struct *cell)
{
    double rh_grnd = 1.0;
    double qsat_soil = 0.0;
    double qsat_snow = 0.0;
    double qsdT_soil = 0.0;
    double qsdT_snow = 0.0;
    double alpha_soil = 0.0;
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
        rh_grnd = 1.0;
    }
    else if (cell->IS_GLAC || cell->IS_WET) {
        rh_grnd = 1.0;
    }
    // compute humidities individually for snow, soil for vegetated
    if (cell->IS_VEG) {
        svp_flags(soil_T[0], pressure, 
                NULL, &qsat_soil, 
                NULL, &qsdT_soil, QSAT | QSDT);

        if (qsat_soil > Qair && Qair > alpha_soil * qsat_soil) {
            qsat_soil = Qair;
            qsdT_soil = 0.0;
        }
        // soil humidity
        cell->Qair_soil = qsat_soil * alpha_soil;
        cell->QsdT_soil = qsdT_soil * alpha_soil;
        /* 2) snow 顶层饱和比湿 */
        if (snow->Nsnow > 0) {
            svp_flags(pack_T[0], pressure, 
                      NULL, &qsat_snow, 
                      NULL, &qsdT_snow, QSAT | QSDT);
            cell->QsdT_snow = qsdT_snow;     
            cell->QsdT_grnd = coverage * qsdT_snow + (1.0 - coverage) * alpha_soil * qsdT_soil;
            cell->Qair_snow = qsat_snow;
            cell->Qair_grnd = coverage * cell->Qair_snow + (1.0 - coverage) * cell->Qair_soil; 
        } else {
            cell->QsdT_snow = alpha_soil * qsdT_soil;
            cell->Qair_snow = cell->Qair_soil;
            cell->QsdT_grnd = alpha_soil * qsdT_soil;
            cell->Qair_grnd = cell->Qair_soil;
        }      
    }
    else {
        svp_flags(Tgrnd, pressure, 
                  NULL, &qsat_soil,
                  NULL, &qsdT_soil, QSAT | QSDT);
        cell->Qair_grnd = rh_grnd * qsat_soil;
        cell->QsdT_grnd = rh_grnd * qsdT_soil;
        if (qsat_soil > Qair && Qair > rh_grnd * qsat_soil) {
            cell->Qair_grnd = Qair;
            cell->QsdT_grnd = 0.0;
        }
        cell->Qair_soil = cell->Qair_grnd;
        cell->Qair_snow = cell->Qair_grnd;
        cell->QsdT_snow = cell->QsdT_grnd;
        cell->QsdT_soil = cell->QsdT_grnd;
    }
    
    return (0);
}