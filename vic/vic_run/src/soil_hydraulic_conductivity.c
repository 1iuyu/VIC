/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine calculate the hydraulic conductivity of the soil layer.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    Calculate the hydraulic conductivity of the soil layer.
 *****************************************************************************/
void
soil_hydraulic_conductivity(cell_data_struct *cell,
                            soil_con_struct  *soil_con)
{
    size_t Nsoil = cell->Nsoil;
    size_t i, j;
    double ice_param = -6.0;
    double *ice = cell->ice;
    double *liq = cell->liq;
    double frac_ice[MAX_SOILS];
    double frac_liq[MAX_SOILS];
    double *matric = cell->matric;
    double *zc_soil = soil_con->zc_soil;
    double *Zsum_soil = soil_con->Zsum_soil;
    double *Wsat_node = soil_con->Wsat_node;
    double *conductivity = cell->conductivity;
    double *conduct_int = cell->conduct_int;
    double *liquid_flux = cell->liquid_flux;
    double *soil_imped = cell->soil_imped;
    // 初始化变量
    for (j = 0; j < MAX_SOILS; j++) {
        frac_ice[j] = 0.0;
        frac_liq[j] = 0.0;
    }

    for (i = 0; i < Nsoil; i++) {
        //  计算土壤水和冰的相对饱和度
        frac_ice[i] = max(0.01, min(1.0, ice[i] / Wsat_node[i]));
        frac_liq[i] = max(0.01, min(1.0, liq[i] / Wsat_node[i]));
        // 计算冻土导致的不透水率
        if (i == Nsoil - 1) {
            soil_imped[i] = pow(10.0, frac_ice[i] * ice_param);
        }
        else {
            soil_imped[i] = pow(10.0, 0.5 * (frac_ice[i] + frac_ice[i+1]) * ice_param);
        }
        // 修正冻土导致的不透水性
        conductivity[i] = soil_imped[i] * SoilWaterRetentionCurve(CONDUCT_FLAG, i, 
                                                                  liq[i], 0.0, soil_con);
    }

    // Calculates the liquid moisture flux between soil nodes.
    for (i = 0; i < Nsoil - 1; i++) {
        double dzp = zc_soil[i+1] - zc_soil[i];
        double k_int = conductivity[i] * conductivity[i+1] * dzp / 
                                        (conductivity[i] * (zc_soil[i+1] - Zsum_soil[i]) +
                                         conductivity[i+1] * (Zsum_soil[i] - zc_soil[i])); // 调和平均
        conduct_int[i] = k_int / dzp;
        liquid_flux[i] = conduct_int[i] * (matric[i] - matric[i+1] + dzp);
    }
    conduct_int[Nsoil-1] = 0.0;
    liquid_flux[Nsoil-1] = 0.0;
}