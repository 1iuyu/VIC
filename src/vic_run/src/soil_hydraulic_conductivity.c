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
    size_t i, Nsoil = cell->Nsoil;
    double ice_param = -6.0;
    double *ice = cell->ice;
    double *liq = cell->liq;
    double frac_ice[MAX_SOILS] = {0};
    double *matric = cell->matric;
    double *Wpwp_node = soil_con->Wpwp_node;
    double *zc_soil = soil_con->zc_soil;
    double *Zsum_soil = soil_con->Zsum_soil;
    double *conductivity = cell->conductivity;
    double *conduct_int = cell->conduct_int;
    double *liquid_flux = cell->liquid_flux;
    double *soil_imped = cell->soil_imped;
    //  计算土壤水和冰的相对饱和度
    for (i = 0; i < Nsoil; i++) {
        double ice_equiv = ice[i] * (CONST_RHOICE / CONST_RHOFW);
        double denom = liq[i] + ice_equiv - Wpwp_node[i];
        if (denom > 1.0e-10) {  // 防止除零
            frac_ice[i] = min(1.0, ice_equiv / denom);
        }
        else {
            frac_ice[i] = 0.0;  // 极端干燥情况（无可流动水）
        }
    }
    // 计算冻土导致的不透水率
    for (i = 0; i < Nsoil; i++) {
        if (i == Nsoil - 1) {
            soil_imped[i] = pow(10.0, frac_ice[i] * ice_param);
        }
        else {
            soil_imped[i] = pow(10.0, 0.5 * (frac_ice[i] + frac_ice[i+1]) * ice_param);
        }
        // 修正冻土导致的不透水性
        conductivity[i] = soil_imped[i] * SoilWaterRetentionCurve(CONDUCT_FLAG, i, 
                                                                  0.0, matric[i], soil_con);
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