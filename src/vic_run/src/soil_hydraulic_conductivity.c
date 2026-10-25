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
soil_hydraulic_conductivity(energy_bal_struct *energy,
                            cell_data_struct  *cell,
                            snow_data_struct  *snow,
                            soil_con_struct   *soil_con)
{
    extern parameters_struct param;
    extern option_struct     options;

    size_t i, Nsoil = cell->Nsoil;
    size_t Nsnow = snow->Nsnow;
    double ice_param = -6.0;
    double liquid_capacity;
    double conduct_mp;
    double dimped_dT;
    double *dpsiliq_dT = cell->dpsiliq_dT;
    double dTK_dpsi_noice;
    double dmacro_dliq, dmacro_dpsi, dmacro_dT;
    double dmicro_dT, dmicro_dpsi;
    double *ice = energy->ice;
    double *liq = energy->liq;
    double *soil_T = cell->soil_T;
    double *pack_liq = snow->pack_liq;
    double *pack_ice = snow->pack_ice;
    double *theta_ice = snow->theta_ice;
    double *theta_liq = snow->theta_liq;
    double *liq_deriv = snow->liq_deriv;
    double *snow_density = snow->density;
    double *snow_porosity = snow->porosity;
    double frac_ice[MAX_SOILS] = {0};
    double *matric = cell->matric;
    double *dWdQ_left = energy->dWdQ_left;
    double *dWdQ_right = energy->dWdQ_right;
    double *dWdE_left = energy->dWdE_left;
    double *dWdE_right = energy->dWdE_right;
    double *Wmac_node = soil_con->Wmac_node;
    double *Wsat_node = soil_con->Wsat_node;
    double *Wpwp_node = soil_con->Wpwp_node;
    double *Ksat_node = soil_con->Ksat_node;
    double *Kmac_node = soil_con->Kmac_node;
    double *dz_soil = soil_con->dz_soil;
    double *zc_soil = soil_con->zc_soil;
    double *Zsum_soil = soil_con->Zsum_soil;
    double *conductivity = cell->conductivity;
    double *conduct_int = cell->conduct_int;
    double *liquid_flux = cell->liquid_flux;
    double *soil_imped = cell->soil_imped;
    double *dtheta_dT = energy->dtheta_dT;
    double *dtheta_dpsi = cell->dtheta_dpsi;
    double *dconduct_dT = cell->dconduct_dT;
    double *dconduct_dpsi = cell->dconduct_dpsi;
    double *dpsiliq_dpsi = cell->dpsiliq_dpsi;

    for (i = 0; i < Nsnow; i++) {
        if (snow_density[i] >= 200.0) {
            liquid_capacity = param.SNOW_LIQUID_WATER_CAPACITY;
        }
        else {
            liquid_capacity = param.SNOW_LIQUID_WATER_CAPACITY + 0.07 * (200.0 - snow_density[i]) / 200.0;
        }
        double theta_resid = liquid_capacity * snow_porosity[i];
        double avail_cap = snow_porosity[i] - theta_resid;
        double snow_conduct = 0.0;
        if (theta_liq[i] > theta_resid && avail_cap > 1.0e-10) {
            double rel_sat = min(1.0, max(0.0, (theta_liq[i] - theta_resid) / avail_cap));
            liquid_flux[i] = snow_conduct * pow(rel_sat, param.SNOW_EXPMELT);
            liq_deriv[i] = snow_conduct * param.SNOW_EXPMELT / avail_cap * pow(rel_sat, param.SNOW_EXPMELT - 1.0);
        }
        else {
            liquid_flux[i] = 0.0;
            liq_deriv[i] = 0.0;
        }
    }

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
        dtheta_dpsi[i] = SoilWaterRetentionCurve(DERIVM_FLAG, i, liq[i], 
                                                 matric[i], soil_T[i], soil_con);

        double conduct_noice = SoilWaterRetentionCurve(CONDUCT_FLAG, i, liq[i], 
                                                       matric[i], soil_T[i], soil_con);
        // derivative of ice impedance factor with respect to temperature
        dimped_dT = log(10.0) * ice_param * pow(10.0, ice_param * frac_ice[i]) * dtheta_dT[i];
        // derivative for macropores
        if (liq[i] > Wmac_node[i]) {
            double relSatMP = (liq[i] - Wmac_node[i]) / (Wsat_node[i] - Wmac_node[i]);
            conduct_mp = (Kmac_node[i] - Ksat_node[i]) * pow(relSatMP, param.SOIL_EXPMAC);
            dmacro_dliq = ((Kmac_node[i] - Ksat_node[i]) / (Wsat_node[i] - Wmac_node[i])) * 
                            param.SOIL_EXPMAC * pow(relSatMP, param.SOIL_EXPMAC - 1.0);
            dmacro_dpsi = dmacro_dliq * dtheta_dpsi[i];
            dmacro_dT = dmacro_dliq * dtheta_dT[i];    
        }
        else {
            conduct_mp = 0.0;
            dmacro_dliq = 0.0;
            dmacro_dpsi = 0.0;
            dmacro_dT = 0.0;
        }
        // derivative for micropores
        if (ice[i] > 1.0e-12) {
            dTK_dpsi_noice = SoilWaterRetentionCurve(DERIVK_FLAG, i, liq[i], 
                                                     matric[i], soil_T[i], soil_con);
            dmicro_dT = dpsiliq_dT[i] * dTK_dpsi_noice;
            dmicro_dpsi = conduct_noice * dimped_dT * dtheta_dpsi[i] +
                                    dTK_dpsi_noice * soil_imped[i];
        }
        else {
            dmicro_dT = 0.0;
            dmicro_dpsi = SoilWaterRetentionCurve(DERIVK_FLAG, i, liq[i], 
                                                  matric[i], soil_T[i], soil_con);
        } 
        
        // combine matric derivatives
        dconduct_dpsi[i] = dmicro_dpsi + dmacro_dpsi;

        // combine temperature derivatives
        dconduct_dT[i] = conduct_noice * dimped_dT + dmicro_dT * soil_imped[i] + dmacro_dT;
        
        // 修正冻土导致的不透水性
        conductivity[i] = conduct_mp + soil_imped[i] * conduct_noice;
    }

    // Calculates the liquid moisture flux between soil nodes.
    for (i = 0; i < Nsoil - 1; i++) {
        // 设置临时变量
        double dzp = zc_soil[i+1] - zc_soil[i];
        double zdn = Zsum_soil[i] - zc_soil[i];          /* 界面到下节点距离 */
        double zdp = zc_soil[i+1] - Zsum_soil[i];        /* 界面上节点距离 */
        double denom = conductivity[i] * zdp + conductivity[i+1] * zdn;
        double k_int = conductivity[i] * conductivity[i+1] * dzp / denom; // 调和平均

        // 中间因子 dK_int/dK_L, dK_int/dK_R （调和平均）
        double dTKdK_left  = conductivity[i+1] * conductivity[i+1] * dzp * zdn / (denom * denom);
        double dTKdK_right = conductivity[i] * conductivity[i] * dzp * zdp / (denom * denom);

        // 界面导水率对左侧/右侧基质吸力的导数
        double dTKdpsi_left  = dTKdK_left  * dconduct_dpsi[i];
        double dTKdpsi_right = dTKdK_right * dconduct_dpsi[i+1];

        // 界面导水率对左侧/右侧温度的导数
        double dTKdT_left  = dTKdK_left  * dconduct_dT[i];
        double dTKdT_right = dTKdK_right * dconduct_dT[i+1];

        // 界面水通量
        conduct_int[i] = k_int / dzp;
        liquid_flux[i] = conduct_int[i] * (matric[i] - matric[i+1] + dzp);

        double delt_psi = matric[i+1] - matric[i];

        // 对左侧（上层）基质吸力的导数
        dWdQ_left[i] = -dTKdpsi_left * delt_psi / dzp + conduct_int[i] + dTKdpsi_left * dpsiliq_dpsi[i];

        // 对右侧（下层）基质吸力的导数
        dWdQ_right[i] = -dTKdpsi_right * delt_psi / dzp - conduct_int[i] + dTKdpsi_right * dpsiliq_dpsi[i+1];

        // 对左侧（上层）温度的导数
        dWdE_left[i] = -dTKdT_left * delt_psi / dzp + conduct_int[i] * dpsiliq_dT[i] + dTKdT_left;

        // 对右侧（下层）温度的导数
        dWdE_right[i] = -dTKdT_right * delt_psi / dzp - conduct_int[i] * dpsiliq_dT[i+1] + dTKdT_right;
    }
    // 
    if (options.NOFLUX) {
        conduct_int[Nsoil-1] = 0.0;
        liquid_flux[Nsoil-1] = 0.0;
        dWdQ_left[Nsoil-1] = 0.0;
        dWdQ_right[Nsoil-1] = 0.0;
        dWdE_left[Nsoil-1] = 0.0;
        dWdE_right[Nsoil-1] = 0.0;
    }
}