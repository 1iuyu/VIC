/******************************************************************************
 * @section DESCRIPTION
 *
 * Compute spatial average water table position (zwt) and discharge.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    Compute spatial average water table position (zwt).  Water table
 *           position is measured in m and is positive below the soil surface.
 *****************************************************************************/
double
wrap_compute_zwt(double            step_dt,
                 cell_data_struct *cell,
                 soil_con_struct  *soil_con)
{
    extern parameters_struct param;
    size_t  i, j, lidx;
    size_t  Nfrost;
    double  decay_fator;
    double  discharge;
    double  baseflow = 0.0;
    double  recharge;
    double  recharge_layer;
    double  aqf_yield;  // aquifer yield [m]
    double  aqf_yield1; // 含水层给水度 [m]
    double  soil_excess;
    double  frac_ice[MAX_SOILS];
    double *moist = cell->moist;
    double *ice = cell->ice;
    double *liq = cell->liq;
    double *soil_T = cell->soil_T;
    double *matric = cell ->matric;
    double *zc_soil = soil_con->zc_soil;
    double *dz_soil = soil_con->dz_soil;
    double *Zsum_soil = soil_con->Zsum_soil;
    double *expt_node = soil_con->expt_node;
    double *soil_imped = cell->soil_imped;
    double *ksat_node = soil_con->Ksat_node;
    double *Wsat_node = soil_con->Wsat_node;
    double *Wpwp_node = soil_con->Wpwp_node;
    double *porosity = cell->porosity;
    double *alpha_node = soil_con->alpha_node;
    double *mpar_node = soil_con->mpar_node;
    double *conduct_int = cell->conduct_int;

    size_t Nsoil = cell->Nsoil;
    double zwt = cell->zwt; // 地下水位
    double storage_aqf = cell->storage_aqf;

    for (i = 0; i < Nsoil; i++) {
        frac_ice[i] = max(0.01, min(1.0, ice[i] / Wsat_node[i]));
    }
    // 寻找第一个非饱和层
    size_t zwt_lidx = Nsoil;
    for (i = 0; i < Nsoil; i++) {
        if (zwt <= Zsum_soil[i]) {
            zwt_lidx = i;
            break;
        }
    }
    double waterhead = 0.0;
    double aqf_conduct = 0.0;
    if (zwt_lidx < Nsoil) {
        // 非饱和带最底层的水头
        waterhead = matric[zwt_lidx] + zc_soil[zwt_lidx];
        aqf_conduct = soil_imped[zwt_lidx] * conduct_int[zwt_lidx];
        recharge = aqf_conduct * (waterhead - zwt) / (zwt - zc_soil[zwt_lidx]);
        double max_flux = 0.01 / step_dt;  // 最大允许通量
        recharge = max(-max_flux, min(max_flux, recharge)); // m/s
    }
    // 保存下泄量
    cell->recharge = recharge;
    lidx = Nsoil - 1; // 土层下标最后一层
    aqf_yield = (Wsat_node[lidx] - Wpwp_node[lidx]) * (1.0 - pow(1.0 + 
                            pow(alpha_node[lidx] * zwt, expt_node[lidx]), -mpar_node[lidx]));
    aqf_yield = max(0.02, aqf_yield);
    double recharge_flux = recharge * step_dt;
    // 地下水位位于土层下方
    if (zwt_lidx == Nsoil) {
        storage_aqf += recharge_flux;
        zwt -= recharge_flux / aqf_yield;
    }
    else {
        if (recharge_flux > 0.0) { // 地下水位上升
            for (i = zwt_lidx + 1; i > 0; i--) { // 防止i下溢到SIZE_MAX
                j = i - 1; // j是当前层的下标
                aqf_yield1 = (Wsat_node[j] - Wpwp_node[j]) * (1.0 - pow(1.0 + 
                                        pow(alpha_node[j] * zwt, expt_node[j]), -mpar_node[j]));
                aqf_yield1 = max(0.02, aqf_yield1);
                double layer_top;
                if (j > 0) {
                    layer_top = Zsum_soil[j-1];  // 第j-1层的顶部
                } else {
                    layer_top = 0.0;  // 地表
                }
                recharge_layer = min(recharge_flux, (aqf_yield1 * (zwt - layer_top)));
                recharge_layer = max(0.0, recharge_layer);
                if (aqf_yield1 > 0.0) {
                    zwt -= recharge_layer / aqf_yield1;
                }
                recharge_flux -= recharge_layer;
                if (recharge_flux <= 0.0) {
                    break;
                }
            }        
        }
        else {  // 地下水位下降
            for (i = zwt_lidx; i < Nsoil; i++) {
                aqf_yield1 = (Wsat_node[i] - Wpwp_node[i]) * (1.0 - pow(1.0 + 
                                        pow(alpha_node[i] * zwt, expt_node[i]), -mpar_node[i]));
                aqf_yield1 = max(0.02, aqf_yield1);
                recharge_layer = max(recharge_flux, -(aqf_yield1 * (Zsum_soil[i] - zwt)));
                recharge_layer = min(0.0, recharge_layer);
                recharge_flux -= recharge_layer;
                if (recharge_flux >= 0.0) {
                    zwt -= recharge_layer / aqf_yield1;
                    break;
                }
                else {
                    zwt = Zsum_soil[i];
                }
            }
            if (recharge_flux > 0.0) {
                zwt -= recharge_flux / aqf_yield;
            }
        }
        // 重新计算地下水位
        zwt_lidx = Nsoil;
        for (i = 0; i < Nsoil; i++) {
            if (zwt <= Zsum_soil[i]) {
                zwt_lidx = i;
                break;
            }
        }
    }
    // 计算BASEFLOW
    if (cell->soil_T[0] > CONST_TKFRZ) {
        Nfrost = Nsoil - 1;
    }
    else {
        Nfrost = 0;
    }
    for (i = 1; i < Nsoil; i++) {
        if (soil_T[i-1] > CONST_TKFRZ && soil_T[i] <= CONST_TKFRZ) {
            Nfrost = i;
            break;
        }
    }
    double frost_zwt = zc_soil[Nfrost];
    double zwt_perch = frost_zwt;
    double aqf_perch = 0.0;
    double q_perch = 0.0;
    double tmp_imped = 0.0;
    double drain_perch = 0.0;
    double mass_liq = 0.0;
    double max_drain = 1.0e-05 * sin(soil_con->slope * (CONST_PI / 180.0));
    // water table above frost table
    if (zwt < frost_zwt && soil_T[Nfrost] <= CONST_TKFRZ) {
        for (i = zwt_lidx; i < Nfrost; i++) {
            lidx = min(Nsoil-1, i+1);
            tmp_imped = pow(10.0, 0.5 * (frac_ice[i] + frac_ice[lidx]) * -6.0);
            q_perch += ksat_node[i] * tmp_imped * dz_soil[i];
            aqf_perch += dz_soil[i];
        }
        if (aqf_perch > 0.0) {
            q_perch /= aqf_perch;
        }
        drain_perch = max_drain * q_perch * (frost_zwt - zwt);
        double sub_drain = -drain_perch * step_dt;
        double drain_layer = 0.0;
        for (i = zwt_lidx; i < Nfrost; i++) {
            mass_liq = liq[i] * dz_soil[i];
            drain_layer = max(sub_drain, -(mass_liq - MIN_SOILMOIST));
            drain_layer = min(0.0, drain_layer);
            sub_drain -= drain_layer;
            liq[i] += drain_layer / dz_soil[i];
            if (sub_drain >= 0.0) {
                zwt -= drain_layer / porosity[i];
                break;
            }
            else {
                zwt = Zsum_soil[i];
            }
        }
        drain_perch += sub_drain / step_dt;
        // 重新计算地下水位
        zwt_lidx = Nsoil;
        for (i = 0; i < Nsoil; i++) {
            if (zwt <= Zsum_soil[i]) {
                zwt_lidx = i;
                break;
            }
        }
    }
    else {   // water table below frost table
        size_t k_perch = 0;
        double sub_drain = 0.0;
        double drain_layer = 0.0;
        for (i = Nfrost + 1; i > 0; i--) { // 防止i下溢到SIZE_MAX
            j = i - 1; // j是当前层的下标
            if (moist[j] / Wsat_node[j] <= 0.9) {
                k_perch = j;
            }
        }
        if (soil_T[Nfrost] > CONST_TKFRZ) {
            k_perch = Nfrost;
        }
        if (Nfrost > k_perch && k_perch < Nsoil - 1) {
            double frac1 = (liq[k_perch] + ice[k_perch]) / Wsat_node[k_perch+1];
            double frac2 = (liq[k_perch+1] + ice[k_perch+1]) / Wsat_node[k_perch+1];
            if (fabs(frac2 - frac1) < param.TOL_A) {
                zwt_perch = zc_soil[k_perch];
            }
            else {
                double m = (zc_soil[k_perch+1] - zc_soil[k_perch]) / (frac2 - frac1);
                double b = zc_soil[k_perch+1] - m * frac2;
                zwt_perch = max(0.9 * m + b, 0.0);
            }
            for (i = k_perch; i < Nfrost; i++) {
                lidx = min(Nsoil-1, i+1);
                tmp_imped = pow(10.0, 0.5 * (frac_ice[i] + frac_ice[lidx]) * -6.0);
                q_perch += ksat_node[i] * tmp_imped * dz_soil[i];
                aqf_perch += dz_soil[i];
            }
            if (aqf_perch > 0.0) {
                q_perch /= aqf_perch;
            }
            drain_perch = max_drain * q_perch * (frost_zwt - zwt_perch);
            sub_drain = -drain_perch * step_dt;
            for (i = k_perch; i < Nfrost; i++) {
                mass_liq = liq[i] * dz_soil[i];
                drain_layer = max(sub_drain, -(mass_liq - MIN_SOILMOIST));
                drain_layer = min(0.0, drain_layer);
                sub_drain -= drain_layer;
                liq[i] += drain_layer / dz_soil[i];
                if (sub_drain > 0.0) {
                    zwt_perch -= drain_layer / porosity[i];
                    break;
                }
                else {
                    zwt_perch = Zsum_soil[i];
                }
            }
            drain_perch += sub_drain / step_dt;
        }
        else {
            drain_perch = 0.0;
        }
        // Topographic runoff
        double total_ice = 0.0;
        double tmp_liq = 0.0;
        double dz_sum = 0.0;
        if (zwt_lidx < Nsoil) {
            for (i = zwt_lidx; i < Nsoil; i++) {
                dz_sum += dz_soil[i];
                total_ice += frac_ice[i] * dz_soil[i];
            }
            tmp_imped = pow(10.0, (total_ice / dz_sum * -6.0));
        }
        else {
            tmp_imped = pow(10.0, -6.0 * frac_ice[Nsoil-1]);
        }
        max_drain = min(0.01 * sin((CONST_PI / 180.0) * soil_con->slope), 0.01);
        if (zwt_lidx == Nsoil) {
            decay_fator = expt_node[lidx] / 3.0;
        }
        else {
            decay_fator = expt_node[zwt_lidx] / 3.0;
        }
        // coef_baseflow = conductivity[layer] * exp(3.0);  // ???
        baseflow = tmp_imped * max_drain * exp(-decay_fator * zwt);
        // baseflow = soil_imped * coef_baseflow * exp(-decay_fator * zwt);
        aqf_yield = (Wsat_node[lidx] - Wpwp_node[lidx]) * (1.0 - pow(1.0 + 
                                pow(alpha_node[lidx] * zwt, expt_node[lidx]), -mpar_node[lidx]));
        aqf_yield = max(0.02, aqf_yield);
        double baseflow_flux = baseflow * step_dt;
        // 地下水位在土层下方
        if (zwt_lidx == Nsoil) {
            storage_aqf -= baseflow_flux;
            zwt += baseflow_flux / aqf_yield;
            tmp_liq = max(0.0, storage_aqf - 5.0);
            liq[lidx] += tmp_liq / dz_soil[lidx];
            storage_aqf = min(5.0, storage_aqf);
        }
        else { // 地下水位在土层内
            sub_drain = -baseflow_flux;
            if (sub_drain > 0.0) {
                log_err("Error: baseflow exceeds soil moisture capacity.");
            }
            else {
                for (i = zwt_lidx; i < Nsoil; i++) {
                    aqf_yield1 = (Wsat_node[i] - Wpwp_node[i]) * (1.0 - pow(1.0 + 
                                        pow(alpha_node[i] * zwt, expt_node[i]), -mpar_node[i]));
                    aqf_yield1 = max(0.02, aqf_yield1);
                    drain_layer = max(sub_drain, -(aqf_yield1 * (Zsum_soil[i] - zwt)));
                    drain_layer = min(0.0, drain_layer);
                    sub_drain -= drain_layer;
                    liq[i] += drain_layer / dz_soil[i];
                    if (sub_drain >= 0.0) {
                        zwt -= drain_layer / aqf_yield1;
                        break;
                    }
                    else {
                        zwt = Zsum_soil[i];
                    }
                }
                zwt -= sub_drain / aqf_yield;
                storage_aqf += sub_drain;
            }
            // 重新计算地下水位
            zwt_lidx = Nsoil;
            for (i = 0; i < Nsoil; i++) {
                if (zwt <= Zsum_soil[i]) {
                    zwt_lidx = i;
                    break;
                }
            }
        }
        zwt = max(0.0, zwt);
        zwt = min(80.0, zwt);
    }
    cell->zwt = zwt;
    cell->storage_aqf = storage_aqf;
    // 检查是否存在过度饱和的层
    for (i = Nsoil - 1; i >= 1; i--) {
        porosity[i] = max(1.0e-4, Wsat_node[i] - ice[i]);
        soil_excess = max(0.0, liq[i] - porosity[i]) * dz_soil[i];
        liq[i] = min(liq[i], porosity[i]);
        liq[i-1] += soil_excess / dz_soil[i-1];
    }
    // 对于第一层，额外检查是否小于MIN_SOILMOIST
    porosity[0] = max(1.0e-4, Wsat_node[0] - ice[0] - MIN_SOILMOIST);
    soil_excess = max(max(0.0, liq[0] - MIN_SOILMOIST) - porosity[0], 0.0) * dz_soil[0];
    liq[0] -= soil_excess / dz_soil[0];
    double sub_excess = soil_excess / step_dt;
    // 额外检查是否小于MIN_SOILMOIST
    double soil_loss = 0.0;
    for (i = 0; i <= Nsoil - 2; i++) {
        if (liq[i] < MIN_SOILMOIST) {
            soil_loss = (MIN_SOILMOIST - liq[i]) * dz_soil[i];
            if (i == zwt_lidx) {
                zwt += soil_loss / porosity[i];
            }
        }
        else {
            soil_loss = 0.0;
        }
        liq[i] += soil_loss / dz_soil[i];
        liq[i+1] -= soil_loss / dz_soil[i+1];
    }
    lidx = Nsoil - 1;
    if (liq[lidx] < MIN_SOILMOIST) {
        soil_loss = (MIN_SOILMOIST - liq[lidx]) * dz_soil[lidx];
        for (i = Nsoil; i > 0; i--) { // 防止i下溢到SIZE_MAX
            j = i - 1; // j是当前层的下标
            mass_liq = max(0.0, (liq[j] - MIN_SOILMOIST) * dz_soil[j]);
            if (mass_liq >= soil_loss) {
                liq[lidx] += soil_loss / dz_soil[lidx];
                liq[j] -= soil_loss / dz_soil[j];
                soil_loss = 0.0;
                break;
            }
            else {
                liq[lidx] += mass_liq / dz_soil[lidx];
                liq[j] -= mass_liq / dz_soil[j];
                soil_loss -= mass_liq;
            }
        }
        // 如果还不够，从地下水补充
        if (soil_loss > 0.0) {
            double draw = min(soil_loss / step_dt, baseflow);
            liq[lidx] += draw * step_dt / dz_soil[lidx];
            baseflow -= draw;
        }
    }
    for (i = 0; i < Nsoil; i++) {
        moist[i] = liq[i] + ice[i];
    }
    // 计算排泄量
    discharge = baseflow + sub_excess;
    cell->soil_excess = soil_excess;

    return(discharge);
}

