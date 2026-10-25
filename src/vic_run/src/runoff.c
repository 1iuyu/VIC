/******************************************************************************
* @section DESCRIPTION
*
* Calculate infiltration and runoff from the surface, gravity driven drainage
* between all soil layers, and generates baseflow from the bottom layer.
******************************************************************************/

#include "vic_run.h"

/******************************************************************************
* @brief    Calculate infiltration and runoff from the surface, gravity driven
*           drainage between all soil layers, and generates baseflow from the
*           bottom layer.
******************************************************************************/
int
runoff(double             step_dt,
       double             inflow,
       energy_bal_struct *energy,
       cell_data_struct  *cell,
       soil_con_struct   *soil_con)
{
    size_t i;
    double ex = 0.0;
    double top_moist = 0.0;
    double top_max_moist = 0.0;
    double infil_area, sat_area;
    // 指针赋值 
    double *liq = energy->liq;
    double *ice = energy->ice;
    double *moist = energy->moist;
    double *soil_T = cell->soil_T;
    double *dz_soil = soil_con->dz_soil;
    double *matric = cell->matric;
    double *dtheta_dT = energy->dtheta_dT;
    double *dtheta_dpsi = cell->dtheta_dpsi;
    double *Zsum_soil = soil_con->Zsum_soil;
    double *excess_ice = cell->excess_ice;
    double *dWdQ_surf = energy->dWdQ_surf;
    double *dWdE_surf = energy->dWdE_surf;
    double *Ksat_node = soil_con->Ksat_node;
    double *Wsat_node = soil_con->Wsat_node;
    double *eff_porosity = cell->eff_porosity;  // 有效孔隙度
    // 初始化
    size_t Nsoil = cell->Nsoil;
    size_t Nroot = cell->Nroot;
    double root_depth = cell->root_depth;
    double b_infilt = soil_con->b_infilt;
    double capil_drive = soil_con->capil_drive;
    cell->soil_excess = 0.0; // 土壤超出水分[m]
    inflow *= step_dt;
    /*************************
       Initialize Variables
    *************************/
    double baseflow = 0.0;
    // 计算土壤超出水分
    for (i = 0; i < Nsoil; i++) {
        eff_porosity[i] = max(0.0, Wsat_node[i] - ice[i]);
        /** Set Layer Liquid Moisture Content **/
        if (liq[i] > eff_porosity[i]) {
            cell->soil_excess += (liq[i] - eff_porosity[i]) * dz_soil[i];
            liq[i] = eff_porosity[i];
        }
        moist[i] = liq[i] + (ice[i] + excess_ice[i]) *
                         CONST_RHOICE / CONST_RHOFW;
        // 土壤水分诊断
        if (moist[i] > Wsat_node[i] + 1.0e-12) {
            log_err("Soil moisture exceeds saturation: layer = %zu, moist = %.4f, Wsat = %.4f, ",
                    i, moist[i], Wsat_node[i]);
        }
    }
    // 计算土层水量
    for (i = 0; i < Nroot; i++) {
        top_moist += liq[i] * dz_soil[i];
        top_max_moist += Wsat_node[i] * dz_soil[i];
    }
    top_moist = min(top_moist, top_max_moist);
    double base = max(0.0, 1.0 - (top_moist / top_max_moist));
    if (top_moist <= 0.0) {
        infil_area = 1.0;
        sat_area = 0.0;
    }
    else {
        ex = 1.0 / (1.0 + b_infilt);
        infil_area = pow(base, ex);
        sat_area = 1.0 - infil_area;
    }
    double i_m = (1.0 + b_infilt) * top_max_moist;
    double i_0 = i_m * sat_area;
    // 计算根区容量
    double rootliq = 0.0, rootice = 0.0, rootCapacity = 0.0;
    for (i = 0; i < Nroot; i++) {
        rootliq += dz_soil[i] * liq[i];
        rootice += dz_soil[i] * ice[i];
        rootCapacity += dz_soil[i] * Wsat_node[i];
    }
    double compHeadRootZone = 0.0;
    double compHeadCutoff = 1.0;
    double compHeadWidth = 0.02;
    double headSmooth = 1.0e-4;
    double infil_area_prev = infil_area;
    for (i = 0; i < Nroot; i++) {
        double denom = sqrt(matric[i] * matric[i] + headSmooth * headSmooth);
        double posHead = 0.5 * (matric[i] + denom);
        compHeadRootZone += posHead * dz_soil[i];
    }
    compHeadRootZone /= Zsum_soil[Nsoil-1];
    double compArg = (compHeadRootZone - compHeadCutoff) / compHeadWidth;
    double compLimiter = 1.0 / (1.0 + exp(2.0 * compArg));
    double dCompLimiter_dHead = -(2.0 / compHeadWidth) * compLimiter * (1.0 - compLimiter);
    infil_area = infil_area * compLimiter;
    double availCapacity = rootCapacity - rootice;
    double wetfront = rootliq / availCapacity * root_depth;
    double f_0 = Ksat_node[0] * (1.0 + (1.0 - wetfront / Zsum_soil[Nsoil-1]) * capil_drive / wetfront);
    double f_m = f_0 * (soil_con->b_dynamic + 1.0);
    double dliq_dWat, dice_dWat;
    for (i = 0; i < Nroot; i++) {
        double soil_tkfrz = calc_tkfrz(i, NULL, energy, cell, soil_con);
        double dInfilArea_dS1 = -ex * pow(base, ex - 1.0) / top_max_moist;
        if (soil_T[i] < soil_tkfrz) {
            dliq_dWat = 0.0;
            dice_dWat = dtheta_dpsi[i];
        }
        else {
            dice_dWat = 0.0;
            dliq_dWat = dtheta_dpsi[i];
        }
        // 根区导数
        double drootliq_dWat = dliq_dWat * dz_soil[i];
        double drootliq_dTk = dtheta_dT[i] * dz_soil[i];
        double drootice_dWat = dice_dWat * dz_soil[i];
        double drootice_dTk = -dtheta_dT[i] * dz_soil[i];
        // 湿润锋深度导数
        double dwetfront_dWat = (drootliq_dWat * root_depth +
                                 drootice_dWat * wetfront) / availCapacity;
        double dwetfront_dTk = (drootliq_dTk * root_depth +
                                drootice_dTk * wetfront) / availCapacity;
        // 最大下渗率（Green-Ampt）的导数
        double dinfilrate_dWat = -Ksat_node[i] * capil_drive * dwetfront_dWat / (wetfront * wetfront);
        double dinfilrate_dTk = -Ksat_node[i] * capil_drive * dwetfront_dTk / (wetfront * wetfront);
        // 下渗面积导数
        double dInfilArea_dWat;
        double dInfilArea_dTk = dInfilArea_dS1 * dz_soil[i] * dtheta_dT[i];
        if (soil_T[i] < soil_tkfrz) {
            dInfilArea_dWat = 0.0;
        } else {
            dInfilArea_dWat = dInfilArea_dS1 * dz_soil[i] * dtheta_dpsi[i];
        }
        
        // 压缩限制因子修正
        double denom = sqrt(matric[i] * matric[i] + headSmooth * headSmooth);
        double dPosHead_dPsi = 0.5 * (1.0 + matric[i] / max(denom, 1.0e-12));
        double dCompHead_dWat = (dz_soil[i] / root_depth) * dPosHead_dPsi;
        dInfilArea_dWat = dInfilArea_dWat * compLimiter +
                          infil_area_prev * dCompLimiter_dHead * dCompHead_dWat;
        dInfilArea_dTk  = dInfilArea_dTk * compLimiter;
        // 地表下渗对于水量和能量的导数
        dWdQ_surf[i] = dInfilArea_dWat * f_0 + infil_area * dinfilrate_dWat;
        dWdE_surf[i] = dInfilArea_dTk * f_0 + infil_area * dinfilrate_dTk;          
    }

    /*************************************************
       Compute surface runoff
    *************************************************/
    // Compute infiltration rate, saturation excess runoff and infiltration excess runoff.
    calc_dynamicVIC(inflow, i_0, i_m,
                    f_0, f_m, step_dt,
                    cell, soil_con);

    /************************************************
      Compute Baseflow and watertable
    ************************************************/
    baseflow = wrap_compute_zwt(step_dt, energy,
                                cell, soil_con);

    /** Store tile-wide values **/    
    cell->baseflow = baseflow; // [m/s]

    return (0);
}

/******************************************************************************
* @brief    Calculate the saturated area, surface infiltration, saturation 
*           excess runoff and infiltration excess runoff based on moist.
******************************************************************************/
int
calc_dynamicVIC(double            inflow, // m
                double            i_0,
                double            i_m,
                double            f_0,
                double            f_m,
                double            step_dt,   // 时间步长[s]
                cell_data_struct *cell,
                soil_con_struct  *soil_con)
{
    extern parameters_struct param;
    bool converged = false;
    double runoff_sat = 0.0;            // saturation excess runoff [m]
    double runoff_infil = 0.0;          // infiltration excess runoff [m]
    double dsat_dy = 0.0;
    double dinfil_dy = 0.0;
    double b_infilt = soil_con->b_infilt;
    double b_dynamic = soil_con->b_dynamic;
    
    if (inflow <= 0.0) {
        cell->runoff = 0.0;
        return (0); 
    }

    size_t iter = 0;
    double y = 0.5 * inflow; // 初始猜测y

    do {
        // 计算蓄满产流和导数
        calc_sat_runoff(i_0, i_m, y, 
                        b_infilt, 
                        &runoff_sat, 
                        &dsat_dy);

        // 计算超渗产流和导数
        calc_infil_runoff(step_dt, inflow, 
                        f_m, f_0, 
                        dsat_dy, b_dynamic, 
                        runoff_sat, 
                        &dinfil_dy, 
                        &runoff_infil);

        double F_val = y + runoff_infil - inflow;
        double dF_val = 1.0 + dinfil_dy;
        // 检查导数是否接近零（避免除零）
        if (fabs(F_val) < param.TOL_A) {
            converged = true;
            break;
        }
        if (fabs(dF_val) < 1.0e-12) {
            break;
        }
        double y_new = y - F_val / dF_val;

        if (y_new <= 0.0 || y_new >= inflow) {
            break;
        }
        // 检查收敛条件
        if (fabs(y_new - y) < param.TOL_A) {
            y = y_new;
            calc_sat_runoff(i_0, i_m, y, 
                            b_infilt, 
                            &runoff_sat, 
                            &dsat_dy);
            calc_infil_runoff(step_dt, inflow, 
                              f_m, f_0, dsat_dy, 
                              b_dynamic, 
                              runoff_sat, 
                              &dinfil_dy, 
                              &runoff_infil);
            if (fabs(y + runoff_infil - inflow) < param.TOL_A) {
                converged = true;
            }
            break;
        }       

        // 更新迭代数值
        y = y_new;

        iter++;
    }
    while (iter < param.MAX_ITER_MOST);

    // 备选二分法
    if (!converged) {
        double y_low = 0.0;
        double y_high = inflow;
        double F_low;
        double F_high;

        calc_sat_runoff(i_0, i_m, y_low, b_infilt, &runoff_sat, &dsat_dy);
        calc_infil_runoff(step_dt, inflow, f_m, f_0, dsat_dy, b_dynamic, runoff_sat, &dinfil_dy, &runoff_infil);
        F_low = y_low + runoff_infil - inflow;

        calc_sat_runoff(i_0, i_m, y_high, b_infilt, &runoff_sat, &dsat_dy);
        calc_infil_runoff(step_dt, inflow, f_m, f_0, dsat_dy, b_dynamic, runoff_sat, &dinfil_dy, &runoff_infil);
        F_high = y_high + runoff_infil - inflow;

        if (F_low > 0.0 || F_high < 0.0) {
            log_err("Cannot bracket root in calc_dynamicVIC: "
                    "F_low = %.6f, F_high = %.6f, inflow = %.6f", F_low, F_high, inflow);
            return (1);
        }

        for (iter = 0; iter < param.MAX_ITER_MOST; iter++) {
            y = 0.5 * (y_low + y_high);

            calc_sat_runoff(i_0, i_m, y, b_infilt, &runoff_sat, &dsat_dy);
            calc_infil_runoff(step_dt, inflow, f_m, f_0, dsat_dy, b_dynamic, runoff_sat, &dinfil_dy, &runoff_infil);

            double F_val = y + runoff_infil - inflow;

            if (fabs(F_val) < param.TOL_A || 0.5 * (y_high - y_low) < param.TOL_A) {
                converged = true;
                break;
            }

            if (F_val > 0.0) {
                y_high = y;
            }
            else {
                y_low = y;
            }
        }
    }

    if (!converged) {
        log_err("Newton iteration and bisection failed in calc_dynamicVIC: inflow = %.6f, y = %.6f", inflow, y);
        return (ERROR);
    }

    /* Recalculate runoff using the final y */
    calc_sat_runoff(i_0, i_m, y, b_infilt, &runoff_sat, &dsat_dy);
    calc_infil_runoff(step_dt, inflow, f_m, f_0, dsat_dy, b_dynamic, runoff_sat, &dinfil_dy, &runoff_infil);

    // 更新地表产流和入渗量
    cell->runoff_infil = runoff_infil / step_dt; // [m/s]
    cell->runoff_sat = runoff_sat / step_dt;     // [m/s]
    cell->infiltration = (y - runoff_sat) / step_dt;         // [m/s]
    cell->runoff = (runoff_sat + runoff_infil) / step_dt;    // [m/s]

    return (0);

}

/******************************************************************************
* @brief    Calculate the saturated area and runoff
******************************************************************************/
void
calc_sat_runoff(double    i_0,
                double    i_m,
                double    y,
                double    b_infilt,
                double   *runoff_sat,
                double   *dsat_dy)  // 蓄满产流对y的导数
{
    double  zwt;
    zwt = i_0 + y;
    if (zwt > i_m) {
        zwt = i_m;
    }
    // 饱和面积未覆盖全区域
    if (y <= i_m - i_0) {
        (*runoff_sat) = y - (i_m / (b_infilt + 1.0)) * 
                        ((pow(1.0 - (i_0 / i_m), b_infilt + 1.0)) - 
                        (pow(1.0 - (zwt / i_m), b_infilt + 1.0)));
        (*dsat_dy) = 1.0 - pow(1.0 - (i_0 + y) / i_m, b_infilt);
    }
    // 饱和面积覆盖全区域
    else {
        (*runoff_sat) = i_m - i_0 - (i_m / (b_infilt + 1.0)) * 
                        pow(1.0 - i_0 / i_m, b_infilt + 1.0) + y - (i_m - i_0);
        (*dsat_dy) = 1.0;
    }

    if (*runoff_sat < 0.0) {
        *runoff_sat = 0.0;
    }
}

/******************************************************************************
* @brief    Calculate the saturated area and runoff
******************************************************************************/
void
calc_infil_runoff(double    step_dt,
                  double    inflow,
                  double    f_m,
                  double    f_0,
                  double    dsat_dy,
                  double    b_dynamic,
                  double    runoff_sat,
                  double   *dinfil_dy,
                  double   *runoff_infil)
{
    double inflow_eff = max(0.0, inflow - runoff_sat);
    double ratio = inflow_eff / (f_m * step_dt);
    if (inflow_eff <= 0.0) {
        *runoff_infil = 0.0;
        *dinfil_dy = 0.0;
    }
    // 有效降雨强度 ≤ 最大入渗率
    else if (ratio <= 1.0) {
        (*runoff_infil) = inflow_eff - (f_0 * step_dt) * 
                        (1.0 - pow(1.0 - ratio, b_dynamic + 1.0));
        (*dinfil_dy) = -dsat_dy * (1.0 - pow(1.0 - ratio, b_dynamic));
    }
    // 有效降雨强度 ≥ 最大入渗率
    else {
        (*runoff_infil) = inflow_eff - (f_0 * step_dt);
        (*dinfil_dy) = -dsat_dy;
    }

    // 数值保护
    if (*runoff_infil < 0.0) {
        *runoff_infil = 0.0;
    }
    if (*runoff_infil > inflow) {
        *runoff_infil = inflow;
    }
}