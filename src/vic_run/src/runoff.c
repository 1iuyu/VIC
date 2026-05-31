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
       cell_data_struct  *cell,
       soil_con_struct   *soil_con)
{
    size_t i, Nsoil;
    double infil_rate = 0.0;
    // 指针赋值 
    double *liq = cell->liq;
    double *ice = cell->ice;
    double *dz_soil = soil_con->dz_soil;
    double *Wsat_node = soil_con->Wsat_node;
    double *porosity = cell->porosity;  // 有效孔隙度
    // 初始化
    Nsoil = cell->Nsoil;
    /*************************
       Initialize Variables
    *************************/
    double infil_acc = 1.0e-06;
    double runoff = 0.0;
    double baseflow = 0.0;
    // 计算土壤超出水分
    for (i = 0; i < Nsoil; i++) {
        porosity[i] = max(1.0e-4, Wsat_node[i] - ice[i]);
        /** Set Layer Liquid Moisture Content **/
        cell->soil_excess += max(0.0, liq[i] - porosity[i]) * dz_soil[i];
        if (liq[i] > porosity[i]) {
            liq[i] = porosity[i];
        }
    }
    /**************************************************
       Compute Flow Between Soil Layers
    **************************************************/
    // 计算入渗率、蓄满产流和超渗产流
    calc_dynamicVIC(inflow, &infil_acc,
                    step_dt, &runoff,
                    &infil_rate, 
                    cell, soil_con);

    /************************************************
      Compute Baseflow and watertable
    ************************************************/
    baseflow = wrap_compute_zwt(step_dt,
                                cell, soil_con);

    /** Store tile-wide values **/    
    cell->runoff = runoff;     // [m/s]
    cell->baseflow = baseflow; // [m/s]

    return (0);
}

/******************************************************************************
* @brief    Calculate the saturated area, surface infiltration, saturation 
*           excess runoff and infiltration excess runoff based on moist.
******************************************************************************/
void
calc_dynamicVIC(double            dt_inflow, // 平均地表入流[m/s]
                double           *infil_acc, // 累计入渗率[m/s]
                size_t            iter_dt,   // 时间步长[s]
                double           *runoff,    // 地表径流[m/s]
                double           *infil_rate,
                cell_data_struct *cell,
                soil_con_struct  *soil_con)
{
    size_t  lidx;
    size_t  iter;
    bool    begin_flag = false;
    bool    middle_flag = false;
    bool    finally_flag = false;
    double  A, ex, i_0;  // i_0临时入渗率[m/s]
    double  max_infil;   // maximum infiltration rate [m/s]
    double  tmp_var1;              // temporary variable
    double  runoff_sat = 0.;            // saturation excess runoff [m/s]
    double  runoff_infil = 0.;          // infiltration excess runoff [m/s]
    double  *moist = cell->moist;
    double  *dz_soil = soil_con->dz_soil;
    double  *Wsat_node = soil_con->Wsat_node;
    double  last_depth;
    double  infil_capacity;        // infiltration capacity, [m]
    double  max_infil_capacity;    // maximum infiltration capacity [m]

    // initialization
    double top_moist = 0.;
    double top_max_moist = 0.;
    double init_depth = 0.;
    double tmp_sat = 0.;
    double tmp_sat1 = 0.;
    double tmp_infil = 0.;
    size_t last_iter = 20;
    double tmp_depth = 0.;
    double tol_error = 1.388889e-07 * iter_dt; // 0.5 mm per hour time step
    double capil_drive = soil_con->capil_drive;
    double b_infilt = soil_con->b_infilt;
    double b_dynamic = soil_con->b_dynamic;
    double *conductivity = cell->conductivity;

    for (lidx = 0; lidx < 3; lidx++) {
        top_moist += moist[lidx] * dz_soil[lidx];
        top_max_moist += Wsat_node[lidx] * dz_soil[lidx];
    }
    if (top_moist > top_max_moist) {
        top_moist = top_max_moist;
    }
    // 进入土层表面的水量[m]
    double inflow = dt_inflow * iter_dt;
    
    /** A as in Wood et al. in JGR 97, D3, 1992 equation (1) **/
    ex = 1.0 / (1.0 + b_infilt);
    A = 1.0 - (1.0 - pow((top_moist / top_max_moist), ex));
    
    max_infil_capacity = (1.0 + b_infilt) * top_max_moist;
    infil_capacity = max_infil_capacity * A;

    /** compute surface infiltration **/
    calc_soil_infil(moist[0], Wsat_node[0],
                    capil_drive,
                    dz_soil[0], conductivity[0],
                    soil_con->Ksat_node[0],
                    dt_inflow, &i_0, infil_acc);

    /** equation (3a) Wood et al. **/
    max_infil = i_0 * (b_dynamic + 1.0);

    if (inflow <= 0.0) {
        runoff_sat = 0.0;
        runoff_infil = 0.0;
        tmp_infil = 0.0;
        finally_flag = true;       
    }
    else {
        if (top_moist > top_max_moist && infil_capacity >= max_infil_capacity) {
            top_moist = top_max_moist;
            infil_capacity = max_infil_capacity;
            runoff_sat = inflow;
            runoff_infil = 0.0;
            tmp_infil = 0.0;
            finally_flag = true;       
        }
        else {
            ex = 1.0 / (1.0 + b_infilt);
            A = 1.0 - (1.0 - pow((top_moist / top_max_moist), ex));
            infil_capacity = max_infil_capacity * A;
            if (inflow + infil_capacity > max_infil_capacity) {
                if (max_infil * iter_dt >= inflow) {
                    tmp_depth = max_infil_capacity - infil_capacity;
                    tmp_sat = 0.0;
                    calc_sat_runoff(infil_capacity, max_infil_capacity, tmp_depth, b_infilt, &tmp_sat);
                    tmp_var1 = max_infil_capacity - infil_capacity - tmp_sat - (i_0 * iter_dt) * 
                              (1.0 - (1.0 - pow((inflow - tmp_sat) / (max_infil * iter_dt), b_dynamic + 1.0)));
                    if (tmp_var1 <= 0.0) {
                        tmp_depth = max_infil_capacity - infil_capacity;
                        tmp_infil = top_max_moist - top_moist;
                        runoff_sat = inflow - tmp_infil;
                        runoff_infil = 0.0;
                        top_moist = top_max_moist;
                        infil_capacity = max_infil_capacity;
                        finally_flag = true;
                    }
                    else {
                        tmp_depth = 0.0;
                        iter = 0;
                        do {
                            last_depth = tmp_depth;
                            tmp_sat = 0.0;
                            calc_sat_runoff(infil_capacity, max_infil_capacity, tmp_depth, b_infilt, &tmp_sat);
                            tmp_depth = tmp_sat + ((i_0 * iter_dt) * (1.0 - (1.0 - pow((inflow - tmp_sat) / 
                                                            (max_infil * iter_dt), b_dynamic + 1.0))));
                            iter++;
                        }
                        while (fabs(tmp_depth - last_depth) > tol_error && 
                                   iter < last_iter);
                        begin_flag = true;
                    }
                }
                else {
                    tmp_sat = 0.0;
                    calc_sat_runoff(infil_capacity, max_infil_capacity, tmp_depth, b_infilt, &tmp_sat);
                    if (tmp_sat + max_infil * iter_dt <= inflow) {
                        if (max_infil_capacity - infil_capacity - tmp_sat - (max_infil * iter_dt) <= 0.0) {
                            tmp_depth = max_infil_capacity - infil_capacity;
                            tmp_infil = top_max_moist - top_moist;
                            runoff_sat = inflow - tmp_infil;
                            runoff_infil = 0.0;
                            top_moist = top_max_moist;
                            infil_capacity = max_infil_capacity;
                            finally_flag = true;
                        }
                        else {
                            tmp_depth = 0.0;
                            iter = 0;
                            do {
                                last_depth = tmp_depth;
                                tmp_sat = 0.0;
                                calc_sat_runoff(infil_capacity, max_infil_capacity, tmp_depth, b_infilt, &tmp_sat);
                                tmp_depth = tmp_sat + (i_0 * iter_dt);
                                iter++;
                            }
                            while (fabs(tmp_depth - last_depth) > tol_error && 
                                        iter < last_iter);
                            begin_flag = true;
                        }
                    }
                    else {
                        tmp_depth = inflow / 2.0;
                        iter = 0;
                        do {
                            last_depth = tmp_depth;
                            tmp_sat = 0.0;
                            calc_sat_runoff(infil_capacity, max_infil_capacity, tmp_depth, b_infilt, &tmp_sat);
                            tmp_depth -= tmp_sat - (i_0 * iter_dt) + inflow;
                            if (tmp_depth < 0.0) {
                                tmp_depth = 0.0;
                            }
                            if (tmp_depth > inflow) {
                                tmp_depth = inflow;
                            }
                            iter++;
                        }
                        while (fabs(tmp_depth - last_depth) > tol_error && 
                                    iter < last_iter); 
                        init_depth = tmp_depth;
                        iter = 0;
                        do {
                            last_depth = tmp_depth;
                            tmp_sat = 0.0;
                            tmp_infil = 0.0;
                            calc_sat_runoff(infil_capacity, max_infil_capacity, 
                                            tmp_depth, b_infilt, &tmp_sat);
                            calc_infil_runoff(tmp_depth, init_depth, iter_dt, inflow, 
                                              &tmp_sat, max_infil, b_dynamic, &tmp_infil);
                            tmp_depth = inflow - tmp_infil;
                            iter++;
                        }
                        while (fabs(tmp_depth - last_depth) > tol_error && 
                                    iter < last_iter);
                        begin_flag = true;
                    }
                }
                if (begin_flag) {
                    if (tmp_depth < 0.0) {
                        tmp_depth = 0.0;
                    }
                    if (tmp_depth > inflow) {
                        tmp_depth = inflow;
                    }
                    calc_sat_runoff(infil_capacity, max_infil_capacity, tmp_depth, b_infilt, &tmp_sat1);
                    runoff_sat = tmp_sat1;
                    runoff_infil = inflow - tmp_depth;
                    tmp_infil = tmp_depth - runoff_sat;
                    top_moist += tmp_infil;
                    tmp_depth += infil_capacity;
                    if (top_moist < 0.0) {
                        top_moist = 0.0;
                    }
                    if (top_moist > top_max_moist) {
                        top_moist = top_max_moist;
                    }
                    ex = 1.0 / (1.0 + b_infilt);
                    A = 1.0 - (1.0 - pow((top_moist / top_max_moist), ex));
                    infil_capacity = max_infil_capacity * A;
                    finally_flag = true;
                }      
            }
            else {
                if (max_infil * iter_dt >= inflow) {
                    tmp_depth = inflow / 2.0;
                    iter = 0;
                    do {
                        last_depth = tmp_depth;
                        tmp_sat = 0.0;
                        calc_sat_runoff(infil_capacity, max_infil_capacity, tmp_depth, b_infilt, &tmp_sat);
                        tmp_depth = tmp_sat + (i_0 * iter_dt) * (1.0 - (1.0 - pow((inflow - tmp_sat) / 
                                                        (max_infil * iter_dt), b_dynamic + 1.0)));
                        iter++;
                    }
                    while (fabs(tmp_depth - last_depth) > tol_error &&
                                iter < last_iter);
                    middle_flag = true;
                }
                else {
                    tmp_sat = 0.0;
                    calc_sat_runoff(infil_capacity, max_infil_capacity, tmp_depth, b_infilt, &tmp_sat);
                    if (tmp_sat + (max_infil * iter_dt) <= inflow) {
                        tmp_depth = inflow / 2.0;
                        iter = 0;
                        do {
                            last_depth = tmp_depth;
                            tmp_sat = 0.0;
                            calc_sat_runoff(infil_capacity, max_infil_capacity, tmp_depth, b_infilt, &tmp_sat);
                            tmp_depth = tmp_sat + (i_0 * iter_dt);
                            iter++;
                        }
                        while (fabs(tmp_depth - last_depth) > tol_error && 
                                    iter < last_iter);
                        middle_flag = true;          
                    }
                    else {
                        tmp_depth = 0.0;
                        iter = 0;
                        do {
                            last_depth = tmp_depth;
                            tmp_sat = 0.0;
                            calc_sat_runoff(infil_capacity, max_infil_capacity, tmp_depth, b_infilt, &tmp_sat);
                            tmp_depth += (inflow - max_infil * iter_dt) - tmp_sat;
                            if (tmp_depth < 0.0) {
                                tmp_depth = 0.0;
                            }
                            if (tmp_depth > inflow) {
                                tmp_depth = inflow;
                            }
                            tmp_sat = 0.0;
                            calc_sat_runoff(infil_capacity, max_infil_capacity, tmp_depth, b_infilt, &tmp_sat);
                            iter++;
                        }
                        while (fabs(tmp_sat + (max_infil * iter_dt) - inflow) > tol_error &&
                                    iter < last_iter);
                        init_depth = tmp_depth;
                        iter = 0;
                        do {
                            last_depth = tmp_depth;
                            tmp_sat = 0.0;
                            tmp_infil = 0.0;
                            calc_sat_runoff(infil_capacity, max_infil_capacity,
                                            tmp_depth, b_infilt, &tmp_sat);
                            calc_infil_runoff(tmp_depth, init_depth, iter_dt, 
                                              inflow, &tmp_sat, max_infil, 
                                              b_dynamic, &tmp_infil);
                            tmp_depth = inflow - tmp_infil;
                            iter++;
                        }
                        while (fabs(tmp_depth - last_depth) > tol_error && 
                                    iter < last_iter);
                        middle_flag = true;
                    }
                }
                if (middle_flag) {
                    if (tmp_depth < 0.0) {
                        tmp_depth = 0.0;
                    }
                    if (tmp_depth > inflow) {
                        tmp_depth = inflow;
                    }
                    tmp_sat1 = 0.0;
                    calc_sat_runoff(infil_capacity, max_infil_capacity, tmp_depth, b_infilt, &tmp_sat1);
                    runoff_sat = tmp_sat1;
                    runoff_infil = inflow - tmp_depth;
                    tmp_infil = tmp_depth - runoff_sat;
                    top_moist += tmp_infil;
                    if (top_moist < 0) {
                        top_moist = 0;
                    }
                    if (top_moist > top_max_moist) {
                        top_moist = top_max_moist;
                    }
                    ex = 1.0 / (1.0 + b_infilt);
                    A = 1.0 - (1.0 - pow((top_moist / top_max_moist), ex));
                    infil_capacity = max_infil_capacity * A;
                    finally_flag = true;
                }
            }
        }
    }
    if (finally_flag) {
        *runoff = (runoff_sat + runoff_infil) / iter_dt;
        if (*runoff > dt_inflow) {
            *runoff = dt_inflow;
        }
        if (*runoff < 0) {
            *runoff = 0;
        }
        *infil_rate = dt_inflow - *runoff;
        cell->asat = A;
    }
}

/******************************************************************************
* @brief    Calculate the saturated area and runoff
******************************************************************************/
void
calc_sat_runoff(double    infil_capacity,
                double    max_infil_capacity,
                double    tmp_depth,
                double    b_infilt,
                double   *runoff_sat)
{
    double  zwt;

    zwt = infil_capacity + tmp_depth;
    if (zwt > max_infil_capacity) {
        zwt = max_infil_capacity;
    }

    *runoff_sat = tmp_depth - (max_infil_capacity / (b_infilt + 1.0)) * 
                    ((pow(1.0 - (infil_capacity / max_infil_capacity), b_infilt + 1.0)) - 
                    (pow(1.0 - (zwt / max_infil_capacity), b_infilt + 1.0)));
    if (*runoff_sat < 0.0) {
        *runoff_sat = 0.0;
    }
}

/******************************************************************************
* @brief    Calculate the saturated area and runoff
******************************************************************************/
void
calc_infil_runoff(double    tmp_depth,
                  double    init_depth,
                  double    step_dt,
                  double    inflow,
                  double   *runoff_sat,
                  double    max_infil,
                  double    b_dynamic,
                  double   *runoff_infil)
{
    if (tmp_depth >= init_depth) {
        *runoff_infil = inflow - *runoff_sat - (max_infil * step_dt) * 
                        (1.0 - pow((1.0 - (inflow - *runoff_sat) / 
                        (max_infil * step_dt)), b_dynamic + 1.0));
    }
    else {
        *runoff_infil = inflow - *runoff_sat - (max_infil * step_dt);
    }

    if (*runoff_infil < 0.0) {
        *runoff_infil = 0.0;
    }
}

/******************************************************************************
* @brief   Compute  soil surface infiltration rate based on Green-Ampt equa.
******************************************************************************/
void
calc_soil_infil(double       moist,
                double       Wsat,
                double       capil_drive,
                double       dz_node,
                double       Ksat,
                double       inflow,
                double       conductivity,
                double      *i_0,
                double      *infil_acc)
{

    // Maximum infiltrability based on the Eq. 6.25. (m/s)
    double tmp_infil = capil_drive * max(0.01, Wsat - moist) * dz_node;
    (*i_0) = Ksat + (tmp_infil / *infil_acc) * (Ksat - conductivity);

    if (Ksat < inflow) {
        (*i_0) = min(inflow, *i_0);
    }
    else {
        (*i_0) = inflow;
    }
    (*infil_acc) += *i_0;
}
