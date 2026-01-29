/******************************************************************************
* @section DESCRIPTION
*
* Calculate infiltration and runoff from the surface, gravity driven drainage
* between all soil layers, and generates baseflow from the bottom layer.
******************************************************************************/

#include <vic_run.h>

/******************************************************************************
* @brief    Calculate infiltration and runoff from the surface, gravity driven
*           drainage between all soil layers, and generates baseflow from the
*           bottom layer.
******************************************************************************/
int
runoff(double             ppt,
       cell_data_struct  *cell,
       soil_con_struct   *soil_con)
{
    extern option_struct       options;
    extern global_param_struct global_param;
    extern parameters_struct   param;

    size_t      i, lidx;
    size_t      time_step;
    size_t      last_iter;
    size_t      iter_dt;
    size_t      iter;
    size_t      runoff_steps_per_dt;
    double      tmp_moist;
    double      infil_rate;
    double      dt_inflow;
    double      min_liq;
    double      mat_A[MAX_NODES];   // lower diag
    double      mat_B[MAX_NODES];   // main diag
    double      mat_C[MAX_NODES];   // upper diag
    double      mat_RHS[MAX_NODES]; // RHS (source).
    double      tmp_ice[MAX_NODES];
    double      tmp_liq[MAX_NODES];
    // 指针赋值 
    double *liq = cell->liq;
    double *ice = cell->ice;
    double *moist = cell->moist;
    double *dz_soil = soil_con->dz_soil;
    double *Wsat_node = soil_con->Wsat_node;
    double *soil_frost = cell->soil_frost;
    double *porosity = soil_con->porosity_node;
    // 初始化
    size_t Nnode = options.Nnode;
    double infil_acc = 1.0e-06;
    double sat_excess = 0.0;
    double runoff = 0.0;
    double baseflow = 0.0;
    double dt_soil_drain = 0.0;
    double dt_runoff = 0.0;
    double soil_drain = 0.0;
    runoff_steps_per_dt = global_param.runoff_steps_per_day /
                          global_param.model_steps_per_day;
    for (i = 0; i < MAX_NODES; i++) {
        tmp_ice[i] = 0.0;
        tmp_liq[i] = 0.0;
    }
    for (i = 0; i < MAX_NODES; i++) {
        mat_A[i] = 0.0;
        mat_B[i] = 0.0;
        mat_C[i] = 0.0;
        mat_RHS[i] = 0.0;
    }
    
    /**************************************************
       Initialize Variables
    **************************************************/
    for (i = 0; i < Nnode; i++) {
        porosity[i] = max(1.0e-4, Wsat_node[i] - ice[i]);
        /** Set Layer Liquid Moisture Content **/
        liq[i] = moist[i] - ice[i];
        if (liq[i] > porosity[i]) {
            liq[i] = porosity[i];
        }
        sat_excess += max(0.0, liq[i] - porosity[i]) * dz_soil[i];
    }                        
    // impermeable fraction due to frozen soil
    for (i = 0; i < Nnode; i++) {
        tmp_ice[i] = min(1.0, ice[i] / Wsat_node[i]);
        soil_frost[i] = max(0.0, exp(-param.SOIL_FROST * 
                        (1.0 - tmp_ice[i])) - exp(-param.SOIL_FROST)) / 
                                (1.0 - exp(-param.SOIL_FROST));
    }

    double max_frost = 0.0;
    for (i = 0; i < Nnode; i++) {
        if (soil_frost[i] > max_frost) {
            max_frost = soil_frost[i];
        }
    }

    /** ppt = amount of liquid water coming to the surface [m/s] **/
    infil_rate = 0.;
    dt_inflow = ppt;
    calc_dynamic_runoff(dt_inflow, &infil_acc,     
                        global_param.model_steps_per_day,
                        &runoff, &infil_rate,
                        cell, soil_con);
    double runoff_dt = global_param.runoff_dt;
    /** determine iteration times  to solve soil water diffusion and moisture **/
    if (infil_rate * runoff_dt > Wsat_node[0] * dz_soil[0]) {

        last_iter = param.MAX_ITER_INFLI_RUNOFF * 2;
        iter_dt = SEC_PER_DAY / (global_param.model_steps_per_day *
                                param.MAX_ITER_INFLI_RUNOFF * 2);
    }
    else {
        last_iter = param.MAX_ITER_INFLI_RUNOFF;
        iter_dt = SEC_PER_DAY / (global_param.model_steps_per_day *
                                param.MAX_ITER_INFLI_RUNOFF);                          
    }

    /**************************************************
       Compute Flow Between Soil Layers
    **************************************************/
    infil_acc = 1.0e-06;
    runoff = 0.0;

    for (time_step = 0; time_step < runoff_steps_per_dt; time_step++) {

        for (iter = 0; iter < last_iter; iter++) {
            if (dt_inflow > 0.0) {
                calc_dynamic_runoff(dt_inflow, &infil_acc,
                                   iter_dt, &dt_runoff,
                                   &infil_rate, cell, soil_con);
            }
            /** Compute the right hand side of the time tendency term of the soil
                water diffusion equation. prepare the matrix coefficients for the
                tri-diagonal matrix of the implicit time scheme. **/
            calc_soilmoist_Richards(mat_A, mat_B,
                                    mat_C, mat_RHS,
                                    infil_rate,
                                    &dt_soil_drain,
                                    cell, soil_con);
            /**  Compute soil moisture content using based on 
                 Richards diffusion & tri-diagonal matrix **/
            solve_soilmoist_Richards(mat_A, mat_B,
                                     mat_C, mat_RHS,
                                     iter_dt, Wsat_node, liq, ice,
                                     moist, dz_soil,
                                     &cell->soil_excess,
                                     porosity);
            runoff += dt_runoff;
            soil_drain += dt_soil_drain;
            sat_excess += cell->soil_excess;
        }
        runoff /= (double) last_iter;
        runoff *= MM_PER_M + sat_excess * MM_PER_M / runoff_dt;  // convert to mm/s
        soil_drain /= (double) last_iter;
        soil_drain *= MM_PER_M;  // convert to mm/s

        /** Adjust liquid soil moisture to account for minimum moisture **/
/*      for (lidx = 0; lidx < Nnode; lidx++) {
            tmp_liq[lidx] = liq[lidx] * dz_node[lidx] * MM_PER_M;
        }
        min_liq = 0.01;  // mm
        for (lidx = 0; lidx < Nnode - 1; lidx++) {
            if (tmp_liq[lidx] < 0.0) {
                tmp_moist = min_liq - tmp_liq[lidx];
            }
            else {
                tmp_moist = 0.0;
            }
            tmp_liq[lidx] += tmp_moist;
            tmp_liq[lidx + 1] -= tmp_moist;
        }
        if (tmp_liq[Nnode - 1] < min_liq) {
            tmp_moist = min_liq - tmp_liq[Nnode - 1];
        }
        else {
            tmp_moist = 0.0;
        }
        tmp_liq[Nnode - 1] += tmp_moist;
        baseflow -= tmp_moist / runoff_dt;

        for (lidx = 0; lidx < Nnode; lidx++) {
            liq[lidx] = tmp_liq[lidx] / (dz_node[lidx] * MM_PER_M);
        }

        baseflow += soil_drain; */

    } /* end of sub-dt time step loop */

    /************************************************
      Compute Baseflow
    ************************************************/
    baseflow = compute_zwt(runoff_dt, 
                           max_frost, 
                           cell, soil_con); // mm/s

    /** Store tile-wide values **/
    for (lidx = 0; lidx < Nnode; lidx++) {
        moist[lidx] = (liq[lidx] + ice[lidx]);
    }
    //cell->asat += A * frost_fract[fidx];
    cell->runoff = runoff * runoff_dt;

    cell->baseflow = baseflow * runoff_dt;


    return (0);
}

/******************************************************************************
* @brief    Calculate the saturated area and runoff
******************************************************************************/
void
calc_dynamic_runoff(double            dt_inflow,
                    double           *infil_acc,
                    size_t            time_step,
                    double           *runoff,
                    double           *infil_rate,
                    cell_data_struct *cell,
                    soil_con_struct  *soil_con)
{
    extern option_struct options;

    size_t  lidx;
    size_t  iter;
    bool    begin_flag = false;
    bool    middle_flag = false;
    bool    finally_flag = false;
    double  A, ex, dt, i_0;
    double  max_infil;
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
    double tol_error = 0.5 * HOURS_PER_DAY / (double) time_step / MM_PER_M; // 0.5 mm per hour time step
    double cap_drive = soil_con->cap_drive;
    double b_infilt = soil_con->b_infilt;
    double b_dynamic = soil_con->b_dynamic;
    double *b_exp = soil_con->bexp_node;

    for (lidx = 0; lidx < options.Nnode - 7; lidx++) {
        top_moist += moist[lidx] * dz_soil[lidx];
        top_max_moist += Wsat_node[lidx] * dz_soil[lidx];
    }
    if (top_moist > top_max_moist) {
        top_moist = top_max_moist;
    }

    /** A as in Wood et al. in JGR 97, D3, 1992 equation (1) **/
    ex = 1.0 / (1.0 + b_infilt);
    A = 1.0 - (1.0 - pow((top_moist / top_max_moist), ex));

    max_infil_capacity = (1.0 + b_infilt) * top_max_moist;
    infil_capacity = max_infil_capacity * A;
    /* precipitation depth [m] */
    dt = SEC_PER_DAY / time_step;
    double inflow = dt_inflow * dt;

    /** compute surface infiltration **/
    calc_soil_infil(moist[0], Wsat_node[0],
                    b_exp[0], cap_drive,
                    dz_soil[0],
                    soil_con->Ksat_node[0],
                    inflow, &i_0, infil_acc);

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
                if (max_infil * dt >= inflow) {
                    tmp_depth = max_infil_capacity - infil_capacity;
                    tmp_sat = 0.0;
                    calc_sat_runoff(infil_capacity, max_infil_capacity, tmp_depth, b_infilt, &tmp_sat);
                    tmp_var1 = max_infil_capacity - infil_capacity - tmp_sat - (i_0 * dt) * 
                              (1.0 - (1.0 - pow((inflow - tmp_sat) / (max_infil * dt), b_dynamic + 1.0)));
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
                            tmp_depth = tmp_sat + ((i_0 * dt) * (1.0 - (1.0 - pow((inflow - tmp_sat) / 
                                                            (max_infil * dt), b_dynamic + 1.0))));
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
                    if (tmp_sat + max_infil * dt <= inflow) {
                        if (max_infil_capacity - infil_capacity - tmp_sat - (max_infil * dt) <= 0.0) {
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
                                tmp_depth = tmp_sat + (i_0 * dt);
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
                            tmp_depth -= tmp_sat - (i_0 * dt) + inflow;
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
                            calc_infil_runoff(tmp_depth, init_depth, dt, inflow, 
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
                if (max_infil * dt >= inflow) {
                    tmp_depth = inflow / 2.0;
                    iter = 0;
                    do {
                        last_depth = tmp_depth;
                        tmp_sat = 0.0;
                        calc_sat_runoff(infil_capacity, max_infil_capacity, tmp_depth, b_infilt, &tmp_sat);
                        tmp_depth = tmp_sat + (i_0 * dt) * (1.0 - (1.0 - pow((inflow - tmp_sat) / 
                                                        (max_infil * dt), b_dynamic + 1.0)));
                        iter++;
                    }
                    while (fabs(tmp_depth - last_depth) > tol_error &&
                                iter < last_iter);
                    middle_flag = true;
                }
                else {
                    tmp_sat = 0.0;
                    calc_sat_runoff(infil_capacity, max_infil_capacity, tmp_depth, b_infilt, &tmp_sat);
                    if (tmp_sat + (max_infil * dt) <= inflow) {
                        tmp_depth = inflow / 2.0;
                        iter = 0;
                        do {
                            last_depth = tmp_depth;
                            tmp_sat = 0.0;
                            calc_sat_runoff(infil_capacity, max_infil_capacity, tmp_depth, b_infilt, &tmp_sat);
                            tmp_depth = tmp_sat + (i_0 * dt);
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
                            tmp_depth += (inflow - max_infil * dt) - tmp_sat;
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
                        while (fabs(tmp_sat + (max_infil * dt) - inflow) > tol_error &&
                                    iter < last_iter);
                        init_depth = tmp_depth;
                        iter = 0;
                        do {
                            last_depth = tmp_depth;
                            tmp_sat = 0.0;
                            tmp_infil = 0.0;
                            calc_sat_runoff(infil_capacity, max_infil_capacity,
                                            tmp_depth, b_infilt, &tmp_sat);
                            calc_infil_runoff(tmp_depth, init_depth, dt, 
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
        *runoff = (runoff_sat + runoff_infil) / dt;
        if (*runoff > dt_inflow) {
            *runoff = dt_inflow;
        }
        if (*runoff < 0) {
            *runoff = 0;
        }
        *infil_rate = dt_inflow - *runoff;
    }
}

/******************************************************************************
* @brief    Calculate the saturated area and runoff
******************************************************************************/
void
calc_soilmoist_Richards(double           *mat_A,
                        double           *mat_B,
                        double           *mat_C,
                        double           *mat_RHS,
                        double            infil_rate,
                        double           *soil_drain,
                        cell_data_struct *cell,
                        soil_con_struct  *soil_con)
{
    extern option_struct options;

    size_t  i, Nnode;
    double  soil_drain_slope;           // slope index for soil drainage      
    double  rel_sat[MAX_NODES];
    double  grad_theta[MAX_NODES];
    double  WaterExcess[MAX_NODES];
    double  depthInv[MAX_NODES];
    double *moist = cell->moist;
    double *Ksat = soil_con->Ksat_node;
    double *Dsat = soil_con->Dsat_node;
    double *Wsat_node = soil_con->Wsat_node;
    double *bexp = soil_con->bexp_node;
    double *dz_soil = soil_con->dz_soil;
    double *Zsum_soil = soil_con->Zsum_soil;
    double *soil_transp = cell->soil_transp;
    double *soil_frost = cell->soil_frost;
    double *diffusivity = cell->diffusivity;    // soil water diffusivity [m2/s]
    double *conductivity = cell->conductivity;   // soil hydraulic conductivity [m/s]

    Nnode = options.Nnode;
    soil_drain_slope = 0.5;

    for (i = 0; i < MAX_NODES; i++) {
        rel_sat[i] = 0.0;
        depthInv[i] = 0.0;
        grad_theta[i] = 0.0;
        WaterExcess[i] = 0.0;
    }

    /** compute soil hydraulic conductivity and diffusivity **/
    for (i = 0; i < Nnode; i++) {
        rel_sat[i] = max(0.01, moist[i] / Wsat_node[i]);
        /** soil water diffusivity **/
        diffusivity[i] = Dsat[i] * pow(rel_sat[i], bexp[i] + 2.0);
        diffusivity[i] *= (1.0 - soil_frost[i]);
        /** soil hydraulic conductivity **/
        conductivity[i] = Ksat[i] * pow(rel_sat[i], 2.0 * bexp[i] + 3.0);
        conductivity[i] *= (1.0 - soil_frost[i]);
    }

    /*************************************
       Compute Drainage between Sublayers
    *************************************/
    // 构造水分梯度项
    for (i = 0; i < Nnode; i++) {
        if (i == 0) {
            depthInv[i] = 2.0 / (dz_soil[i] + dz_soil[i + 1]);
            grad_theta[i] = 2.0 * (moist[i] - moist[i + 1]) /
                                           (dz_soil[i] + dz_soil[i + 1]);

            WaterExcess[i] = diffusivity[i] * grad_theta[i] + conductivity[i] -
                                  infil_rate + soil_transp[i] + cell->esoil;
        }
        else if (i < Nnode - 1) {
            depthInv[i] = 2.0 / (dz_soil[i] + dz_soil[i + 1]);
            grad_theta[i] = 2.0 * (moist[i] - moist[i + 1]) /
                                            (dz_soil[i] + dz_soil[i + 1]);

            WaterExcess[i] = diffusivity[i] * grad_theta[i] + conductivity[i] - 
                                  diffusivity[i - 1] * grad_theta[i - 1] - conductivity[i - 1] +
                                  soil_transp[i];
        }
        else {
            (*soil_drain) = conductivity[i] * soil_drain_slope;
            WaterExcess[i] = -diffusivity[i - 1] * grad_theta[i - 1] - conductivity[i - 1]
                             + soil_transp[i] + (*soil_drain);
        }

    }
    // 构建三对角矩阵 A (下), B (中), C (上)
    for (i = 0; i < Nnode; i++) {
        if (i == 0) {
            mat_A[i] = 0.0;
            mat_B[i] = diffusivity[i] * depthInv[i] / dz_soil[i];
            mat_C[i] = -mat_B[i];
        } 
        else if (i < Nnode - 1) {
            mat_A[i] = -diffusivity[i - 1] * depthInv[i - 1] / dz_soil[i];
            mat_C[i] = -diffusivity[i] * depthInv[i] / dz_soil[i];
            mat_B[i] = -(mat_A[i] + mat_C[i]);
        } 
        else {
            mat_A[i] = -diffusivity[i - 1] * depthInv[i - 1] / dz_soil[i];
            mat_C[i] = 0;
            mat_B[i] = -(mat_A[i] + mat_C[i]);
        }
        mat_RHS[i] = -WaterExcess[i] / dz_soil[i];
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
                  double    dt,
                  double    inflow,
                  double   *runoff_sat,
                  double    max_infil,
                  double    b_dynamic,
                  double   *runoff_infil)
{
    if (tmp_depth >= init_depth) {
        *runoff_infil = inflow - *runoff_sat - (max_infil * dt) * 
                        (1.0 - pow((1.0 - (inflow - *runoff_sat) / 
                        (max_infil * dt)), b_dynamic + 1.0));
    }
    else {
        *runoff_infil = inflow - *runoff_sat - (max_infil * dt);
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
                double       bexp,
                double       cap_drive,
                double       dz_node,
                double       Ksat,
                double       inflow,
                double      *i_0,
                double      *infil_acc)
{
    double  conductivity; // soil water conductivity[m/s]
    double  tmp_infil;    // temporary infiltrability variable
    double  rel_sat;      // relative soil moisture

    rel_sat = max(0.01, moist / Wsat);

    /** calculating soil hydraulic conductivity **/
    conductivity = Ksat * pow(rel_sat, 2.0 * bexp + 3.0);

    // Maximum infiltrability based on the Eq. 6.25. (m/s)
    tmp_infil = cap_drive * max(0.01, Wsat - moist) * dz_node;
    (*i_0) = Ksat + (tmp_infil / *infil_acc) * (Ksat - conductivity);

    if (Ksat < inflow) {
        (*i_0) = min(inflow, *i_0);
    }
    else {
        (*i_0) = inflow;
    }
    (*infil_acc) += *i_0;
}

/******************************************************************************
* @brief    Calculate the saturated area and runoff
******************************************************************************/
void
solve_soilmoist_Richards(double   *mat_A,
                         double   *mat_B,
                         double   *mat_C,
                         double   *mat_RHS,
                         size_t    time_step,
                         double   *Wsat,
                         double   *ice,
                         double   *liq,
                         double   *moist,
                         double   *dz_node,
                         double   *soil_excess,
                         double   *eff_porosity)
{
    extern option_struct       options;

    size_t       i, lidx;
    size_t       Nnode;

    Nnode = options.Nnode;

    for (lidx = 0; lidx < Nnode; lidx++) {
        mat_A[lidx] = mat_A[lidx] * (double) time_step;
        mat_B[lidx] = 1.0 + mat_B[lidx] * (double) time_step;
        mat_C[lidx] = mat_C[lidx] * (double) time_step;
        mat_RHS[lidx] = mat_RHS[lidx] * (double) time_step;
    }
    
    tridiag(mat_A, mat_B, mat_C, mat_RHS, Nnode);

    for (lidx = 0; lidx < Nnode; lidx++) {
        liq[lidx] += mat_RHS[lidx];
    }
    
    for (i = Nnode - 1; i > 0; i--) {
        eff_porosity[i] = max(1.0e-4, Wsat[i] - ice[i]);
        *soil_excess = max(0.0, liq[i] - eff_porosity[i]) * dz_node[i];
        liq[i] = min(liq[i], eff_porosity[i]);
        liq[i - 1] += *soil_excess / dz_node[i - 1];

    }
    eff_porosity[0] = max(1.0e-4, Wsat[0] - ice[0]);
    *soil_excess = max(liq[0] - eff_porosity[0], 0.0) * dz_node[0];
    liq[0] = min(eff_porosity[0], liq[0]);
    if (*soil_excess > 0.0) {
        liq[1] += *soil_excess / dz_node[1];
        for (i = 1; i < Nnode - 1; i++) {
            eff_porosity[i] = max(1.0e-4, Wsat[i] - ice[i]);
            *soil_excess = max(liq[i] - eff_porosity[i], 0.0) * dz_node[i];
            liq[i] = min(eff_porosity[i], liq[i]);
            liq[i + 1] += *soil_excess / dz_node[i + 1];
        }
        eff_porosity[Nnode] = max(1.0e-4, Wsat[Nnode] - ice[Nnode]);
        *soil_excess = max(liq[Nnode] - eff_porosity[Nnode], 0.0) * dz_node[Nnode];
        liq[Nnode] = min(eff_porosity[Nnode], liq[Nnode]);
    }

    for (i = 0; i < Nnode; i++) {
        moist[i] = liq[i] + ice[i];
    }

}
