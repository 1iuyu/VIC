/******************************************************************************
* @section DESCRIPTION
*
* This routine computes all surface fluxes, and solves the snow accumulation
* and ablation algorithm. Solutions are for the current snow band and
* vegetation type.
******************************************************************************/

#include "vic_run.h"

/******************************************************************************
* @brief        This routine computes all surface fluxes.
******************************************************************************/
int
surface_fluxes(double             step_dt,
               force_data_struct *force,
               energy_bal_struct *energy,
               cell_data_struct  *cell,
               snow_data_struct  *snow,
               soil_con_struct   *soil_con,
               veg_var_struct    *veg_var,
               veg_lib_struct    *veg_lib)
{
    extern parameters_struct param;
    int    ErrorFlag;
    
    // Structures holding values for current iteration
    cell_data_struct  iter_cell;
    energy_bal_struct iter_energy;
    veg_var_struct    iter_veg_var;
    snow_data_struct  iter_snow;
    double fScale[MAX_STATES] = {0};
    double xScale[MAX_STATES] = {0};
    double mat_RHS[MAX_STATES] = {0};
    double mat_LHS[MAX_STATES*10] = {0};
    double shortwave_dir[MAX_SWBANDS] = {0};
    double shortwave_dfs[MAX_SWBANDS] = {0};

    /*******************************************
       Set-up sub-time step controls
    *******************************************/
    size_t iter = 0;
    size_t nidx;
    size_t i, j, k;
    size_t col, row;
    size_t jmin, jmax;
    size_t Nnode = cell->Nnode;
    size_t Nstate = Nnode * 2;
    int IPIV[MAX_STATES] = {0};
    double grad[MAX_STATES] = {0};
    double mat_state[MAX_STATES] = {0};
    double time_accum = 0.0;
    double iter_dt = step_dt;
    double shortwave = force->shortwave;

    /**************************************
      计算直射和漫射辐射通量
    **************************************/
    shortwave_dir[BAND_VIS] = shortwave * param.RAD_DIR_F * param.RAD_VIS_F;                 // 直射-可见光
    shortwave_dir[BAND_NIR] = shortwave * param.RAD_DIR_F * (1.0 - param.RAD_VIS_F);         // 直射-近红外
    shortwave_dfs[BAND_VIS] = shortwave * (1.0 - param.RAD_DIR_F) * param.RAD_VIS_F;         // 漫射-可见光
    shortwave_dfs[BAND_NIR] = shortwave * (1.0 - param.RAD_DIR_F) * (1.0 - param.RAD_VIS_F); // 漫射-近红外

    /*******************************
      Advected heat flux from Prec
    *******************************/
    AdvectedEnergy(force, energy, 
                   cell, veg_var);
    
    /*******************************
      Update surface water fraction
    *******************************/
    calc_surf_water(step_dt, cell,
                    snow, soil_con);
                    
    /***************************
      Surface shortwave albedo
    ***************************/
    surface_albedo(step_dt, force, 
                   energy, 
                   cell, snow, 
                   soil_con, 
                   veg_var, veg_lib);

    /***************************
      Surface radiative fluxes
    ***************************/
    surface_radiation(shortwave_dir,
                      shortwave_dfs,
                      energy, cell,
                      snow, veg_var);

    /******************************
      Compute longwave emissivity
    ******************************/
    surface_emissivity(energy, cell,
                       snow, veg_var);
                  
    /**************************************************
       Begin iterative solution of surface fluxes
    **************************************************/
    do {
        // 将当前迭代的结构体值赋给临时结构体，以便在迭代过程中更新
        iter_cell = (*cell);
        iter_energy = (*energy);
        iter_veg_var = (*veg_var);
        iter_snow = (*snow);

        // 尝试找到一个可接受的步长
        while (!iter_energy.energy_flag || !iter_energy.moist_flag) {

            for (i = 0, j = 0; i < Nnode; i++, j+=2) {
                if (i < iter_snow.Nsnow) {
                    mat_state[j] = iter_snow.pack_T[i];
                    mat_state[j+1] = iter_snow.theta_moist[i];
                }
                else {
                    nidx = i - iter_snow.Nsnow;
                    mat_state[j] = iter_cell.soil_T[nidx];
                    mat_state[j+1] = iter_cell.moist[nidx];
                }
            }

            /************************************
              Update energy balance variables
            ************************************/
            update_fluxes(step_dt, force,
                          &iter_energy, 
                          &iter_cell,
                          &iter_snow, soil_con,
                          &iter_veg_var, veg_lib);
            
            /*******************************************
                Solve Water Balance Components
            ********************************************/
            ErrorFlag = calc_resid_vec(iter_dt, mat_RHS,
                                       &iter_energy, 
                                       &iter_cell, 
                                       &iter_snow,
                                       soil_con);

            if (ErrorFlag == ERROR) {
                return (ERROR);
            }

            /**********************************************
                Solve Energy Balance Components
            **********************************************/
            ErrorFlag = calc_Jacob_mat(iter_dt, mat_LHS,
                                       &iter_energy,
                                       &iter_cell,
                                       &iter_snow, soil_con);

            if (ErrorFlag == ERROR) {
                return (ERROR);
            }

            /*******************************************************
                Scale residual and Jacobian before Newton solve
            *******************************************************/
            for (i = 0, j = 0; i < Nnode; i++, j+=2) {
                // Energy equation
                fScale[j] = 1.0e-6;
                xScale[j] = 1.0;
                // Water equation
                fScale[j+1] = 100.0;
                xScale[j+1] = 1.0;
            }

            /* Scale RHS: R_s = F * R */
            double fOld = 0.0;
            for (i = 0; i < Nstate; i++) {
                mat_RHS[i] *= fScale[i];
                fOld += 0.5 * mat_RHS[i] * mat_RHS[i];
            }

            /* Scale band Jacobian: J_s = F * J * X */
            for (col = 0; col < Nstate; col++) {
                if (col > LAPACK_KU) {
                    jmin = col - LAPACK_KU;
                } else {
                    jmin = 0;
                }
                jmax = min(Nstate - 1, col + LAPACK_KL);

                for (row = jmin; row <= jmax; row++) {
                    k = get_band_idx(row, col);
                    mat_LHS[k] *= fScale[row] * xScale[col];
                }
            }

            // 计算函数向量的梯度
            calc_gradient(mat_LHS, mat_RHS, grad, Nstate);

            // 调用LAPACK函数求解线性方程组
            int info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, Nstate,
                                     LAPACK_KL, LAPACK_KU, 1, 
                                     mat_LHS, LAPACK_LDAB, IPIV, mat_RHS, Nstate);

            if (info != 0) {
                return ERROR;
            }

            // 计算初始斜率
            double slopeInit = 0.0;
            for (i = 0; i < Nstate; i++) {
                slopeInit += grad[i] * mat_RHS[i];
            }

            // Scale solution: X_s = X * X_scale
            for (i = 0; i < Nstate; i++) {
                mat_RHS[i] *= xScale[i];
            }

            /* Refine Newton increment */
            ErrorFlag = physical_constraints(Nstate, iter_dt, 
                                             fOld, slopeInit,
                                             fScale, mat_RHS, 
                                             mat_state, force, 
                                             &iter_energy, 
                                             &iter_cell, &iter_snow, 
                                             soil_con, 
                                             &iter_veg_var, veg_lib);

            if (ErrorFlag == ERROR) {
                return (ERROR);
            }

            double diff_t = 0.0, diff_q = 0.0;
            double max_energy = 0.0, max_moist = 0.0;
            // 物理约束能量向量
            for (i = 0, j = 0; i < Nnode; i++, j+=2) {
                if (i < iter_snow.Nsnow) {
                    diff_t = iter_snow.pack_T[i] - mat_state[j];
                    diff_q = iter_snow.theta_moist[i] - mat_state[j+1];
                }
                else {
                    nidx = i - iter_snow.Nsnow;
                    diff_t = iter_cell.soil_T[nidx] - mat_state[j];
                    diff_q = iter_cell.moist[nidx] - mat_state[j+1];
                }

                if (fabs(diff_t) > max_energy) {
                    max_energy = fabs(diff_t);
                }
                if (fabs(diff_q) > max_moist) {
                    max_moist = fabs(diff_q);
                }
            }

            if (max_energy < 0.01 && max_moist < 0.01) {
                iter_energy.energy_flag = true;
                iter_energy.moist_flag = true;
            }

            iter++;
            // 迭代控制
            if (iter >= param.MAX_ITER_OVER) {
                if (iter_dt > MAX_SUBDAILY_STEPS_PER_DAY * SEC_PER_DAY) {
                    // 时间步长减半，重新开始
                    iter_dt = iter_dt / 2.0;
                    iter_energy.energy_flag = false;
                    iter_energy.moist_flag = false;
                    for (i = 0; i < MAX_SNOWS; i++) {
                        iter_snow.phase_snow[i] = 0.0;
                        iter_snow.pack_frze[i] = 0.0;
                        iter_snow.pack_melt[i] = 0.0;
                    }
                    iter = 0;
                    // 回滚状态至初始值
                    iter_cell = (*cell);
                    iter_energy = (*energy);
                    iter_veg_var = (*veg_var);
                    iter_snow = (*snow);
                    continue;
                } 
                else {
                    // 已达到最小时间步长，强制收敛并继续
                    iter_energy.energy_flag = true;
                    iter_energy.moist_flag = true;
                    break;
                }
            }
        }

        // Update the original structures with the 
        // converged values from the last iteration.
        (*cell) = iter_cell;
        (*energy) = iter_energy;
        (*veg_var) = iter_veg_var;
        (*snow) = iter_snow;

        // 应用最终相变
        apply_phase_change(iter_dt, step_dt, energy, snow);
        time_accum += iter_dt;
        /* 保存当前子时间步结束状态 */
        update_last_state(energy, cell, snow);
        
    }
    while(time_accum < step_dt);

    return (0);
}