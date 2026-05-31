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
surface_fluxes(size_t             hidx,
               double             step_dt,
               force_data_struct *force,
               energy_bal_struct *energy,
               cell_data_struct  *cell,
               snow_data_struct  *snow,
               soil_con_struct   *soil_con,
               veg_var_struct    *veg_var,
               veg_lib_struct    *veg_lib)
{
    extern option_struct     options;
    extern parameters_struct param;

    int    ErrorFlag;
    size_t i, hidx; // index of initial element of atmos array
    // Structures holding values for current iteration
    cell_data_struct  iter_cell;
    energy_bal_struct iter_energy;
    veg_var_struct    iter_veg_var;
    snow_data_struct  iter_snow;
    double aPAR_sun;
    double aPAR_sha;
    double transmit_direct;
    double transmit_diffuse;
    double tmp_absorb_grnd;
    double ShortOverDir[MAX_SWBANDS];
    double ShortOverDfs[MAX_SWBANDS];
    double shortwave_dir[MAX_SWBANDS];
    double shortwave_dfs[MAX_SWBANDS];

    /*******************************************
       Set-up sub-time step controls
    *******************************************/
    size_t dt_min = 1800; // seconds (half hour)
    size_t iter = 0;
    double time_accum = 0.0;
    double temp_error = 0.01;    // 0.01K
    double moist_error = 0.0001; // 0.1mm
    /***************************
       Compute surface fluxes
    ***************************/
    double coszen = force->coszen[hidx];
    double snowfall = force->snowf[hidx];
    double rainfall = force->rainf[hidx];
    double pressure = force->pressure[hidx];
    double shortwave = force->shortwave[hidx];
    double Tfoliage = energy->Tfoliage;
    double coverage = snow->coverage;
    double NetLAI = veg_var->NetLAI;
    double NetSAI = veg_var->NetSAI;
    double fcanopy = veg_var->fcanopy;

    /**************************************
      计算直射和漫射辐射通量[0]-vis, [1]-nir.
    **************************************/
    shortwave_dir[0] = shortwave * param.RAD_DIR_F * param.RAD_VIS_F;                 // 直射-可见光
    shortwave_dir[1] = shortwave * param.RAD_DIR_F * (1.0 - param.RAD_VIS_F);         // 直射-近红外
    shortwave_dfs[0] = shortwave * (1.0 - param.RAD_DIR_F) * param.RAD_VIS_F;         // 漫射-可见光
    shortwave_dfs[1] = shortwave * (1.0 - param.RAD_DIR_F) * (1.0 - param.RAD_VIS_F); // 漫射-近红外

    /*******************************
      Advected heat flux from Prec
    *******************************/
    AdvectedEnergy(hidx, force, 
                   energy, veg_var);
                    
    /***************************
      Surface shortwave albedo
    ***************************/
    surface_albedo(step_dt, coszen, 
                   snowfall, energy, 
                   cell, snow, 
                   soil_con, 
                   veg_var, veg_lib);

    /***************************
      Surface radiative fluxes
    ***************************/
    double NetShortSub = 0.0;
    double NetShortSurf = 0.0;
    double NetShortGrnd = 0.0;
    double NetShortSoil = 0.0;
    double NetShortSnow = 0.0;
    for (i = 0; i < options.Nswband; i++) {
        // absorbed by canopy
        ShortOverDir[i] = shortwave_dir[i] * energy->AbsSubDir[i];
        ShortOverDfs[i] = shortwave_dfs[i] * energy->AbsSubDfs[i];
        NetShortSub += ShortOverDir[i] + ShortOverDfs[i];
        // transmitted solar fluxes incident on grnd.
        transmit_direct = shortwave_dir[i] * energy->ShortDir2Dir[i];
        transmit_diffuse = shortwave_dir[i] * energy->ShortDfs2Dir[i] +
                                shortwave_dfs[i] * energy->ShortDfs2Dfs[i];
        // solar radiation absorbed by ground surface.
        tmp_absorb_grnd = transmit_direct * (1.0 - energy->AlbedoGrndDir[i]) + 
                            transmit_diffuse * (1.0 - energy->AlbedoGrndDfs[i]);
        NetShortGrnd += tmp_absorb_grnd;

        // calculate absorbed solar by soil/snow separately
        NetShortSoil += transmit_direct * (1.0 - energy->AlbedoSoilDir[i]) + 
                                        transmit_diffuse * 
                                            (1.0 - energy->AlbedoSoilDfs[i]);
        NetShortSnow += transmit_direct * (1.0 - energy->AlbedoSnowDir[i]) + 
                                        transmit_diffuse * 
                                            (1.0 - energy->AlbedoSnowDfs[i]);
    }
    NetShortSurf = NetShortGrnd + NetShortSub;
    energy->NetShortGrnd = NetShortGrnd;
    energy->shortwave = NetShortSurf;
    energy->NetShortSub = NetShortSub;
    energy->NetShortSoil = NetShortSoil;
    energy->NetShortSnow = NetShortSnow;
   
    /* Compute surface radiative fluxes */
    double f_sun = veg_var->f_sun;
    double f_shade = veg_var->f_shade;
    double leaf_sun = veg_var->leaf_sun;
    double leaf_sha = veg_var->leaf_sha;
    double NetLAI_frac = NetLAI / max(NetLAI + NetSAI, param.TOL_A);
    if (f_sun > 0.0) {
        aPAR_sun = (ShortOverDir[0] + f_sun * ShortOverDfs[0]) *
                        NetLAI_frac / max(leaf_sun, param.TOL_A);
        aPAR_sha = (ShortOverDfs[0] * f_shade) *
                        NetLAI_frac / max(leaf_sha, param.TOL_A);                      
    }
    else {
        aPAR_sun = 0.0;
        aPAR_sha = (ShortOverDir[0] + ShortOverDfs[0]) *
                        NetLAI_frac / max(leaf_sha, param.TOL_A);
    }
    veg_var->aPAR_sun = aPAR_sun;
    veg_var->aPAR_sha = aPAR_sha;
    /* reflected solar radiation */
    double refl_vis = energy->AlbedoSurfDir[0] * shortwave_dir[0] + 
                        energy->AlbedoSurfDfs[0] * shortwave_dfs[0];
    double refl_nir = energy->AlbedoSurfDir[1] * shortwave_dir[1] + 
                        energy->AlbedoSurfDfs[1] * shortwave_dfs[1];
    energy->ReflShortSurf = refl_vis + refl_nir;
    energy->ReflShortGrnd = energy->ReflGrndDir[0] * shortwave_dir[0] + 
                    energy->ReflGrndDfs[0] * shortwave_dfs[0] +
                        energy->ReflGrndDir[1] * shortwave_dir[1] + 
                            energy->ReflGrndDfs[1] * shortwave_dfs[1];
    energy->ReflShortSub = energy->ReflSubDir[0] * shortwave_dir[0] +
                    energy->ReflSubDfs[0] * shortwave_dfs[0] +
                        energy->ReflSubDir[1] * shortwave_dir[1] + 
                            energy->ReflSubDfs[1] * shortwave_dfs[1];

    /******************************
      Compute longwave emissivity
    ******************************/
    double EmissLongSub = 1.0 - exp(-(NetLAI + NetSAI) / 1.0);
    double EmissLongGrnd = param.EMISS_SNOW * coverage +
                           param.EMISS_ICE * (1.0 - coverage);
    double EmissLongSurf = fcanopy * (EmissLongGrnd * (1.0 - EmissLongSub) + EmissLongSub + 
                        EmissLongSub * (1.0 - EmissLongSub) * 
                            (1.0 - EmissLongGrnd)) + (1.0 - fcanopy) * EmissLongGrnd;
    energy->EmissLongSub = EmissLongSub;
    energy->EmissLongGrnd = EmissLongGrnd;
    energy->EmissLongSurf = EmissLongSurf;
    // 更新体积热容量
    for (i = 0; i < cell->Nnode; i++) {
        energy->last_T[i] = energy->T[i];
        energy->last_Cs[i] = energy->Cs_node[i];
    }
    // 更新土层体积水和冰分数
    for (i = 0; i < cell->Nsoil; i++) {
        cell->last_ice[i] = cell->ice[i];
        cell->last_liq[i] = cell->liq[i];
        cell->last_matric[i] = cell->matric[i];
    }
                  
    /**************************************************
       Begin iterative solution of surface fluxes
    **************************************************/
    do {
        // 迭代控制标志
        bool energy_flag = false;
        bool moist_flag = false;
        size_t iter_dt = step_dt;
        // 将当前迭代的结构体值赋给临时结构体，以便在迭代过程中更新
        iter_cell = (*cell);
        iter_energy = (*energy);
        iter_veg_var = (*veg_var);
        iter_snow = (*snow);

        // 尝试找到一个可接受的步长
        while (energy_flag == false || moist_flag == false) {

            /**********************************************
             Solve Energy Balance Components
            **********************************************/
            ErrorFlag = calc_energy_bal(hidx, step_dt, force, 
                                        &iter_energy, &iter_cell,
                                        &iter_snow, soil_con, 
                                        &iter_veg_var, veg_lib);

            if (ErrorFlag == ERROR) {
                return (ERROR);
            }

            // 检查温度是否收敛
            if (iter_energy.delt_T <= temp_error) {
                energy_flag = true;
            } 
            else {
                energy_flag = false;
            }
            
            /*******************************************
             Solve Water Balance Components
            ********************************************/
            ErrorFlag = calc_water_bal(step_dt, pressure,
                                       &iter_energy, 
                                       &iter_cell, soil_con);

            if (ErrorFlag == ERROR) {
                // Return error flag to skip rest of grid cell
                return (ERROR);
            }
            // 检查土壤水是否收敛
            if (iter_energy.delt_Q <= moist_error) {
                moist_flag = true;
            }
            else {
                moist_flag = false;
            }
            iter++;
            // 迭代控制
            if (iter >= param.MAX_ITER_OVER) {
                if (iter_dt > dt_min) {
                    // 时间步长减半，重新开始
                    iter_dt = iter_dt / 2.0;
                    energy_flag = false;
                    moist_flag = false;
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
                    energy_flag = true;
                    moist_flag = true;
                    break;
                }
            }
        }
        time_accum += iter_dt;
    }
    while(time_accum < step_dt);

    // Update the original structures with the 
    // converged values from the last iteration.
    (*cell) = iter_cell;
    (*energy) = iter_energy;
    (*veg_var) = iter_veg_var;
    (*snow) = iter_snow;

    /********************************************************
      Compute Runoff, Baseflow, and Soil Moisture Transport
    ********************************************************/
    snow_hydrology(hidx, step_dt, 
                   force, 
                   energy, cell, 
                   snow, soil_con);
    /* 计算进入土层表面的水量 */
    // pack_melt: 雪层融化的水量[mm]
    // pack_comb: 薄雪层合并中pack_liq导致的表面积水[mm]
    // pack_transp: 从多层积雪转变为无积雪时由薄雪造成的地表积水深度[mm]
    double snow_outflow = 0.0;
    if (snow->Nsnow > 0) {
        size_t Nsnow = snow->Nsnow;
        snow_outflow = snow->pack_outflow[Nsnow-1] / step_dt;
    }
    double soil_inflow = (snow->pack_melt + snow->pack_comb +
                            snow->pack_transp) / step_dt / MM_PER_M; // 转换为m/s
    // 添加雪层多余水分,露水,降水
    soil_inflow += (snow_outflow + cell->dewsoil + 
                    rainfall * (1.0 - snow->coverage)) / MM_PER_M;
    soil_inflow -= cell->esoil;
    
    cell->soil_inflow = soil_inflow;

    ErrorFlag = runoff(step_dt, soil_inflow, cell, soil_con);

    cell->baseflow += snow->glac_excess;

    return (0);
}
