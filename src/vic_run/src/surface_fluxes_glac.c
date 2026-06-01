/******************************************************************************
* @section DESCRIPTION
*
* This routine computes all surface fluxes, and solves the snow accumulation
* and ablation algorithm. Solutions are for the current snow band and
* vegetation type.
******************************************************************************/

#include "vic_run.h"

/******************************************************************************
* @brief  This is a modified version of surface_fluxes.c specific for glaciers.
******************************************************************************/
int
surface_fluxes_glac(size_t             hidx,
                    double             step_dt,
                    double             air_temp,
                    double             snowfall,
                    double             rainfall,
                    force_data_struct *force,
                    energy_bal_struct *energy,
                    cell_data_struct  *cell,
                    snow_data_struct  *snow,
                    soil_con_struct   *soil_con)
{
    extern parameters_struct param;
    extern option_struct     options;
    cell_data_struct  iter_cell;
    energy_bal_struct iter_energy;
    snow_data_struct  iter_snow;
    int    ErrorFlag;
    size_t i;
    size_t iter = 0;
    size_t dt_min = 1800; // seconds (half hour)
    double time_accum = 0.0;
    double tol_error = 0.1;    // 0.01K
    double f_snowage = 0.0;
    double shortwave_dir[MAX_SWBANDS];
    double shortwave_dfs[MAX_SWBANDS];
    double Tgrnd = energy->Tgrnd;
    double coverage = snow->coverage;
    double coszen = force->coszen[hidx];
    double pressure = force->pressure[hidx];
    double shortwave = force->shortwave[hidx];
    double *AlbedoGrndDir = energy->AlbedoGrndDir;
    double *AlbedoGrndDfs = energy->AlbedoGrndDfs;

    /***************
      MAIN ROUTINE
    ***************/                             
    shortwave_dir[0] = shortwave * param.RAD_DIR_F * param.RAD_VIS_F;                 // 直射-可见光
    shortwave_dir[1] = shortwave * param.RAD_DIR_F * (1.0 - param.RAD_VIS_F);         // 直射-近红外
    shortwave_dfs[0] = shortwave * (1.0 - param.RAD_DIR_F) * param.RAD_VIS_F;         // 漫射-可见光
    shortwave_dfs[1] = shortwave * (1.0 - param.RAD_DIR_F) * (1.0 - param.RAD_VIS_F); // 漫射-近红外

    /*******************************
      Advected heat flux from Prec
    *******************************/
    AdvectedEnergyGlac(air_temp, snowfall, 
                       rainfall, energy);
    
    /*******************************
    Thermal properties for the layer
    *******************************/
	prepare_full_energy(cell->IS_GLAC,
                        cell, energy,
						snow, soil_con);

    /***************************
      Surface shortwave albedo
    ***************************/
    f_snowage = snow_aging(step_dt, Tgrnd,
                           snowfall, snow);

    /** compute understory albedo and net shortwave radiation **/
    if (coszen > 0.) {
        // age snow albedo if no new snowfall
        // solar radiation process is only done if there is light
        snow_albedo(coszen, f_snowage, energy);

        for (i = 0; i < options.Nswband; i++) {
            AlbedoGrndDir[i] =
                (param.GLAC_ALBEDO[i] * (1. - coverage) + 
                                energy->AlbedoSnowDir[i] * coverage);
            AlbedoGrndDfs[i] =
                (param.GLAC_ALBEDO[i] * (1. - coverage) + 
                                energy->AlbedoSnowDfs[i] * coverage);
        }
    }
    else {
        for (i = 0; i < options.Nswband; i++) {
            AlbedoGrndDir[i] = 0.0;
            AlbedoGrndDfs[i] = 0.0;
            energy->AlbedoSnowDir[i] = 0.0;
            energy->AlbedoSnowDfs[i] = 0.0;
        }
    }

    /***************************
      Surface radiative fluxes
    ***************************/
    double NetShortGrnd = 0.0;
    double ReflShortSurf = 0.0;
    for (i = 0; i < options.Nswband; i++) {
        NetShortGrnd += 
            shortwave_dir[i] * (1.0 - AlbedoGrndDir[i]) +
                shortwave_dfs[i] * (1.0 - AlbedoGrndDfs[i]);
        ReflShortSurf +=
            shortwave_dir[i] * AlbedoGrndDir[i] +
                shortwave_dfs[i] * AlbedoGrndDfs[i];
        energy->NetShortGrnd = NetShortGrnd;
        energy->ReflShortSurf = ReflShortSurf;
        energy->NetShortSurf = NetShortGrnd;
    }

    /******************************
      Compute longwave emissivity
    ******************************/
    double EmissLongGrnd = param.EMISS_ICE * (1.0 - coverage) + 
                                    param.EMISS_SNOW * coverage;
    energy->EmissLongGrnd = EmissLongGrnd;
    double Ra_evap = 1.0;
    double rh_grnd = 1.0;
    cell->Ra_evap = Ra_evap;
    cell->rh_grnd = rh_grnd;
    
    /******************************
      Set psychrometric variables
    ******************************/
    energy->LatentVapGrnd = CONST_LATSUB;

    /*********************************
      Compute glacier surface fluxes
    *********************************/
    do {
        // 迭代控制标志
        bool iter_flag = false;
        double iter_dt = step_dt;

        // 将当前迭代的结构体值赋给临时结构体，以便在迭代过程中更新
        iter_cell = (*cell);
        iter_energy = (*energy);
        iter_snow = (*snow);

        // 尝试找到一个可接受的步长
        while (iter_flag == false) {        
            /** Solve energy balence processes **/
            ErrorFlag = calc_energy_bal_glac(hidx, iter_dt,
                                             air_temp, force,
                                            &iter_energy, 
                                            &iter_cell,
                                            &iter_snow, soil_con);

            if (ErrorFlag == ERROR) {
                return (ERROR);
            }

            // 检查能量是否收敛
            if (iter_energy.delt_T <= tol_error) {
                iter_flag = true;
            } 
            else {
                iter_flag = false;
            }
            
            // 检查温度变化幅度是否过大
            if (fabs(iter_energy.delt_T) > 25.0) {
                // 时间步长减半，重新开始
                iter_dt = iter_dt / 2.0;
                iter_flag = false;
                iter = 0;
                // 回滚状态至初始值
                iter_cell = (*cell);
                iter_energy = (*energy);
                iter_snow = (*snow);
                continue;
            }

            iter++;
            // 迭代控制
            if (iter >= param.MAX_ITER_OVER) {
                if (iter_dt > dt_min) {
                    // 时间步长减半，重新开始
                    iter_dt = iter_dt / 2.0;
                    iter_flag = false;
                    iter = 0;
                    // 回滚状态至初始值
                    iter_cell = (*cell);
                    iter_energy = (*energy);
                    iter_snow = (*snow);
                    continue;
                } 
                else {
                    // 已达到最小时间步长，强制收敛并继续
                    iter_flag = true;
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
    (*snow) = iter_snow;

    /** Solve water balence processes **/
    ErrorFlag = calc_water_bal_glac(step_dt, air_temp,
                                    snowfall, rainfall,
                                    pressure,
                                    energy, cell, 
                                    snow, soil_con);
    if (ErrorFlag == ERROR) {

        return (ERROR);
    }

    return(0);
}
