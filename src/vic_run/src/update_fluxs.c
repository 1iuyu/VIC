/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes the fluxes for each cell in the grid. It sets the
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    This routine initializes roughness lengths and 
 *           forcing heights for each cell in the grid.
 *****************************************************************************/
int
update_fluxes(double             step_dt,
              force_data_struct *force,
              energy_bal_struct *energy,
              cell_data_struct  *cell,
              snow_data_struct  *snow,
              soil_con_struct   *soil_con,
              veg_var_struct    *veg_var,
              veg_lib_struct    *veg_lib)
{
    extern option_struct     options;
    int    ErrorFlag;
    size_t Nsnow = snow->Nsnow;

    /********************************
       Updata surface temperature
    ********************************/
    double coverage = snow->coverage;
    double pressure = force->pressure;
    double Tfoliage = energy->Tfoliage;
    double frac_h2o = cell->frac_h2o;
    double *soil_T = cell->soil_T;
    double *pack_T = snow->pack_T;
    double h2osfc_T = cell->h2osfc_T;
    if (Nsnow > 0) {
        if (frac_h2o > 0.0) {
            energy->Tgrnd = (1.0 - coverage) * h2osfc_T + coverage * pack_T[0];
        } else {
            energy->Tgrnd = (1.0 - coverage) * soil_T[0] + coverage * pack_T[0];
        }
    }
    else {
        if (frac_h2o > 0.0) {
            energy->Tgrnd = (1.0 - frac_h2o) * soil_T[0] + frac_h2o * h2osfc_T;
        } else {
            energy->Tgrnd = soil_T[0];
        }
    }

    /********************************
      Calculate active layer depth
    ********************************/
    if (options.ACTIVE_LAYER) {
        ActiveLayer(energy, cell, soil_con);
    }

    /**********************************
      Compute psychrometric variables
    **********************************/
    if (cell->IS_VEG) {
        if (Tfoliage > CONST_TKFRZ) {
            energy->LatentVapOver = CONST_LATVAP;
            energy->FrozenOver = false;
        }
        else {
            energy->LatentVapOver = CONST_LATSUB;
            energy->FrozenOver = true;
        }
    }
    // 假设只有当liq[0] = 0.0 时发生升华
    if ((energy->ice[Nsnow] > 0.0 && energy->Tgrnd < CONST_TKFRZ) || cell->IS_GLAC) {
        energy->LatentVapGrnd = CONST_LATSUB;
        energy->FrozenGrnd = true;
    }
    else {
        energy->LatentVapGrnd = CONST_LATVAP;
        energy->FrozenGrnd = false;
    }
    /*****************************
      Compute surface resistance
    *****************************/
    compute_soil_resis(coverage, energy, cell, soil_con);

    /***************************************
      Compute the vapor flux between nodes
    ***************************************/
    calc_vapor_flux(pressure, energy, cell, 
                    snow, soil_con);

    /***************************************
      Compute the hydraulic conductivity
    ***************************************/
    soil_hydraulic_conductivity(energy, cell, snow, soil_con);

    /***************************************
      Compute layer thermal properties
    ***************************************/
    prepare_full_energy(pressure, 
                        cell, energy,
                        snow, soil_con);
    
    /**********************************
      Bare ground surface energy flux
    **********************************/
    ErrorFlag = func_surf_energy_bal(force,
                                     energy,
                                     cell, snow);

    if (ErrorFlag == ERROR) {
        return (ERROR);
    }

    /**********************************
      Vegetation surface energy flux
    **********************************/
    if (cell->IS_VEG || cell->IS_URBAN) {
        ErrorFlag = func_canopy_energy_bal(step_dt, force,
                                           energy, cell, 
                                           snow, soil_con,
                                           veg_var, veg_lib);

        if (ErrorFlag == ERROR) {
            return (ERROR);
        }
    }

    /************************************
      Compute plant transpiration sink
    ************************************/
    if (cell->IS_VEG || cell->IS_URBAN) {
        calc_transp_sink(cell, soil_con, veg_var);
    }

    /********************************************************
      Compute Runoff, Baseflow, and Soil Moisture Transport
    ********************************************************/
    snow_hydrology(step_dt, force,
                   energy, cell, 
                   snow, soil_con);
    /* 计算进入土层表面的水量 */
    // pack_comb: 薄雪层合并中pack_liq导致的表面积水[mm]
    // snow_outflow：流出雪层的水[mm/s]
    // soil_inflow: 裸土区流入的水[mm/s]
    double soil_inflow = (snow->pack_comb / step_dt + 
                          snow->snow_outflow + cell->soil_inflow) / MM_PER_M; // 转换为m/s
    cell->soil_inflow = soil_inflow;

    ErrorFlag = runoff(step_dt, soil_inflow, energy, cell, soil_con);

    cell->baseflow += snow->glac_excess;
    
    return (0);
}