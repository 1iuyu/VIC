/******************************************************************************
* @section DESCRIPTION
*
* This routine was written to handle the various calls and data
* handling needed to solve the various components of the new VIC
* snow code for both the full_energy and water_balance models.
******************************************************************************/

#include "vic_run.h"

/******************************************************************************
* @brief        This routine was written to handle the various calls and data
*               handling needed to solve the various components of the new VIC
*               snow code for both the full_energy and water_balance models.
******************************************************************************/
int
calc_energy_bal(size_t             hidx,
                double             step_dt,
                double             air_temp,
                force_data_struct *force,
                energy_bal_struct *energy,
                cell_data_struct  *cell,
                snow_data_struct  *snow,
                soil_con_struct   *soil_con,
                veg_var_struct    *veg_var,
                veg_lib_struct    *veg_lib)
{
    int    ErrorFlag;
    /***************
     MAIN ROUTINE
    ***************/
    double pressure = force->pressure[hidx];
    double Tfoliage = energy->Tfoliage;
    double fcanopy = veg_var->fcanopy;
    double coverage = snow->coverage;
    /************************************
      Update energy balance variables
    ************************************/
    update_fluxes(energy, cell,
                  snow, 
                  soil_con);

    /*****************************
      Compute surface resistance
    *****************************/
    compute_soil_resis(coverage, cell, soil_con);

    /***************************************
      Compute the vapor flux between nodes
    ***************************************/
    soil_vapor_flux(pressure, cell, soil_con);

    /***************************************
      Compute the hydraulic conductivity
    ***************************************/
    soil_hydraulic_conductivity(cell, soil_con);

    /**********************************
    Compute psychrometric variables
    **********************************/
    if (Tfoliage > CONST_TKFRZ) {
        energy->LatentVapOver = CONST_LATVAP;
        energy->FrozenOver = false;
    }
    else {
        energy->LatentVapOver = CONST_LATSUB;
        energy->FrozenOver = true;
    }

    // 假设只有当liq[0] = 0.0 时发生升华
    if (cell->ice[0] > 0.0 && energy->Tgrnd < CONST_TKFRZ) {
        energy->LatentVapGrnd = CONST_LATSUB;
        energy->FrozenGrnd = true;
    }
    else {
        energy->LatentVapGrnd = CONST_LATVAP;
        energy->FrozenGrnd = false;
    }
    
    /**********************************
      Bare ground surface energy flux
    **********************************/
    ErrorFlag = func_surf_energy_bal(hidx, air_temp, 
                                     force,
                                     energy,
                                     cell, snow, 
                                     soil_con);

    if (ErrorFlag == ERROR) {
        return (ERROR);
    }

    /**********************************
      Vegetation surface energy flux
    **********************************/
    if (cell->IS_VEG == true) {
        ErrorFlag = func_canopy_energy_bal(hidx, step_dt, 
                                           air_temp, force,
                                           energy, cell, 
                                           snow, soil_con,
                                           veg_var, veg_lib);

        if (ErrorFlag == ERROR) {
            return (ERROR);
        }
    }

    /********************************
      Compute grid mean quantities
    ********************************/
    if (cell->IS_VEG) {
        energy->longwave = fcanopy * energy->NetLongSub +
                            (1.0 - fcanopy) * energy->NetLongGrnd;
        energy->sensible = fcanopy * energy->SensibleSub + 
                            (1.0 - fcanopy) * energy->SensibleGrnd;
        energy->latent = fcanopy * energy->LatentSub +
                        (1.0 - fcanopy) * energy->LatentGrnd;
        energy->advection = fcanopy * energy->AdvectSub +
                            (1.0 - fcanopy) * energy->AdvectGrnd + 
                            energy->AdvectOver;
        energy->deriv_terms = fcanopy * energy->deriv_sub + 
                            (1.0 - fcanopy) * energy->deriv_grnd;
        cell->esoil = fcanopy * cell->esoil_sub + 
                            (1.0 - fcanopy) * cell->esoil_grnd;
        energy->deriv_evap = fcanopy * energy->deriv_esub +
                            (1.0 - fcanopy) * energy->deriv_egrnd;
    }
    else {
        energy->longwave = energy->NetLongGrnd;
        energy->sensible = energy->SensibleGrnd;
        energy->latent = energy->LatentGrnd;
        energy->advection = energy->AdvectGrnd;
        energy->Tsurf = energy->Tgrnd;
        energy->SensibleLeaf = 0.0;
        energy->LatentLeaf = 0.0;
        energy->deriv_terms = energy->deriv_grnd;
        cell->esoil = cell->esoil_grnd;
        energy->deriv_evap = energy->deriv_egrnd;
    }

    /************************************
      Compute snow and soil temperature
    ************************************/
    ErrorFlag = SoilTemperature(step_dt, 
                                pressure, cell,
                                energy, 
                                snow, soil_con);

    if (ErrorFlag == ERROR) {
        return (ERROR);
    }

    /************************************
      Compute plant transpiration sink
    ************************************/
    calc_transp_sink(cell, soil_con, veg_var);

    return (0);
}