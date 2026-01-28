/******************************************************************************
* @section DESCRIPTION
*
* This routine was written to handle the various calls and data
* handling needed to solve the various components of the new VIC
* snow code for both the full_energy and water_balance models.
******************************************************************************/

#include <vic_run.h>

/******************************************************************************
* @brief        This routine was written to handle the various calls and data
*               handling needed to solve the various components of the new VIC
*               snow code for both the full_energy and water_balance models.
******************************************************************************/
int
calc_energy_bal(size_t             hidx,
                double             step_dt,
                double             air_temp,
                double             rainfall,
                double             snowfall,
                double             ref_height,
                double            *roughness,
                double            *displacement,
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

    int             ErrorFlag;
    size_t          i;
    double          proj_area;
    double          APAR_sunlit;
    double          APAR_shade;
    double          wet_fact;
    double          f_sun;
    double          transmit_direct;
    double          transmit_diffuse;
    double          tmp_absorb_grnd;
    double          psi_soil[MAX_NODES];
    double          ShortOverDir[MAX_SWBANDS];
    double          ShortOverDfs[MAX_SWBANDS];
    double          shortwave_dir[MAX_SWBANDS];
    double          shortwave_dfs[MAX_SWBANDS];

    /**************************
       Prepare node variables
    **************************/    
    double *liq = cell->liq;
    double *dz_soil = soil_con->dz_soil;
    double *Zsum_soil = soil_con->Zsum_soil;
    double *Wfc_node = soil_con->Wfc_node;
    double *Wsat_node = soil_con->Wsat_node;
    double *Wpwp_node = soil_con->Wpwp_node;
    double *bexp_node = soil_con->bexp_node;
    double *transp_fact = cell->transp_fact;
    double *psi_sat = soil_con->psi_sat_node;

    /***************
      MAIN ROUTINE
    ***************/
    double vp = force->vp[hidx];
    double Qair = force->Qair[hidx];
    double wind = force->wind[hidx];
    double coszen = force->coszen[hidx];
    double pressure = force->pressure[hidx];
    double longwave = force->longwave[hidx];
    double shortwave = force->shortwave[hidx];
    double air_density = force->density[hidx];
    double Tgrnd = energy->Tgrnd;
    double Tcanopy = energy->Tcanopy;
    double coverage = snow->coverage;
    double NetLAI = veg_var->NetLAI;
    double NetSAI = veg_var->NetSAI;
    double wetFrac = veg_var->wetFrac;
    double fcanopy = veg_var->fcanopy;
    double theta = air_temp * pow((pressure / pressure), CONST_RDAIR / CONST_CPDAIR);

    /* initialize diffuse and direct partitioned
       from shortwave. [0]-vis, [1]-nir.*/
    shortwave_dir[0] = shortwave * param.RAD_DIR_F *
                                                param.RAD_VIS_F;   // 直射-可见光
    shortwave_dir[1] = shortwave * param.RAD_DIR_F * 
                                          (1 - param.RAD_VIS_F);   // 直射-近红外
    shortwave_dfs[0] = shortwave * (1 - param.RAD_DIR_F) *
                                                param.RAD_VIS_F;   // 漫射-可见光
    shortwave_dfs[1] = shortwave * (1 - param.RAD_DIR_F) * 
                                          (1 - param.RAD_VIS_F);   // 漫射-近红外
    /*******************************
      Advected heat flux from Prec
    *******************************/
    AdvectedEnergy(fcanopy, Tcanopy,
                   Tgrnd, air_temp,
                   rainfall, snowfall,
                   veg_var->RainThroughFall,
                   veg_var->SnowThroughFall,
                   veg_var->RainDrip,
                   veg_var->SnowDrip,
                   &energy->AdvectOver,
                   &energy->AdvectGrnd,
                   &energy->AdvectSub);

    /***************************
      Surface shortwave albedo
    ***************************/

    /* snow grain size and aging for SNICAR */
    double SnowAge_fact = snow_aging(step_dt, snow->swq,
                                     &snow->old_swq,
                                     Tgrnd, snowfall,
                                     &snow->SnowAge);

    /** compute understory albedo and net shortwave radiation **/
    if (coszen > 0.) {

        double leaf_frac = NetLAI / max(NetLAI + NetSAI, param.TOL_A);
        double stem_frac = NetSAI / max(NetLAI + NetSAI, param.TOL_A);

        for (i = 0; i < options.Nswband; i++) {
            energy->ReflectVeg[i] = max(veg_lib->reflleaf[i] * leaf_frac + 
                                        veg_lib->reflstem[i] * 
                                            stem_frac, param.TOL_A);
            energy->TransmitVeg[i] = max(veg_lib->transleaf[i] * leaf_frac +
                                        veg_lib->transstem[i] * 
                                            stem_frac, param.TOL_A);
        }
        // age snow albedo if no new snowfall
        // solar radiation process is only done if there is light
        snow_albedo(coszen, SnowAge_fact,
                    energy->AlbedoSnowDir,
                    energy->AlbedoSnowDfs);

        /* Compute ground albedo based on soil and snow albedo */
        GroundAlbedo(cell->moist[0],
                     coverage,
                     soil_con->AlbedoSat,
                     soil_con->AlbedoDry,
                     energy->AlbedoSoilDir,
                     energy->AlbedoSoilDfs,
                     energy->AlbedoSnowDir,
                     energy->AlbedoSnowDfs,
                     energy->AlbedoGrndDir,
                     energy->AlbedoGrndDfs);

        /* Compute canopy radiative transfer 
           using two-stream approximation */
        for (i = 0; i < options.Nswband; i++) {
            // direct
            canopy_two_stream(i, 0, Tcanopy,
                              NetLAI, NetSAI,
                              fcanopy,
                              coszen, &proj_area,
                              wetFrac, 
                              energy, veg_lib);
            // diffuse
            canopy_two_stream(i, 1, Tcanopy,
                              NetLAI, NetSAI,
                              fcanopy,
                              coszen, &proj_area, 
                              wetFrac, 
                              energy, veg_lib);
        }
        double tau_beam = proj_area / coszen *
                        sqrt(1.0 - energy->ReflectVeg[0] -
                                            energy->TransmitVeg[0]);
        f_sun = (1.0 - exp(-tau_beam *
                            (NetLAI + NetSAI))) /
                                        max(tau_beam * 
                                                (NetLAI + NetSAI), 
                                                        param.TOL_A);
        tau_beam = f_sun;
        if (tau_beam < 0.01) {
            leaf_frac = 0.0;
        }
        else {
            leaf_frac = tau_beam;
        }
        f_sun = leaf_frac;
    } // end of coszen > 0.
    else {
        f_sun = 0.0;

    }

    /* shaded canopy fraction */
    double f_shade = 1.0 - f_sun;
    double leaf_sun = NetLAI * f_sun;
    double leaf_shade = NetLAI * f_shade;
    veg_var->leaf_sun = leaf_sun;
    veg_var->leaf_shade = leaf_shade;

    /***************************
      Surface radiative fluxes
    ***************************/
    double NetShortSub = 0.;
    double NetShortSurf = 0.;
    
    for (i = 0; i < options.Nswband; i++) {
        // absorbed by canopy
        ShortOverDir[i] = shortwave_dir[i] * energy->AbsSubDir[i];
        ShortOverDfs[i] = shortwave_dfs[i] * energy->AbsSubDfs[i];
        NetShortSub += ShortOverDir[i] + ShortOverDfs[i];
        NetShortSurf += ShortOverDir[i] + ShortOverDfs[i];
        // transmitted solar fluxes incident on grnd.
        transmit_direct = shortwave_dir[i] * energy->ShortDir2Dir[i];
        transmit_diffuse = shortwave_dir[i] * energy->ShortDfs2Dir[i] +
                                shortwave_dfs[i] * energy->ShortDfs2Dfs[i];
        // solar radiation absorbed by ground surface.
        tmp_absorb_grnd = transmit_direct * (1.0 - energy->AlbedoGrndDir[i]) + 
                            transmit_diffuse * (1.0 - energy->AlbedoGrndDfs[i]);
        NetShortSurf += tmp_absorb_grnd;
        // calculate absorbed solar by soil/snow separately
        energy->NetShortSoil += transmit_direct * 
                                    (1.0 - energy->AlbedoSoilDir[i]) + 
                                        transmit_diffuse * 
                                            (1.0 - energy->AlbedoSoilDfs[i]);
        energy->NetShortSnow += transmit_direct * 
                                    (1.0 - energy->AlbedoSnowDir[i]) + 
                                        transmit_diffuse * 
                                            (1.0 - energy->AlbedoSnowDfs[i]);
        energy->NetShortGrnd += tmp_absorb_grnd;
        energy->shortwave = NetShortSurf;
        energy->NetShortSub = NetShortSub;
        if (snow->Nsnow == 0) {
            energy->NetShortSoil = NetShortSurf;
            energy->NetShortSnow = NetShortSurf;
        }
    }

    /* Compute surface radiative fluxes */
    double NetLAI_frac = NetLAI / max(NetLAI + NetSAI, param.TOL_A);
    if (f_sun > 0.0) {
        APAR_sunlit = (ShortOverDir[0] + f_sun * ShortOverDfs[0]) *
                        NetLAI_frac / max(leaf_sun, param.TOL_A);
        APAR_shade = (ShortOverDfs[0] * f_shade) *
                        NetLAI_frac / max(leaf_shade, param.TOL_A);                      
    }
    else {
        APAR_sunlit = 0.0;
        APAR_shade = (ShortOverDir[0] + ShortOverDfs[0]) *
                        NetLAI_frac / max(leaf_shade, param.TOL_A);
    }
    energy->APAR_sunlit = APAR_sunlit;
    energy->APAR_shade = APAR_shade;
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
    double EmissLongGrnd = param.EMISS_SNOW * coverage + param.EMISS_ICE * (1.0 - coverage);
    double EmissLongSurf = fcanopy * (EmissLongGrnd * (1.0 - EmissLongSub) + EmissLongSub + 
                        EmissLongSub * (1.0 - EmissLongSub) * 
                            (1.0 - EmissLongGrnd)) + (1.0 - fcanopy) * EmissLongGrnd;
    energy->EmissLongSub = EmissLongSub;
    energy->EmissLongGrnd = EmissLongGrnd;
    energy->EmissLongSurf = EmissLongSurf;

    /************************************
      Compute soil transpiration factor
    ************************************/
    double psi_wp = param.SRESP_PSIWILT;
    double total_transp = 0.;
    for (i = 0; i < options.Nroot; i++) {
        if (options.SOIL_TRANSP == NOAH) {
            wet_fact = (liq[i] - Wpwp_node[i]) / 
                        (Wfc_node[i] - Wpwp_node[i]);
        }
        else if (options.SOIL_TRANSP == CLM) {
            psi_soil[i] = max(psi_wp, -psi_sat[i] * 
                                pow(max(0.01, liq[i]) /
                                        Wsat_node[i], -bexp_node[i]));
            wet_fact = (1.0 - psi_soil[i] / psi_wp) / 
                            (1.0 + psi_sat[i] / psi_wp);
        }
        else {
            log_err("Unknown SOIL_TRANSP option");
        }
        wet_fact = min(1.0, max(0.0, wet_fact));

        transp_fact[i] = max(param.TOL_A, dz_soil[i] / Zsum_soil[options.Nroot] * wet_fact);

        total_transp += transp_fact[i];
    }
    total_transp = min(max(total_transp, param.TOL_A), 1.0);
    cell->total_transp = total_transp;
    for (i = 0; i < options.Nroot; i++) {
        transp_fact[i] /= total_transp;
    }

    /*****************************
      Compute surface resistance
    *****************************/
    double esoil_frac = max(0.0, liq[0] / Wsat_node[0]);
    double dry_depth = dz_soil[0] * (exp(pow(1.0 - min(1.0, liq[0] / 
                            Wsat_node[0]), param.SRESP_EXP))
                                             - 1.0) / (2.71828 - 1.0);
    double Dvap_red = 2.2e-5 * Wsat_node[0] * Wsat_node[0] *
                        pow(1.0 - (Wpwp_node[0] / Wsat_node[0]), 
                                      2.0 + 3.0 / bexp_node[0]);
    double Ra_evap = dry_depth / Dvap_red;

    if (liq[0] < 0.01 && snow->snow_depth == 0.0) {
        Ra_evap = param.MAX_LIMIT;
    }
//    Ra_evap = 1.0 / (coverage * (1.0 / param.SNOW_RASNOW) + 
//                    (1.0 - coverage) * (1.0 / max(Ra_evap, 0.001)));

    cell->Ra_evap = Ra_evap;
    double psi_surf = -psi_sat[0] *
                    pow((max(0.01, liq[0]) /
                            Wsat_node[0]), -bexp_node[0]);
    double rh_grnd = coverage + (1.0 - coverage) *
                        exp(psi_surf * CONST_G /
                                (CONST_RWV * Tgrnd));

    /**********************************
      Compute psychrometric variables
    **********************************/
    if (Tcanopy > CONST_TKFRZ) {
        energy->LatentVapOver = CONST_LATVAP;
        energy->FrozenOver = false;
    }
    else {
        energy->LatentVapOver = CONST_LATSUB;
        energy->FrozenOver = true;
    }
    /* Calculate psychrometric constant */
    energy->PsyCh_canopy = CONST_CPDAIR * pressure / 
                        (0.622 * energy->LatentVapOver);

    if (Tgrnd > CONST_TKFRZ) {
        energy->LatentVapGrnd = CONST_LATVAP;
        energy->FrozenGrnd = false;
    }
    else {
        energy->LatentVapGrnd = CONST_LATSUB;
        energy->FrozenGrnd = true;
    }
    /* Calculate psychrometric constant */
    energy->PsyCh_grnd = CONST_CPDAIR * pressure /
                        (0.622 * energy->LatentVapGrnd);

    /*********************************
      Vegetation surface energy flux
    **********************************/
    if (NetLAI + NetSAI > 0. && fcanopy > 0.) {
        energy->Tupper = Tgrnd;
        ErrorFlag = func_canopy_energy_bal(step_dt, wind,
                                           Qair, vp, longwave,
                                           air_temp,
                                           air_density,
                                           pressure,
                                           rh_grnd, ref_height,
                                           roughness, displacement,
                                           energy, cell,
                                           snow, soil_con,
                                           veg_var, veg_lib);

        if (ErrorFlag == ERROR) {
            return (ERROR);
        }
    }
    /**********************************
      Bare ground surface energy flux
    **********************************/
    energy->Tlower = Tgrnd;
    ErrorFlag = func_surf_energy_bal(air_density,
                                     Ra_evap, longwave,
                                     air_temp,
                                     rh_grnd, theta,
                                     vp, pressure,
                                     Qair, wind, ref_height,
                                     roughness, displacement, 
                                     energy, cell, 
                                     snow, soil_con);

    if (ErrorFlag == ERROR) {
        return (ERROR);
    }
    /********************************
      Compute grid mean quantities
    ********************************/
    if (NetLAI + NetSAI > 0.0 && fcanopy > 0.0) {
        energy->longwave = fcanopy * energy->NetLongSub + 
                    (1.0 - fcanopy) * energy->NetLongGrnd + 
                                            energy->NetLongOver;
        energy->sensible = fcanopy * energy->SensibleSub + 
                    (1.0 - fcanopy) * energy->SensibleGrnd + 
                                            energy->SensibleOver;
        energy->latent = fcanopy * energy->LatentSub +
                    (1.0 - fcanopy) * energy->LatentGrnd;
        energy->grnd_flux = fcanopy * energy->GroundSub + 
                    (1.0 - fcanopy) * energy->GroundGrnd;
        energy->advection = fcanopy * energy->AdvectSub +
                    (1.0 - fcanopy) * energy->AdvectGrnd + 
                                            energy->AdvectOver;
        energy->Tgrnd = fcanopy * energy->Tupper + 
                              (1.0 - fcanopy) * energy->Tlower;
        energy->Tsurf = fcanopy * energy->Tcanopy + 
                              (1.0 - fcanopy) * energy->Tlower;
    }
    else {
        energy->longwave = energy->NetLongGrnd;
        energy->sensible = energy->SensibleGrnd;
        energy->latent = energy->LatentGrnd;
        energy->grnd_flux = energy->GroundGrnd;
        energy->advection = energy->AdvectGrnd;
        energy->Tgrnd = energy->Tlower;
        energy->Tsurf = energy->Tgrnd;
        energy->LatentCanopy = 0.0;
        energy->LatentTransp = 0.0;
    }
    
    // emitted longwave radiation and physical check
    double lw_emit = energy->longwave + longwave;

    if (lw_emit <= 0.0) {
        log_err("Emitted longwave <= 0 (value: %.4f). Components: lw_out = %.4f, lw_in = %.4f",
            lw_emit, energy->longwave, longwave);
    }

    /************************************
      Compute snow and soil temperature
    ************************************/
    ErrorFlag = SoilTemperature(step_dt, cell, energy,
                                snow, soil_con);

    if (ErrorFlag == ERROR) {
        return (ERROR);
    }
    /**************************************
      Phase change of snow and soil water
    **************************************/
    CalcPhaseChange(step_dt, energy,
                    cell, snow, soil_con);

    return (0);

}
