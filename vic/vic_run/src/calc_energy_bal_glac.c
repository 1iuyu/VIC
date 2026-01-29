/******************************************************************************
* @section DESCRIPTION
*
* This routine was written to handle the various calls and data
* handling needed to solve the various components of the new VIC
* snow code for both the full_energy and water_balance models.
******************************************************************************/

#include <vic_run.h>

/******************************************************************************
* @brief        This function is specific for glaciers and based off of the
*               solve_snow function.
******************************************************************************/
int
calc_energy_bal_glac(size_t             hidx,
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
                     soil_con_struct   *soil_con)
{
    extern option_struct     options;
    extern parameters_struct param;

    int          ErrorFlag;
    size_t       i, lidx;
    double       esat_Tgrnd;
    double       esat_slope;
    double       SVP_liq;
    double       SVP_ice;
    double       liq_slope;
    double       ice_slope; 
    double       Qair_surf;
    double       moisture_flux;    
    double       NetLongGrnd;
    double       SensibleGrnd;
    double       layer_T;
    double       LatentGrnd;
    double       GroundGrnd;
    double       coef_sensible;
    double       coef_latent;
    double       coef_longwave;
    double       coef_ground;
    double       layer_depth;
    double       shortwave_dir[MAX_SWBANDS];
    double       shortwave_dfs[MAX_SWBANDS];

    /**************************
       Prepare node variables
    **************************/  
    double vp = force->vp[hidx];
    double Qair = force->Qair[hidx];
    double Tgrnd = energy->Tgrnd;
    double wind = force->wind[hidx];
    double coverage = snow->coverage;
    double coszen = force->coszen[hidx];
    double pressure = force->pressure[hidx];
    double longwave = force->longwave[hidx];
    double air_density = force->density[hidx];
    double shortwave = force->shortwave[hidx];
    double *ice = cell->ice;
    double *liq = cell->liq;
    double *moist = cell->moist;
    double *T = energy->T;
    double *Ra_grnd = cell->Ra_grnd;
    double *kappa_node = energy->kappa_node;
    double *AlbedoGrndDir = energy->AlbedoGrndDir;
    double *AlbedoGrndDfs = energy->AlbedoGrndDfs;
    
    /***************
      MAIN ROUTINE
    ***************/ 
    /* initialize diffuse and direct partitioned 
       from shortwave. [0]-vis, [1]-nir. */
    shortwave_dir[0] = shortwave * param.RAD_DIR_F *
                                                    param.RAD_VIS_F;
    shortwave_dir[1] = shortwave * param.RAD_DIR_F * 
                                              (1 - param.RAD_VIS_F);
    shortwave_dfs[0] = shortwave * (1 - param.RAD_DIR_F) * 
                                                    param.RAD_VIS_F;
    shortwave_dfs[1] = shortwave * (1 - param.RAD_DIR_F) * 
                                              (1 - param.RAD_VIS_F);

    /*******************************
      Advected heat flux from Prec
    *******************************/
    AdvectedEnergyGlac(Tgrnd, air_temp,
                       rainfall, snowfall, 
                      &energy->AdvectGrnd);

    /***************************
      Surface shortwave albedo
    ***************************/
    double SnowAge_fact = snow_aging(step_dt, snow->swq,
                                     &snow->old_swq, 
                                     Tgrnd, snowfall, 
                                     &snow->SnowAge);

    /** compute understory albedo and net shortwave radiation **/
    if (coszen > 0.) {
        // age snow albedo if no new snowfall
        // solar radiation process is only done if there is light
        snow_albedo(coszen, SnowAge_fact,
                    energy->AlbedoSnowDir,
                    energy->AlbedoSnowDfs);

        for (i = 0; i < options.Nswband; i++) {
            AlbedoGrndDir[i] =
                (param.GLAC_ALBEDO[i] * (1. - coverage) + 
                                energy->AlbedoSnowDir[i] * coverage);
            AlbedoGrndDfs[i] =
                (param.GLAC_ALBEDO[i] * (1. - coverage) + 
                                energy->AlbedoSnowDfs[i] * coverage);
        }
    }

    /***************************
      Surface radiative fluxes
    ***************************/
    double NetShortGrnd = 0.;
    double ReflShortSurf = 0.;
    for (i = 0; i < options.Nswband; i++) {
        NetShortGrnd += 
            shortwave_dir[i] * (1. - AlbedoGrndDir[i]) +
                shortwave_dfs[i] * (1. - AlbedoGrndDfs[i]);
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
    double EmissLongGrnd = param.EMISS_ICE * (1. - coverage) + 
                                    param.EMISS_SNOW * coverage;
    energy->EmissLongGrnd = EmissLongGrnd;
    double Ra_evap = 1.0;
    double rh_grnd = 1.0;
    cell->Ra_evap = Ra_evap;
    /******************************
      Set psychrometric variables
    ******************************/
    energy->LatentVapGrnd = CONST_LATSUB;
    double PsyCh_grnd = CONST_CPDAIR * pressure /
                                (0.622 * energy->LatentVapGrnd);
    energy->PsyCh_grnd = PsyCh_grnd;

    /******************************
      Glacier surface energy flux
    ******************************/
    size_t iter = 0;
    double ustar = 0.06;
    double zL_glac = 1.0;
    double delta_grnd = param.TOL_A + 1.;
    double sensible_glac = 0.;
    double AdvectGrnd = energy->AdvectGrnd;
    if (snow->Nsnow > 0) {
        lidx = snow->Nsnow - 1; // 存在雪层
        layer_depth = snow->dz_snow[lidx];
        layer_T = snow->pack_T[0];
    }
    else {
        layer_depth = soil_con->dz_soil[0]; // 不存在雪层，kappa_node第一层为土壤热力学参数
        layer_T = cell->soil_T[0];
    }
    coef_longwave = EmissLongGrnd * CONST_BOLTZ;
    coef_ground = 2.0 * kappa_node[0] / layer_depth;
    /* begin stability iteration for
       ground temperature and flux */
    do {

        ErrorFlag = CalcAerodynamic(iter, Qair,
                                    air_temp, &zL_glac,
                                    &ustar, Ra_grnd,
                                    &sensible_glac, wind,
                                    air_density,
                                    displacement[2],
                                    ref_height,
                                    roughness[2]);
        if (ErrorFlag == ERROR) {
            return (ERROR);
        }
        
        /* conductance variables for diagnostics */
        // double tmp_CM = 1.0 / Ra_glac[0];
        // double tmp_CH = 1.0 / Ra_glac[1];

        /* Saturated vapor pressure at Tgrnd */
        svp(Tgrnd, &SVP_liq, &SVP_ice);
        svp_slope(Tgrnd, &liq_slope, &ice_slope);
        if (Tgrnd > CONST_TKFRZ) {
            esat_Tgrnd = SVP_liq;
            esat_slope = liq_slope;
        }
        else {
            esat_Tgrnd = SVP_ice;
            esat_slope = ice_slope;
        }
        
        /* Calculate the sensible heat flux */
        coef_sensible = CONST_CPDAIR * air_density / Ra_grnd[1];
        if (snow->snow_depth > 0.) {
            coef_latent = air_density * CONST_CPDAIR / 
                            PsyCh_grnd / (Ra_grnd[2] + Ra_evap);
        }
        else {
            coef_latent = 0.0;
        }

        NetLongGrnd = calc_longwave(Tgrnd, coef_longwave) - EmissLongGrnd * longwave;
        SensibleGrnd = calc_sensible_heat(Tgrnd, air_temp, coef_sensible);
        LatentGrnd = calc_latent_heat(esat_Tgrnd, rh_grnd, vp, coef_latent);
        GroundGrnd = calc_ground_heat(Tgrnd, layer_T, coef_ground);

        double RestTerm = NetShortGrnd - NetLongGrnd - SensibleGrnd -
                         LatentGrnd - GroundGrnd + AdvectGrnd;
        double coef_flux = 4.0 * coef_longwave * pow(Tgrnd, 3) + coef_sensible + 
                                         coef_latent * esat_slope + coef_ground;
        delta_grnd = RestTerm / coef_flux;

        NetLongGrnd += 4.0 * coef_longwave * pow(Tgrnd, 3) * delta_grnd;
        SensibleGrnd += coef_sensible * delta_grnd;
        LatentGrnd += coef_latent * esat_slope * delta_grnd;
        GroundGrnd += coef_ground * delta_grnd;

        /* update ground temperature */
        Tgrnd += delta_grnd;

        /* for computing M-O length */
        sensible_glac = coef_sensible * (Tgrnd - air_temp);

        /* update specific humidity */
        svp(Tgrnd, &SVP_liq, &SVP_ice);

        if (Tgrnd > CONST_TKFRZ) {
            esat_Tgrnd = SVP_liq;
        }
        else {
            esat_Tgrnd = SVP_ice;
        }
        Qair_surf = 0.622 * (esat_Tgrnd * rh_grnd) / 
                            (pressure - 0.378 * (esat_Tgrnd * rh_grnd));
        moisture_flux = (Qair_surf - Qair) * coef_latent *
                                         PsyCh_grnd / CONST_CPDAIR;
        iter++;
    } 
    while (iter < param.MAX_ITER_MOST && fabs(delta_grnd) > 0.01);

    /* if snow on ground and Tgrnd > freezing point: reset 
        Tgrnd = freezing point. reevaluate ground fluxes. */
    double tmp_ice = ice[0];
    for (i = 0; i < options.Nnode; i++) {
        ice[i] = moist[i] - liq[i];
        if (tmp_ice < ice[i]) {
            tmp_ice = ice[i];
        } 
    }
    if ((tmp_ice > 0. || snow->snow_depth > 0.05) && Tgrnd > CONST_TKFRZ) {
        Tgrnd = CONST_TKFRZ;
        /* update specific humidity */
        svp(Tgrnd, &SVP_liq, &SVP_ice);
        esat_Tgrnd = SVP_ice;         
        Qair_surf = 0.622 * (esat_Tgrnd * rh_grnd) / 
                            (pressure - 0.378 * (esat_Tgrnd * rh_grnd));
        moisture_flux = (Qair_surf - Qair) * coef_latent * 
                                         PsyCh_grnd / CONST_CPDAIR;
        NetLongGrnd = calc_longwave(Tgrnd, coef_longwave) - 
                                     longwave * EmissLongGrnd;
        SensibleGrnd = coef_sensible * (Tgrnd - air_temp);
        LatentGrnd = coef_latent * (esat_Tgrnd * rh_grnd - vp);
        GroundGrnd = NetShortGrnd + AdvectGrnd - 
                                NetLongGrnd - SensibleGrnd - 
                                                     LatentGrnd;
    }

    /********************************
      Compute grid mean quantities
    ********************************/
    energy->longwave = NetLongGrnd;
    energy->latent = LatentGrnd;
    energy->grnd_flux = GroundGrnd;
    energy->advection = AdvectGrnd;
    energy->sensible = SensibleGrnd;
    energy->Tsurf = Tgrnd;
    energy->Tgrnd = Tgrnd;
    energy->Tlower = Tgrnd;

    /************************************
      Compute glacier surf temperature
    ************************************/
    ErrorFlag = GlacierTemperature(step_dt, cell, energy,
                                   snow, soil_con);

    if (ErrorFlag == ERROR) {
        return (ERROR);
    }
    /*********************************
      Phase change of glacier ice
    **********************************/
    PhaseChangeGlac(step_dt, energy,
                    cell, snow, soil_con);

    return(0);
}
