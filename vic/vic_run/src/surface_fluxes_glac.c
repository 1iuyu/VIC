/******************************************************************************
* @section DESCRIPTION
*
* This routine computes all surface fluxes, and solves the snow accumulation
* and ablation algorithm. Solutions are for the current snow band and
* vegetation type.
******************************************************************************/

#include <vic_run.h>

/******************************************************************************
* @brief        This routine computes all surface fluxes
******************************************************************************/
int
surface_fluxes_glac(double               BareAlbedo,
//                    double               ice0,
//                    double               moist0,
                    double              *Melt,
                    double              *Le,
                    double              *aero_resist,
                    double              *displacement,
                    double              *gauge_correction,
                    double              *out_prec,
                    double              *out_rain,
                    double              *out_snow,
                    double              *ref_height,
                    double              *roughness,
                    double              *snow_inflow,
                    double              *wind,
                    size_t               Nlayers,
                    size_t               Nveg,
                    unsigned short       band,
                    unsigned short       iveg,
                    unsigned short       veg_class,
                    force_data_struct   *force,
                    dmy_struct          *dmy,
                    energy_bal_struct   *energy,
                    global_param_struct *gp,
                    cell_data_struct    *cell,
                    snow_data_struct    *snow,
                    soil_con_struct     *soil_con,
                    veg_var_struct      *veg_var,
                    glac_data_struct    *glacier,
                    double               lag_one,
                    double               sigma_slope,
                    double               fetch)
{
    extern veg_lib_struct   *vic_run_veg_lib;
    extern option_struct     options;
    extern parameters_struct param;

    int                      ErrorFlag;
    int                      N_steps;
    int                      UnderStory;
    size_t                   hidx; // index of initial element of atmos array
    size_t                   step_inc; // number of atmos array elements to skip per surface fluxes step
    size_t                   endhidx; // index of final element of atmos array
    double                   step_dt; // time length of surface fluxes step (in seconds)
    size_t                   lidx;
    int                      q;
    double                   Ls;
    double                   LongUnderIn; // inmoing LW to ground surface
    double                   LongUnderOut; // outgoing LW from ground surface
    double                   NetLongGlac; // net LW over snowpack
    double                   NetShortGlac; // net SW over understory
    double                   NetShortGrnd; // net SW over snow-free surface
    double                   OldTSurf; // previous snow surface temperature
    double                   ShortUnderIn; // incoming SW to understory
    double                   Tair; // air temperature
    double                   Tcanopy; // canopy air temperature
    double                   Tgrnd; // soil surface temperature
    double                   Tsurf; // ground surface temperature
    double                   VPDcanopy; // vapor pressure deficit in canopy/atmos
    double                   coverage; // mid-step snow cover fraction
    double                   delta_coverage; // change in snow cover fraction
    double                   last_snow_coverage; // previous snow covered area
    double                   ppt; // precipitation/melt reaching soil surface
    double                   rainfall; // rainfall
    double                   snowfall; // snowfall
    double                   snow_flux; // heat flux through snowpack
    double                  *inflow;
    layer_data_struct       *layer;
    double                  *aero_resist_used;


    // Step-specific quantities
    double                   step_Wdew;
    double                   step_melt;
    double                   step_melt_energy; /* energy used to reduce snow coverage */
    double                   step_out_prec = 0;
    double                   step_out_rain = 0;
    double                   step_out_snow = 0;
    double                   step_ppt;
//    double                   step_prec;
    double                   step_melt_glac;

    // Quantities that need to be summed or averaged over multiple snow steps
    // energy structure
    double            store_AlbedoOver;
    double            store_AlbedoUnder;
    double            store_AtmosLatent;
    double            store_AtmosLatentSub;
    double            store_AtmosSensible;
    double            store_LongOverIn;
    double            store_LongUnderIn;
    double            store_LongUnderOut;
    double            store_NetLongAtmos;
    double            store_NetLongOver;
    double            store_NetLongUnder;
    double            store_NetShortAtmos;
    double            store_NetShortGrnd;
    double            store_NetShortOver;
    double            store_NetShortUnder;
    double            store_ShortOverIn;
    double            store_ShortUnderIn;
    double            store_advected_sensible;
    double            store_advection;
    double            store_canopy_advection;
    double            store_canopy_latent;
    double            store_canopy_latent_sub;
    double            store_canopy_sensible;
    double            store_canopy_refreeze;
    double            store_deltaCC;
    double            store_deltaH;
    double            store_fusion;
    double            store_grnd_flux;
    double            store_latent;
    double            store_latent_sub;
    double            store_melt_energy;
    double            store_refreeze_energy;
    double            store_sensible;
    double            store_snow_flux;
    double            store_deltaCC_glac;
    double            store_glacier_flux;
    double            store_glacier_melt_energy;
    // glacier structure
    double            store_melt_glac;
    double            store_vapor_flux_glac;
    double            store_accum_glac;
    // snow structure
    double            store_canopy_vapor_flux;
    double            store_melt;
    double            store_vapor_flux;
    double            store_blowing_flux;
    double            store_surface_flux;
    // veg_var structure
    double            store_canopyevap;
    double            store_throughfall;
    // cell structure
    double            store_layerevap[MAX_LAYERS];
    double            store_ppt;
    double            store_aero_cond_used[2];
    double            store_pot_evap;

    // Structures holding values for current snow step
    energy_bal_struct step_energy;    // energy fluxes at snowpack surface
    veg_var_struct    snow_veg_var;    // veg fluxes/storages in presence of snow
    veg_var_struct    soil_veg_var;    // veg fluxes/storages in soil energy balance
    snow_data_struct  step_snow;
    layer_data_struct step_layer[MAX_LAYERS];
    glac_data_struct  step_glacier;

    // Structures holding values for current iteration
    double            temp_aero_resist[3];
    double            temp_aero_resist_used[2];
//    double            stability_factor[2];
    double            step_pot_evap;

    //step_aero_resist = new AeroResistUsed[N_PET_TYPES];

    /***********************************************************************
       Set temporary variables for convenience
    ***********************************************************************/
    aero_resist_used = cell->aero_resist;
    inflow = &(cell->inflow);
    layer = cell->layer;


    /***********************************************************************
       Set temporary variables - preserves original values until iterations
       are completed
    ***********************************************************************/

    coverage = snow->coverage;
    step_energy = (*energy);
    snow_veg_var = (*veg_var);
    soil_veg_var = (*veg_var);
    step_snow = (*snow);
    step_glacier = (*glacier);
    for (lidx = 0; lidx < Nlayers; lidx++) {
        step_layer[lidx] = layer[lidx];
    }
    for (lidx = 0; lidx < Nlayers; lidx++) {
        step_layer[lidx].evap = 0;
    }
    soil_veg_var.canopyevap = 0;
    snow_veg_var.canopyevap = 0;
    soil_veg_var.throughfall = 0;
    snow_veg_var.throughfall = 0;

    /********************************
       Set-up sub-time step controls
       (May eventually want to set this up so that it is also true
       if frozen soils are present)
    ********************************/

    hidx = 0;
    step_inc = 1;
    endhidx = hidx + NF;
    step_dt = gp->snow_dt;

    /*******************************************
       Initialize sub-model time step variables
    *******************************************/

    // energy structure
    store_AlbedoOver = 0;
    store_AlbedoUnder = 0;
    store_AtmosLatent = 0;
    store_AtmosLatentSub = 0;
    store_AtmosSensible = 0;
    store_LongOverIn = 0;
    store_LongUnderIn = 0;
    store_LongUnderOut = 0;
    store_NetLongAtmos = 0;
    store_NetLongOver = 0;
    store_NetLongUnder = 0;
    store_NetShortAtmos = 0;
    store_NetShortGrnd = 0;
    store_NetShortOver = 0;
    store_NetShortUnder = 0;
    store_ShortOverIn = 0;
    store_ShortUnderIn = 0;
    store_advected_sensible = 0;
    store_advection = 0;
    store_canopy_advection = 0;
    store_canopy_latent = 0;
    store_canopy_latent_sub = 0;
    store_canopy_sensible = 0;
    store_canopy_refreeze = 0;
    store_deltaCC = 0;
    store_deltaH = 0;
    store_fusion = 0;
    store_grnd_flux = 0;
    store_latent = 0;
    store_latent_sub = 0;
    store_melt_energy = 0;
    store_refreeze_energy = 0;
    store_sensible = 0;
    store_snow_flux = 0;
    // snow structure
    last_snow_coverage = snow->coverage;
    store_canopy_vapor_flux = 0;
    store_melt = 0;
    store_vapor_flux = 0;
    store_surface_flux = 0;
    store_blowing_flux = 0;
    // veg_var and cell structures
    store_throughfall = 0.;
    store_canopyevap = 0.;
    for (lidx = 0; lidx < options.Nlayer; lidx++) {
        store_layerevap[lidx] = 0.;
    }
    step_Wdew = veg_var->Wdew;
    // misc
    store_ppt = 0;
    //step_prec = 0;
    store_aero_cond_used[0] = 0;
    store_aero_cond_used[1] = 0;
    (*snow_inflow) = 0;
    store_pot_evap = 0;
    N_steps = 0;
    // glacier structure
    store_deltaCC_glac = 0;
    store_glacier_flux = 0;
    store_glacier_melt_energy = 0;
    store_melt_glac = 0;
    store_vapor_flux_glac = 0;
    store_accum_glac = 0;

    /*************************
       Compute surface fluxes
    *************************/

    do
    {
        /** Solve energy balance for all sub-model time steps **/

        /* set air temperature and precipitation for this snow band */
        Tair = force->air_temp[hidx] + soil_con->Tfactor[band];
        //step_prec = (force->rainf[hidx] + force->snowf[hidx]) * soil_con->Pfactor[band];
        snowfall = gauge_correction[SNOW] * force->snowf[hidx] * soil_con->PADJ_S * soil_con->Pfactor[band];
        rainfall = gauge_correction[RAIN] * force->rainf[hidx] * soil_con->PADJ_R * soil_con->Pfactor[band];

        step_out_prec = snowfall + rainfall;
        step_out_rain = rainfall;
        step_out_snow = snowfall;

        // initialize ground surface temperature: set to glacier temperature
        Tgrnd = GLAC_TEMP;

        // initialize canopy terms
        Tcanopy = 0.;
        VPDcanopy = 0.;

        // Compute mass flux of blowing snow
        if (options.BLOWING && step_snow.swq > 0.) {
            Ls = calc_latent_heat_of_sublimation(step_snow.surf_temp);
            step_snow.blowing_flux = CalcBlowingSnow(step_dt, Tair,
                                                     step_snow.last_snow,
                                                     step_snow.surf_water,
                                                     wind[2], Ls,
                                                     force->density[hidx],
                                                     force->vp[hidx],
                                                     roughness[2],
                                                     ref_height[2],
                                                     step_snow.depth,
                                                     lag_one, sigma_slope,
                                                     step_snow.surf_temp, iveg,
                                                     Nveg, fetch,
                                                     displacement[1],
                                                     roughness[1],
                                                     &step_snow.transport);
            if ((int) step_snow.blowing_flux == ERROR) {
                return (ERROR);
            }
            step_snow.blowing_flux *= step_dt / CONST_RHOFW; /* m/time step */
        }
        else {
            step_snow.blowing_flux = 0.0;
        }
        for (q = 0; q < 3; q++) {
            temp_aero_resist[q] = aero_resist[q];
        }
        temp_aero_resist_used[0] = cell->aero_resist[0];
        temp_aero_resist_used[1] = cell->aero_resist[1];
        step_snow.canopy_vapor_flux = 0;
        step_snow.vapor_flux = 0;
        step_snow.surface_flux = 0;

        LongUnderOut = step_energy.LongUnderOut;

        if (step_snow.swq > 0. || snowfall > 0.) {

            /** Solve snow accumulation, ablation and interception **/
            step_melt = solve_snow_glac(BareAlbedo,
                                        Tgrnd, Tair,
                                        &energy->AlbedoUnder, Le,
                                        &LongUnderIn, &NetLongGlac,
                                        &NetShortGlac,
                                        &ShortUnderIn, &OldTSurf,
                                        temp_aero_resist, temp_aero_resist_used,
                                        &coverage, &delta_coverage,
                                        displacement,
                                        &step_melt_energy,
                                        &step_ppt, &rainfall, ref_height,
                                        roughness, snow_inflow, &snowfall, wind,
                                        iveg, band, step_dt, hidx,
                                        &UnderStory, 
                                        dmy, force, &(step_energy),
                                        &(step_snow), soil_con,
                                        &(step_glacier));

            if (step_melt == ERROR) {
                return (ERROR);
            }

            step_melt_glac = 0.;
            step_glacier.vapor_flux = 0.;
            step_energy.glacier_flux = 0.;
            step_energy.deltaCC_glac = 0.;
            step_energy.glacier_melt_energy = 0.;
            step_energy.snow_flux = -step_energy.grnd_flux;
            step_energy.LongUnderOut = LongUnderIn - NetLongGlac;

        } else {
            /** Solve glacier is bare accumulation, ablation and interception **/
            step_melt_glac = solve_glacier(BareAlbedo, Tgrnd, Tair,
                                           &energy->AlbedoUnder, Le,
                                           &LongUnderIn, &NetLongGlac, 
                                           &NetShortGlac, &ShortUnderIn, 
                                           &OldTSurf, temp_aero_resist, 
                                           temp_aero_resist_used, 
                                           displacement, &step_melt_energy, 
                                           &step_ppt, &rainfall, ref_height, 
                                           roughness, wind, step_dt, 
                                           hidx, &UnderStory, force,
                                           &(step_energy), 
                                           &(step_glacier), soil_con,
                                           &(step_snow));


            if (step_melt_glac == ERROR) {

                return (ERROR);
            }

            step_melt = 0.;
            snow->snow = false;
            snow->store_swq = 0.;
            snow->store_coverage = 1;
            snow->last_snow = 0;
            snow->albedo = soil_con->NEW_SNOW_ALB;
            step_energy.deltaCC = 0.;
            step_energy.refreeze_energy = 0.;
            step_energy.snow_flux = 0.;
            step_energy.advected_sensible = 0.;
            step_energy.glacier_flux = -step_energy.grnd_flux;
            step_energy.LongUnderOut = LongUnderIn - NetLongGlac;
            step_glacier.accumulation = 0.;

        }

        //Put surface fluxes into atmospheric storage
        step_energy.AtmosLatent = step_energy.latent;
        step_energy.AtmosLatentSub = step_energy.latent_sub;
        step_energy.AtmosSensible = step_energy.sensible;
        step_energy.NetLongAtmos = step_energy.NetLongUnder;
        step_energy.NetShortAtmos = step_energy.NetShortUnder;


        /**************************************
           Compute Potential Evap
        **************************************/

        compute_pot_evap(gp->model_steps_per_day,
                         vic_run_veg_lib[veg_class].rmin,
                         soil_veg_var.albedo, force->shortwave[hidx],
                         step_energy.NetLongAtmos,
                         vic_run_veg_lib[veg_class].RGL, Tair, VPDcanopy,
                         soil_veg_var.LAI, soil_con->elevation,
                         &temp_aero_resist_used[0],
                         vic_run_veg_lib[veg_class].overstory,
                         vic_run_veg_lib[veg_class].rarc,
                         soil_veg_var.fcanopy, temp_aero_resist[0],
                         &step_pot_evap);

        /**************************************
           Store sub-model time step variables
        **************************************/

        store_ppt += step_ppt;
        if (temp_aero_resist_used[0] > 0) {
            store_aero_cond_used[0] += 1 / temp_aero_resist_used[0];
        }
        else {
            store_aero_cond_used[0] += param.HUGE_RESIST;
        }
        if (temp_aero_resist_used[1] > 0) {
            store_aero_cond_used[1] += 1 / temp_aero_resist_used[1];
        }
        else {
            store_aero_cond_used[1] += param.HUGE_RESIST;
        }

        store_melt += step_melt;
        store_vapor_flux += step_snow.vapor_flux;
        store_surface_flux += step_snow.surface_flux;
        store_blowing_flux += step_snow.blowing_flux;

        out_prec[0] += step_out_prec;
        out_rain[0] += step_out_rain;
        out_snow[0] += step_out_snow;

        store_AlbedoUnder += step_energy.AlbedoUnder;
        store_AtmosLatent += step_energy.AtmosLatent;
        store_AtmosLatentSub += step_energy.AtmosLatentSub;
        store_AtmosSensible += step_energy.AtmosSensible;
        store_LongUnderIn += LongUnderIn;
        store_LongUnderOut += step_energy.LongUnderOut;
        store_NetLongAtmos += NetLongGlac;
        store_NetLongUnder += NetLongGlac;
        store_NetShortAtmos += NetShortGlac;
        store_NetShortUnder += NetShortGlac;
        store_ShortUnderIn += ShortUnderIn;
        store_latent += step_energy.latent;
        store_latent_sub += step_energy.latent_sub;
        store_melt_energy += step_melt_energy;
        store_sensible += step_energy.sensible;
        store_grnd_flux += step_energy.grnd_flux;

        // glacier
        store_melt_glac += step_melt_glac;
        store_vapor_flux_glac += step_glacier.vapor_flux;
        store_accum_glac += step_glacier.accumulation;
        store_glacier_flux += step_energy.glacier_flux;
        store_deltaCC_glac += step_energy.deltaCC_glac;
        store_glacier_melt_energy += step_energy.glacier_melt_energy;

        store_advected_sensible += step_energy.advected_sensible *
                                    (step_snow.coverage + delta_coverage);
        store_advection += step_energy.advection *
                            (step_snow.coverage + delta_coverage);
        store_deltaCC += step_energy.deltaCC *
                            (step_snow.coverage + delta_coverage);
        store_snow_flux += step_energy.snow_flux *
                            (step_snow.coverage + delta_coverage);
        store_refreeze_energy += step_energy.refreeze_energy *
                                    (step_snow.coverage + delta_coverage);

        store_pot_evap += step_pot_evap;

        /* increment time step */
        N_steps++;
        hidx += step_inc;

    }
    while (hidx < endhidx);

    /************************************************
      Store glacier variables for sub-model time steps
    ************************************************/
    (*glacier) = step_glacier;
    glacier->melt = store_melt_glac;
    glacier->vapor_flux = store_vapor_flux_glac;
    glacier->accumulation = store_accum_glac;

    /************************************************
       Store snow variables for sub-model time steps
    ************************************************/

    (*snow) = step_snow;
    snow->vapor_flux = store_vapor_flux;
    snow->blowing_flux = store_blowing_flux;
    snow->surface_flux = store_surface_flux;
    snow->canopy_vapor_flux = store_canopy_vapor_flux;
    (*Melt) = store_melt;
    snow->melt = store_melt;
    ppt = store_ppt;

    glacier->mass_balance = (*out_prec) / MM_PER_M - ppt - 
                           snow->vapor_flux - glacier->vapor_flux;

    glacier->ice_mass_balance = glacier->accumulation - 
                               glacier->melt - glacier->vapor_flux;

    /******************************************************
       Store energy flux averages for sub-model time steps
    ******************************************************/

    (*energy) = step_energy;
    energy->AlbedoOver = store_AlbedoOver / (double) N_steps;
    energy->AlbedoUnder = store_AlbedoUnder / (double) N_steps;
    energy->AtmosLatent = store_AtmosLatent / (double) N_steps;
    energy->AtmosLatentSub = store_AtmosLatentSub / (double) N_steps;
    energy->AtmosSensible = store_AtmosSensible / (double) N_steps;
    energy->LongOverIn = store_LongOverIn / (double) N_steps;
    energy->LongUnderIn = store_LongUnderIn / (double) N_steps;
    energy->LongUnderOut = store_LongUnderOut / (double) N_steps;
    energy->NetLongAtmos = store_NetLongAtmos / (double) N_steps;
    energy->NetLongOver = store_NetLongOver / (double) N_steps;
    energy->NetLongUnder = store_NetLongUnder / (double) N_steps;
    energy->NetShortAtmos = store_NetShortAtmos / (double) N_steps;
    energy->NetShortGrnd = store_NetShortGrnd / (double) N_steps;
    energy->NetShortOver = store_NetShortOver / (double) N_steps;
    energy->NetShortUnder = store_NetShortUnder / (double) N_steps;
    energy->ShortOverIn = store_ShortOverIn / (double) N_steps;
    energy->ShortUnderIn = store_ShortUnderIn / (double) N_steps;
    energy->advected_sensible = store_advected_sensible / (double) N_steps;
    energy->canopy_advection = store_canopy_advection / (double) N_steps;
    energy->canopy_latent = store_canopy_latent / (double) N_steps;
    energy->canopy_latent_sub = store_canopy_latent_sub / (double) N_steps;
    energy->canopy_refreeze = store_canopy_refreeze / (double) N_steps;
    energy->canopy_sensible = store_canopy_sensible / (double) N_steps;
    energy->deltaH = store_deltaH / (double) N_steps;
    energy->fusion = store_fusion / (double) N_steps;
    energy->grnd_flux = store_grnd_flux / (double) N_steps;
    energy->latent = store_latent / (double) N_steps;
    energy->latent_sub = store_latent_sub / (double) N_steps;
    energy->melt_energy = store_melt_energy / (double) N_steps;
    energy->sensible = store_sensible / (double) N_steps;
    energy->glacier_flux = store_glacier_flux / (double) N_steps;
    energy->deltaCC_glac = store_deltaCC_glac / (double) N_steps;
    energy->glacier_melt_energy = store_glacier_melt_energy / (double) N_steps;
    energy->advection = store_advection / (double) N_steps;
    energy->deltaCC = store_deltaCC / (double) N_steps;
    energy->refreeze_energy = store_refreeze_energy / (double) N_steps;
    energy->snow_flux = store_snow_flux / (double) N_steps;
    energy->Tfoliage = step_energy.Tfoliage;
    energy->Tfoliage_fbflag = step_energy.Tfoliage_fbflag;
    energy->Tfoliage_fbcount = step_energy.Tfoliage_fbcount;
    energy->Tcanopy = Tcanopy;


    /**********************************************************
       Store vegetation variable sums for sub-model time steps
    **********************************************************/

    if (iveg != Nveg) {
        veg_var->throughfall = store_throughfall;
        veg_var->canopyevap = store_canopyevap;
        if (snow->snow) {
            veg_var->Wdew = snow_veg_var.Wdew;
        }
        else {
            veg_var->Wdew = soil_veg_var.Wdew;
        }
    }

    /**********************************************************
       Store soil layer variables for sub-model time steps
    **********************************************************/

    for (lidx = 0; lidx < Nlayers; lidx++) {
        layer[lidx] = step_layer[lidx];
        layer[lidx].evap = store_layerevap[lidx];
    }
    if (store_aero_cond_used[0] > 0 && store_aero_cond_used[0] <
        param.HUGE_RESIST) {
        aero_resist_used[0] = 1 / (store_aero_cond_used[0] / (double) N_steps);
    }
    else if (store_aero_cond_used[0] >= param.HUGE_RESIST) {
        aero_resist_used[0] = 0;
    }
    else {
        aero_resist_used[0] = param.HUGE_RESIST;
    }
    if (store_aero_cond_used[1] > 0 && store_aero_cond_used[1] <
        param.HUGE_RESIST) {
        aero_resist_used[1] = 1 / (store_aero_cond_used[1] / (double) N_steps);
    }
    else if (store_aero_cond_used[1] >= param.HUGE_RESIST) {
        aero_resist_used[1] = 0;
    }
    else {
        aero_resist_used[1] = param.HUGE_RESIST;
    }
    cell->pot_evap = store_pot_evap;


    /********************************************************
       Compute Runoff, Baseflow, and Soil Moisture Transport
    ********************************************************/

    (*inflow) = ppt;

    glacier->outflow_coef = soil_con->GLAC_KMIN + soil_con->GLAC_DK * exp(-soil_con->GLAC_A * snow->swq);
    glacier->water_storage += glacier->inflow;
    glacier->outflow = glacier->outflow_coef * glacier->water_storage;
    glacier->water_storage -= glacier->outflow;

    ErrorFlag = runoff(cell, energy, soil_con, ppt, soil_con->frost_fract,
                       options.Nnode);

    cell->runoff += (glacier->outflow * MM_PER_M); /* convert to mm */

    return(ErrorFlag);

}
