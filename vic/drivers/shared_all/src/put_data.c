/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine converts data units, and stores finalized values in an array
 * for later output to the output files.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    This routine converts data units, and stores finalized values in
 *           an array for later output to the output files.
 *****************************************************************************/
void
put_data(all_vars_struct   *all_vars,
         force_data_struct *force,
         soil_con_struct   *soil_con,
         veg_con_struct    *veg_con,
         veg_lib_struct    *veg_lib,
         double           **out_data,
         save_data_struct  *save_data,
         timer_struct      *timer)
{
    extern global_param_struct global_param;
    extern option_struct       options;
    extern parameters_struct   param;

    size_t             i, veg;
    size_t             lidx;
    size_t             band;
    bool               HasVeg;
    bool               HasGlac;
    double             Cv;
    double             cv_baresoil;
    double             cv_veg;
    double             cv_snow;
    double             cv_glac;
    double             inflow;
    double             outflow;
    double             storage;
    double             dt_sec;

    cell_data_struct  *cell;
    energy_bal_struct *energy;
    snow_data_struct  *snow;
    veg_var_struct    *veg_var;

    cell = all_vars->cell;
    energy = all_vars->energy;
    snow = all_vars->snow;
    veg_var = all_vars->veg_var;

    dt_sec = global_param.dt;

    cv_baresoil = 0;
    cv_veg = 0;
    cv_snow = 0;
    cv_glac = 0;

    // Initialize output data to zero
    zero_output_list(out_data);

    // Set output versions of input forcings
    out_data[OUT_AIR_TEMP][0] = force->air_temp[NR];
    out_data[OUT_DENSITY][0] = force->density[NR];
    out_data[OUT_LWDOWN][0] = force->longwave[NR];
    out_data[OUT_PREC][0] = force->prec[NR] * dt_sec;  // mm over grid cell
    out_data[OUT_PRESSURE][0] = force->pressure[NR] / PA_PER_KPA;
    out_data[OUT_QAIR][0] = force->Qair[NR];
    out_data[OUT_RAINF][0] = force->rainf[NR] * dt_sec;   // mm over grid cell
    out_data[OUT_REL_HUMID][0] = force->rel_humid[NR];
    out_data[OUT_SWDOWN][0] = force->shortwave[NR];
    out_data[OUT_SNOWF][0] = force->snowf[NR] * dt_sec;   // mm over grid cell
    out_data[OUT_VP][0] = force->vp[NR] / PA_PER_KPA;
    out_data[OUT_WIND][0] = force->wind[NR];
    out_data[OUT_COSZEN][0] = force->coszen[NR];
    if (options.CARBON) {
        out_data[OUT_CATM][0] = force->Catm[NR] / PPM_to_MIXRATIO;
        out_data[OUT_FDIR][0] = force->fdir[NR];
        out_data[OUT_PAR][0] = force->par[NR];
    }
    else {
        out_data[OUT_CATM][0] = MISSING;
        out_data[OUT_FDIR][0] = MISSING;
        out_data[OUT_PAR][0] = MISSING;
    }
    /** compute running totals of various landcovers **/
    for (veg = 0; veg <= veg_con[0].vegetat_type_num; veg++) {
        band = veg_con[veg].BandIndex;
        if (veg < veg_con[0].vegetat_type_num) {
            cv_veg += veg_con[veg].Cv;
        }
        else {
            cv_baresoil += veg_con[veg].Cv;
        }
        if (snow[veg].swq > 0.0) {
            cv_snow += veg_con[veg].Cv;
        }
        if (veg_con[veg].IS_GLAC) {
            cv_glac += veg_con[veg].Cv;
        }
    }

    /****************************************
       Store Output for all Vegetation Types (except lakes)
    ****************************************/
    for (veg = 0; veg <= veg_con[0].vegetat_type_num; veg++) {
        Cv = veg_con[veg].Cv;
        if (veg < veg_con[0].vegetat_type_num) {
            HasVeg = true;
        }
        else {
            HasVeg = false;
        }

        HasGlac = veg_con[veg].IS_GLAC;

        if (Cv > 0) {

            band = veg_con[veg].BandIndex;

            /*********************************
               Record Water Balance Terms
            *********************************/
            collect_wb_terms(cell[veg],
                             veg_var[veg],
                             snow[veg],
                             Cv,
                             HasVeg,
                             dt_sec,
                             out_data);

            /**********************************
               Record Energy Balance Terms
            **********************************/
            collect_eb_terms(energy[veg],
                             snow[veg],
                             cell[veg],
                             Cv,
                             HasVeg,
                             dt_sec,
                             band,
                             out_data);
        } // End if Cv > 0
    } // End loop over veg


    /*****************************************
       Finish aggregation of special-case variables
    *****************************************/
    // Normalize quantities that aren't present over entire grid cell
    if (cv_baresoil > 0) {
        out_data[OUT_BARESOILT][0] /= cv_baresoil;
    }
    if (cv_veg > 0) {
        out_data[OUT_VEGT][0] /= cv_veg;
    }
    if (cv_snow > 0) {
        out_data[OUT_SNOW_PACK_TEMP][0] /= cv_snow;
    }

    // Radiative temperature
    out_data[OUT_RAD_TEMP][0] = pow(out_data[OUT_RAD_TEMP][0], 0.25);

    /*****************************************
       Compute derived variables
    *****************************************/
    // Water balance terms
    out_data[OUT_DELSOILMOIST][0] = 0;
    for (lidx = 0; lidx < options.Nnode; lidx++) {
        out_data[OUT_SOIL_MOIST][lidx] =
            out_data[OUT_SOIL_LIQ][lidx] +
            out_data[OUT_SOIL_ICE][lidx];
        out_data[OUT_DELSOILMOIST][0] +=
            out_data[OUT_SOIL_MOIST][lidx];

        out_data[OUT_SMLIQFRAC][lidx] = out_data[OUT_SOIL_LIQ][lidx] /
                                         out_data[OUT_SOIL_MOIST][lidx];
        out_data[OUT_SMFROZFRAC][lidx] = 1 - out_data[OUT_SMLIQFRAC][lidx];
    }
    out_data[OUT_DELSOILMOIST][0] -= save_data->total_soil_moist;
    out_data[OUT_DELSWE][0] = out_data[OUT_SWE][0] +
                              out_data[OUT_SNOW_CANOPY][0] -
                              save_data->swe;
    out_data[OUT_DELINTERCEPT][0] = out_data[OUT_WDEW][0] - save_data->wdew;

    // Save current moisture state for use in next time step
    save_data->total_soil_moist = 0;
    for (lidx = 0; lidx < options.Nnode; lidx++) {
        save_data->total_soil_moist += out_data[OUT_SOIL_MOIST][lidx];
    }
    save_data->swe = out_data[OUT_SWE][0] + out_data[OUT_SNOW_CANOPY][0];
    save_data->wdew = out_data[OUT_WDEW][0];

    // Carbon Terms
    if (options.CARBON) {
        out_data[OUT_RHET][0] *= dt_sec / SEC_PER_DAY;  // convert to gC/m2d
        out_data[OUT_NEE][0] = out_data[OUT_NPP][0] - out_data[OUT_RHET][0];
    }

    /********************
       Check Water Balance
    ********************/
    inflow = out_data[OUT_PREC][0];  // mm over grid cell
    outflow = out_data[OUT_EVAP][0] + out_data[OUT_RUNOFF][0] +
              out_data[OUT_BASEFLOW][0];  // mm over grid cell

    storage = 0.;
    for (lidx = 0; lidx < options.Nnode; lidx++) {
        storage += out_data[OUT_SOIL_LIQ][lidx] +
                   out_data[OUT_SOIL_ICE][lidx];
    }
    storage += out_data[OUT_SWE][0] + out_data[OUT_SNOW_CANOPY][0] +
               out_data[OUT_WDEW][0];

    out_data[OUT_WATER_ERROR][0] = \
        calc_water_balance_error(inflow,
                                 outflow,
                                 storage,
                                 save_data->
                                 total_moist_storage);

    // Store total storage for next timestep
    save_data->total_moist_storage = storage;

    /********************
       Check Energy Balance
    ********************/
    out_data[OUT_ENERGY_ERROR][0] = MISSING;


    // vic_run run time
    out_data[OUT_TIME_VICRUN_WALL][0] = timer->delta_wall;
    out_data[OUT_TIME_VICRUN_CPU][0] = timer->delta_cpu;
}

/******************************************************************************
 * @brief    This routine collects water balance terms.
 *****************************************************************************/
void
collect_wb_terms(cell_data_struct cell,
                 veg_var_struct   veg_var,
                 snow_data_struct snow,
                 double           Cv,
                 bool             HasVeg,
                 double           dt_sec,
                 double         **out_data)
{
    extern option_struct     options;

    size_t      lidx;
    size_t      Nsnow;
    /*********************************
      record evaporation components
    *********************************/

    // Ground (soil/snow) surface
    double tmp_evap = 0.0;
    tmp_evap += cell.evap;
    out_data[OUT_NET_EVAP][0] += cell.evap * Cv * dt_sec;
    out_data[OUT_EVAP_BARE][0] += cell.esoil * Cv * dt_sec;
    out_data[OUT_DEW_SOIL][0] += cell.dewsoil * Cv * dt_sec;
    out_data[OUT_VAPOR_GRND][0] += cell.vapor_grnd * Cv * dt_sec;
    out_data[OUT_CONDEN_GRND][0] += cell.conden_grnd * Cv * dt_sec;
    out_data[OUT_FROST_SNOW][0] += cell.snowfrost * Cv * dt_sec;
    out_data[OUT_SUB_SNOW][0] += cell.snow_sublim * Cv * dt_sec;

    // Canopy evaporation
    if (HasVeg) {
        tmp_evap += cell.canopyevap;
        tmp_evap += cell.transp;
        out_data[OUT_EVAP_CANOP][0] += cell.canopyevap * Cv *dt_sec;
        out_data[OUT_TRANSP_VEG][0] += cell.transp * Cv * dt_sec;
        out_data[OUT_DEW_CANOP][0] += cell.canopydew * Cv * dt_sec;
        out_data[OUT_FROST_CANOP][0] += cell.canopyfrost * Cv * dt_sec;
        out_data[OUT_SUB_CANOP][0] += cell.canopy_sublim * Cv * dt_sec;
        out_data[OUT_VAPOR_CANOP][0] += cell.canopy_vapor * Cv * dt_sec;
    }
    // net evaporation = evap + canopyevap + transp
    out_data[OUT_EVAP][0] += tmp_evap * Cv * dt_sec;  // mm over gridcell

    /** record saturated area fraction **/
    out_data[OUT_ASAT][0] += cell.asat * Cv;

    /** record runoff **/
    out_data[OUT_RUNOFF][0] += cell.runoff * Cv;

    /** record baseflow **/
    out_data[OUT_BASEFLOW][0] += cell.baseflow * Cv;

    /** record soil_inflow **/
    out_data[OUT_INFLOW][0] += cell.soil_inflow * Cv * dt_sec * MM_PER_M;

    /** record canopy interception **/
    if (HasVeg) {
        out_data[OUT_WDEW][0] += veg_var.Wdew * Cv;
    }

    /** record LAI **/
    out_data[OUT_LAI][0] += veg_var.LAI * Cv;

    /** record fcanopy **/
    out_data[OUT_FCANOPY][0] += veg_var.fcanopy * Cv;

    /** record aerodynamic conductance and resistance **/
    for (lidx = 0; lidx < 3; lidx++) {
        out_data[OUT_RA_OVER][0] += cell.Ra_over[lidx];
        out_data[OUT_RA_GRND][0] += cell.Ra_grnd[lidx];
        out_data[OUT_RA_SUB][0] += cell.Ra_sub[lidx];
    }
    out_data[OUT_RA_EVAP][0] += cell.Ra_evap;
    if (HasVeg) {    
        out_data[OUT_RA_LEAF][0] += cell.Ra_leaf;
    }

    /** record nodes moistures **/
    for (lidx = 0; lidx < options.Nnode; lidx++) {
        out_data[OUT_SOIL_LIQ][lidx] += cell.liq[lidx] * Cv;
        out_data[OUT_SOIL_ICE][lidx] += cell.ice[lidx] * Cv;
    }
    out_data[OUT_SOIL_WET][0] += cell.wetness * Cv;
    out_data[OUT_ROOTMOIST][0] += cell.rootmoist * Cv;

    /** record water table position **/
    out_data[OUT_ZWT][0] += cell.zwt * Cv;
    out_data[OUT_ZWT_LUMPED][0] += cell.zwt_lumped * Cv;

    /*****************************
       Record Snow Pack Variables
    *****************************/

    /** record snow water equivalence **/
    out_data[OUT_SWE][0] += snow.swq * Cv;

    /** record snowpack depth **/
    out_data[OUT_SNOW_DEPTH][0] += snow.snow_depth * Cv * CM_PER_M;

    /** record snowpack temperature water and ice content **/
    Nsnow = snow.Nsnow;
    for (lidx = 0; lidx < Nsnow; lidx++) {
        out_data[OUT_SNOW_PACK_TEMP][lidx] += snow.pack_T[lidx] * Cv;
        out_data[OUT_SNOW_PACK_ICE][lidx] += snow.pack_ice[lidx] * Cv;
        out_data[OUT_SNOW_PACK_LIQ][lidx] += snow.pack_liq[lidx] * Cv;
        out_data[OUT_PACK_OUTFLOW][lidx] += snow.pack_outflow[lidx] * Cv;
        out_data[OUT_SNOW_ICEFRAC][lidx] += snow.theta_ice[lidx] * Cv;
        out_data[OUT_SNOW_LIQFRAC][lidx] += snow.theta_liq[lidx] * Cv;
        out_data[OUT_SNOW_POROSITY][lidx] += snow.porosity[lidx] * Cv;
    }
    /** record canopy intercepted snow **/
    if (HasVeg) {
        out_data[OUT_SNOW_CANOPY][0] += (snow.snow_canopy) * Cv *
                                        MM_PER_M;
    }
    /** record snowpack melt **/
    out_data[OUT_SNOW_MELT][0] += snow.pack_melt * Cv;
    /** record snowpack transp **/
    out_data[OUT_SNOW_TRANSP][0] += snow.pack_transp * Cv;
    /** record snowpack combination **/
    out_data[OUT_SNOW_COMB][0] += snow.pack_comb * Cv;
    /** record glacier snow excess flow **/
    out_data[OUT_GLAC_EXCESS][0] += snow.glac_excess * Cv;
    /** record snow cover fraction **/
    out_data[OUT_SNOW_COVER][0] += snow.coverage * Cv;
    /** record snow depth increasing rate **/
    out_data[OUT_DELDEPTH][0] += snow.delta_depth * Cv;
    /** record snow density **/
    out_data[OUT_SNOW_DENSITY][0] += snow.density * Cv;
    /** record NEW snow density **/
    out_data[OUT_NEW_DENSITY][0] += snow.new_snow_density * Cv;
    /** record SnowAge **/
    out_data[OUT_SNOW_AGE][0] += snow.SnowAge * Cv;

    /*****************************
       Record Carbon Cycling Variables
    *****************************/
    if (options.CARBON) {
        out_data[OUT_APAR][0] += veg_var.aPAR * Cv;
        out_data[OUT_GPP][0] += veg_var.GPP * CONST_MWC / MOLE_PER_KMOLE *
                                CONST_CDAY * Cv;
        out_data[OUT_RAUT][0] += veg_var.Raut * CONST_MWC /
                                 MOLE_PER_KMOLE * CONST_CDAY * Cv;
        out_data[OUT_NPP][0] += veg_var.NPP * CONST_MWC / MOLE_PER_KMOLE *
                                CONST_CDAY * Cv;
        out_data[OUT_LITTERFALL][0] += veg_var.Litterfall * Cv;
        out_data[OUT_RHET][0] += cell.RhTot * Cv;
        out_data[OUT_CLITTER][0] += cell.CLitter * Cv;
        out_data[OUT_CINTER][0] += cell.CInter * Cv;
        out_data[OUT_CSLOW][0] += cell.CSlow * Cv;
    }
}

/******************************************************************************
 * @brief    This routine collects energy balance terms.
 *****************************************************************************/
void
collect_eb_terms(energy_bal_struct energy,
                 snow_data_struct  snow,
                 cell_data_struct  cell_wet,
                 double            Cv,
                 bool              HasVeg,
                 double            dt_sec,
                 int               band,
                 double          **out_data)
{
    extern option_struct options;
    double               tmp_fract;
    size_t               lidx;

    /**************************************
       Record Frozen Grnd and Canopy flag
    **************************************/
    out_data[OUT_TRND_FBFLAG][0] += energy.FrozenGrnd;
    out_data[OUT_TCAN_FBFLAG][0] += energy.FrozenOver;

    /**********************************
       Record Energy Balance Variables
    **********************************/

    /** record landcover temperature **/
    if (HasVeg) {
        // landcover is vegetation
        out_data[OUT_VEGT][0] += energy.Tcanopy * Cv;
        out_data[OUT_VEGTAIR][0] += energy.Tair * Cv;
    }
    // landcover is bare soil
    out_data[OUT_BARESOILT][0] += energy.Tgrnd * Cv;
    /** record mean surface temperature [C] **/
    out_data[OUT_SURF_TEMP][0] += energy.Tsurf * Cv;

    /** record NODES temperatures **/
    size_t Nsnow = snow.Nsnow;
    for (lidx = 0; lidx < options.Nnode; lidx++) {
        out_data[OUT_SOIL_TEMP][lidx] += energy.T[lidx + Nsnow] * Cv;
    }
    /** record advective flux from prec **/
    out_data[OUT_ADVECTION][0] += energy.advection * Cv;
    out_data[OUT_ADVECTSUB][0] += energy.AdvectSub * Cv;
    out_data[OUT_ADVECTGRND][0] += energy.AdvectGrnd * Cv;
    out_data[OUT_ADVECTOVER][0] += energy.AdvectOver * Cv;
    /** record net shortwave radiation **/
    out_data[OUT_SWNET][0] += energy.shortwave * Cv;
    out_data[OUT_SWSUB][0] += energy.NetShortSub * Cv;
    out_data[OUT_SWGRND][0] += energy.NetShortGrnd * Cv;
    /** record longwave radiation flux **/
    out_data[OUT_LWNET][0] += energy.longwave * Cv;
    out_data[OUT_LWSUB][0] += energy.NetLongSub * Cv;
    out_data[OUT_LWGRND][0] += energy.NetLongGrnd * Cv;
    out_data[OUT_LWOVER][0] += energy.NetLongOver * Cv;
    /** record latent heat flux **/
    out_data[OUT_LATENT][0] += energy.latent * Cv;
    out_data[OUT_LATENT_SUB][0] += energy.LatentSub * Cv;
    out_data[OUT_LATENT_GRND][0] += energy.LatentGrnd * Cv;
    out_data[OUT_LATENT_CANOP][0] += energy.LatentCanopy * Cv;
    out_data[OUT_LATENT_TRANSP][0] += energy.LatentTransp * Cv;
    /** record sensible heat flux **/
    out_data[OUT_SENSIBLE][0] += energy.sensible * Cv;
    out_data[OUT_SENSIBLE_SUB][0] += energy.SensibleSub * Cv;
    out_data[OUT_SENSIBLE_GRND][0] += energy.SensibleGrnd * Cv;
    out_data[OUT_SENSIBLE_OVER][0] += energy.SensibleOver * Cv;
    /** record ground heat flux **/
    out_data[OUT_GRND_FLUX][0] += energy.grnd_flux * Cv;
    out_data[OUT_GRND_SUB][0] += energy.GroundSub * Cv;
    out_data[OUT_GRND_GRND][0] += energy.GroundGrnd * Cv;
    /** record heat of fusion **/
    //out_data[OUT_FUSION][0] += energy.fusion_flux * Cv;

    /** record radiative effective temperature [K],
        emissivities set = 1.0  **/


    /**********************************
       Record hru-Specific Variables
    **********************************/

    /** record hru snow water equivalent **/
    out_data[OUT_SWE_BAND][band] += snow.swq * Cv * MM_PER_M;

    /** record band snowpack depth **/
    out_data[OUT_SNOW_DEPTH_BAND][band] += snow.snow_depth * Cv *
                                           CM_PER_M;

    /** record band canopy intercepted snow **/
    if (HasVeg) {
        out_data[OUT_SNOW_CANOPY_BAND][band] += snow.snow_canopy * Cv * MM_PER_M;
    }

    /** record band snow melt **/
    out_data[OUT_SNOW_MELT_BAND][band] += snow.pack_melt * Cv;

    /** record band snow coverage **/
    out_data[OUT_SNOW_COVER_BAND][band] += snow.coverage * Cv;

    /** record band advection **/
    out_data[OUT_ADVECTION_BAND][band] += energy.advection * Cv;

    /** record pack layer temperature **/

    /** record latent heat of sublimation **/
    out_data[OUT_LATENT_SUB_BAND][band] += energy.LatentSub * Cv;

    /** record band net downwards longwave radiation **/
    out_data[OUT_LWNET_BAND][band] += energy.longwave * Cv;

    /** record band net latent heat flux **/
    out_data[OUT_LATENT_BAND][band] -= energy.latent * Cv;

    /** record band net sensible heat flux **/
    out_data[OUT_SENSIBLE_BAND][band] -= energy.sensible * Cv;

    /** record band net ground heat flux **/
    out_data[OUT_GRND_FLUX_BAND][band] -= energy.grnd_flux * Cv;
}

/******************************************************************************
 * @brief    Initialize the save data structure.
 *****************************************************************************/
void
initialize_save_data(all_vars_struct   *all_vars,
                     force_data_struct *atmos,
                     soil_con_struct   *soil_con,
                     veg_con_struct    *veg_con,
                     veg_lib_struct    *veg_lib,
                     double           **out_data,
                     save_data_struct  *save_data,
                     timer_struct      *timer)
{
    // Calling put data will populate the save data storage terms
    put_data(all_vars, atmos, soil_con, veg_con, veg_lib,
             out_data, save_data, timer);

    zero_output_list(out_data);
}

/******************************************************************************
 * @brief    This subroutine computes the overall model water balance, and
 *           warns the model user if large errors are found.
 *****************************************************************************/
double
calc_water_balance_error(double inflow,
                         double outflow,
                         double storage,
                         double last_storage)
{
    double error;

    error = inflow - outflow - (storage - last_storage);

    return(error);
}

/******************************************************************************
 * @brief    This subroutine computes the overall model energy balance.
 *****************************************************************************/
double
calc_energy_balance_error(double net_rad,
                          double latent,
                          double sensible,
                          double grnd_flux,
                          double advection)
{
    double error;

    error = net_rad - latent - sensible - grnd_flux + advection;

    return(error);
}
