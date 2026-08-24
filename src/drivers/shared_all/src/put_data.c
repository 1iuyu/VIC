/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine converts data units, and stores finalized values in an array
 * for later output to the output files.
 *****************************************************************************/

#include "vic_driver_shared_all.h"

/******************************************************************************
 * @brief    This routine converts data units, and stores finalized values in
 *           an array for later output to the output files.
 *****************************************************************************/
void
put_data(all_vars_struct   *all_vars,
         force_data_struct *force,
         veg_con_struct    *veg_con,
         double           **out_data,
         timer_struct      *timer)
{
    extern global_param_struct global_param;
    extern option_struct       options;

    size_t             veg;
    size_t             band;
    bool               HasVeg;
    bool               HasGlac;
    double             Cv;
    double             cv_baresoil;
    double             cv_veg;
    double             cv_snow;
    double             cv_glac;
    double             dt_sec;

    cell_data_struct  *cell;
    energy_bal_struct *energy;
    snow_data_struct  *snow;
    veg_var_struct    *veg_var;

    cell = all_vars->cell;
    energy = all_vars->energy;
    snow = all_vars->snow;
    veg_var = all_vars->veg_var;

    dt_sec = global_param.step_dt;

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
        if (veg < veg_con[0].vegetat_type_num && 
            (!cell[veg].IS_GLAC || !cell[veg].IS_URBAN || !cell[veg].IS_WET)) {
            cv_veg += veg_con[veg].Cv;
        }
        else {
            cv_baresoil += veg_con[veg].Cv;
        }
        if (snow[veg].swq > 0.0) {
            cv_snow += veg_con[veg].Cv;
        }
        if (cell[veg].IS_GLAC) {
            cv_glac += veg_con[veg].Cv;
        }
    }

    /****************************************
       Store Output for all Vegetation Types
    ****************************************/
    for (veg = 0; veg <= veg_con[0].vegetat_type_num; veg++) {
        Cv = veg_con[veg].Cv;
        if (veg < veg_con[0].vegetat_type_num && 
            (!cell[veg].IS_GLAC || !cell[veg].IS_URBAN || !cell[veg].IS_WET)) {
            HasVeg = true;
        }
        else {
            HasVeg = false;
        }

        HasGlac = cell[veg].IS_GLAC;

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
                             HasGlac,
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
                             HasGlac,
                             band,
                             out_data);
        } // End if Cv > 0
    } // End loop over veg

    /*****************************************
       Finish aggregation of special-case variables
    *****************************************/
    // Normalize quantities that aren't present over entire grid cell
    if (cv_veg > 0) {
        out_data[OUT_VEGT][0] /= cv_veg;
        out_data[OUT_VEGTAIR][0] /= cv_veg;
        out_data[OUT_ZWT][0] /= cv_veg;
    }
    if (cv_snow > 0) {
        out_data[OUT_SNOW_PACK_TEMP][0] /= cv_snow;
    }
    if (cv_glac > 0.0) {
        out_data[OUT_GLAC_INFLOW][0] /= cv_glac;
        out_data[OUT_GLAC_OUTFLOW][0] /= cv_glac;
    }

    // Radiative temperature
    out_data[OUT_RAD_TEMP][0] = pow(out_data[OUT_RAD_TEMP][0], 0.25);

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
                 bool             HasGlac,
                 double           dt_sec,
                 double         **out_data)
{
    extern option_struct     options;

    size_t i, Nsnow;
    /*********************************
      record evaporation components
    *********************************/

    // Ground (soil/snow) surface
    double tmp_evap = 0.0;
    tmp_evap += cell.evap;
    out_data[OUT_NET_EVAP][0] += cell.evap * Cv * dt_sec;

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
    
    /** record recharge to groundwater storage */
    out_data[OUT_RECHARGE][0] += cell.recharge * Cv;

    /** record soil surface dew rate [mm/s] */
    out_data[OUT_SOIL_DEW][0] += cell.soil_dew * Cv * dt_sec;

    /** record soil surface frost rate [mm/s] */
    out_data[OUT_SOIL_FROST][0] += cell.soil_frost * Cv * dt_sec;

    /** record soil evaporation from soil layer [mm/s] */
    out_data[OUT_SOIL_EVAP][0] += cell.soil_evap * Cv * dt_sec;

    /** record soil surface sublimation rate [mm/s] */
    out_data[OUT_SOIL_SUBLIM][0] += cell.soil_sublim * Cv * dt_sec;

    /** record soil_inflow[mm] **/
    out_data[OUT_INFLOW][0] += cell.soil_inflow * Cv * dt_sec * MM_PER_M;

    /** record canopy interception **/
    if (HasVeg) {
        out_data[OUT_WDEW][0] += veg_var.Wdew * Cv;
        out_data[OUT_INT_SNOW][0] += veg_var.int_snow * Cv;
        out_data[OUT_INT_RAIN][0] += veg_var.int_rain * Cv;
        out_data[OUT_CANOPY_SWQ][0] += veg_var.canopy_swq * Cv;
        out_data[OUT_SNOW_DRIP][0] += veg_var.SnowDrip * Cv;
        out_data[OUT_RAIN_DRIP][0] += veg_var.RainDrip * Cv;
        out_data[OUT_SNOWTHROUGHFALL][0] += veg_var.SnowThroughFall * Cv;
        out_data[OUT_RAINTHROUGHFALL][0] += veg_var.RainThroughFall * Cv;
    }

    /** record LAI **/
    out_data[OUT_LAI][0] += veg_var.LAI * Cv;

    /** record fcanopy **/
    out_data[OUT_FCANOPY][0] += veg_var.fcanopy * Cv;

    /** record aerodynamic conductance and resistance **/
    for (i = 0; i < 3; i++) {
        out_data[OUT_RA_OVER][0] += cell.Ra_over[i];
        out_data[OUT_RA_GRND][0] += cell.Ra_grnd[i];
        out_data[OUT_RA_SUB][0] += cell.Ra_sub[i];
    }
    out_data[OUT_RA_EVAP][0] += cell.Ra_evap;
    if (HasVeg) {    
        out_data[OUT_RA_LEAF][0] += cell.Ra_leaf;
        out_data[OUT_RA_STEM][0] += cell.Ra_stem;
    }

    /** record nodes moistures **/
    for (i = 0; i < cell.Nsoil; i++) {
        out_data[OUT_SOIL_LIQ][i] += cell.liq[i] * Cv;
        out_data[OUT_SOIL_ICE][i] += cell.ice[i] * Cv;
        out_data[OUT_SOIL_MOIST][i] += cell.moist[i] * Cv;
        out_data[OUT_MATRIC][i] += cell.matric[i] * Cv;

    }
    out_data[OUT_ROOTMOIST][0] += cell.rootmoist * Cv;

    /** record water table position **/
    out_data[OUT_ZWT][0] += cell.zwt * Cv;

    /*****************************
       Record Snow Pack Variables
    *****************************/
    /** record snow water equivalence **/
    out_data[OUT_SWE][0] += snow.swq * Cv;
    /** record snowpack depth **/
    out_data[OUT_SNOW_DEPTH][0] += snow.snow_depth * Cv;
    /** record snowpack temperature water and ice content **/
    Nsnow = snow.Nsnow;
    for (i = 0; i < Nsnow; i++) {
        out_data[OUT_SNOW_PACK_TEMP][i] += snow.pack_T[i] * Cv;
        out_data[OUT_SNOW_PACK_ICE][i] += snow.pack_ice[i] * Cv;
        out_data[OUT_SNOW_PACK_LIQ][i] += snow.pack_liq[i] * Cv;
        out_data[OUT_PACK_OUTFLOW][i] += snow.pack_outflow[i] * Cv;
        out_data[OUT_SNOW_ICEFRAC][i] += snow.theta_ice[i] * Cv;
        out_data[OUT_SNOW_LIQFRAC][i] += snow.theta_liq[i] * Cv;
        out_data[OUT_SNOW_POROSITY][i] += snow.porosity[i] * Cv;
        /** record snow density **/
        out_data[OUT_SNOW_DENSITY][i] += snow.density[i] * Cv;
        /** record snowpack melt **/
        out_data[OUT_SNOW_MELT][i] += snow.pack_melt[i] * Cv;
        /** record snow surface freezing rate [mm/s] */
        out_data[OUT_SNOW_FRZE][i] += snow.pack_frze[i] * Cv;
    }

    /** record snow surface evaporation **/
    out_data[OUT_SNOW_EVAP][0] += snow.snow_evap * Cv * dt_sec;
    /** record snowpack combination **/
    out_data[OUT_SNOW_COMB][0] += snow.pack_comb * Cv;
    /** record snow surface frost [mm] */
    out_data[OUT_SNOW_FROST][0] += snow.snow_frost * Cv * dt_sec;
    /** record snow surface dew [mm] */
    out_data[OUT_SNOW_DEW][0] += snow.snow_dew * Cv * dt_sec;
    /** record snow surface sublimation [mm] */
    out_data[OUT_SNOW_SUBLIM][0] += snow.snow_sublim * Cv * dt_sec;
    /** record glacier snow excess flow **/
    out_data[OUT_GLAC_EXCESS][0] += snow.glac_excess * Cv;
    /** record snow cover fraction **/
    out_data[OUT_SNOW_COVER][0] += snow.coverage * Cv;
    /** record snow depth increasing rate **/
    out_data[OUT_DELDEPTH][0] += snow.delta_depth * Cv;
    /** record NEW snow density **/
    out_data[OUT_NEW_DENSITY][0] += snow.new_snow_density * Cv;
    /** record SnowAge **/
    out_data[OUT_SNOW_AGE][0] += snow.snowage * Cv;
    /** outflow of liquid water from the snowpack bottom (m/s) */
    out_data[OUT_SNOW_OUTFLOW][0] += snow.snow_outflow * Cv;

    // Glacier Water Balance Terms
    if (HasGlac) {
        out_data[OUT_GLAC_INFLOW][0] += cell.soil_inflow * Cv;
        out_data[OUT_GLAC_OUTFLOW][0] += (cell.runoff + cell.baseflow) * Cv;
    }
}

/******************************************************************************
 * @brief    This routine collects energy balance terms.
 *****************************************************************************/
void
collect_eb_terms(energy_bal_struct energy,
                 snow_data_struct  snow,
                 cell_data_struct  cell,
                 double            Cv,
                 bool              HasVeg,
                 bool              HasGlac,
                 int               band,
                 double          **out_data)
{
    extern option_struct options;
    size_t i;

    /**************************************
       Record Frozen Grnd and Canopy flag
    **************************************/
    out_data[OUT_TRND_FBFLAG][0] += energy.FrozenGrnd;
    if (HasVeg) {
        out_data[OUT_TCAN_FBFLAG][0] += energy.FrozenOver;
    }

    /**********************************
       Record Energy Balance Variables
    **********************************/
    /** record landcover temperature **/
    if (HasVeg) {
        // landcover is vegetation
        out_data[OUT_VEGT][0] += energy.Tfoliage * Cv;
        out_data[OUT_VEGTAIR][0] += energy.Tcanopy * Cv;
    }
    // landcover is bare soil
    out_data[OUT_BARESOILT][0] += energy.Tgrnd * Cv;
    /** record mean surface temperature [C] **/
    out_data[OUT_SURF_TEMP][0] += energy.Tsurf * Cv;

    /** record NODES temperatures **/
    for (i = 0; i < cell.Nsoil; i++) {
        out_data[OUT_SOIL_TEMP][i] += cell.soil_T[i] * Cv;
    }
    /** record advective flux from prec **/
    out_data[OUT_ADVECTION][0] += energy.advection * Cv;
    out_data[OUT_ADVECTGRND][0] += energy.AdvectGrnd * Cv;
    if (HasVeg) {
        out_data[OUT_ADVECTSUB][0] += energy.AdvectSub * Cv;
        out_data[OUT_ADVECTOVER][0] += energy.AdvectOver * Cv;
    }
    /** record net shortwave radiation **/
    out_data[OUT_SWNET][0] += energy.shortwave * Cv;
    /** record longwave radiation flux **/
    out_data[OUT_LWNET][0] += energy.longwave * Cv;
    /** record latent heat flux **/
    out_data[OUT_LATENT][0] += energy.latent * Cv;
    /** record sensible heat flux **/
    out_data[OUT_SENSIBLE][0] += energy.sensible * Cv;
    /** record ground heat flux **/
    out_data[OUT_GRND_FLUX][0] += energy.grnd_flux * Cv;
    //out_data[OUT_GRND_GRND][0] += energy.GroundGrnd * Cv;
    // 冰川变量
    if (HasGlac) {
        out_data[OUT_H2OSFC_T][0] += energy.Tgrnd * Cv;
    }
    /**********************************
       Record hru-Specific Variables
    **********************************/

    /** record hru snow water equivalent **/
    out_data[OUT_SWE_BAND][band] += snow.swq * Cv; // (mm/H2O)

    /** record band snowpack depth **/
    out_data[OUT_SNOW_DEPTH_BAND][band] += snow.snow_depth * Cv; // (M)

    /** record band snow coverage **/
    out_data[OUT_SNOW_COVER_BAND][band] += snow.coverage * Cv;

    /** record band advection **/
    out_data[OUT_ADVECTION_BAND][band] += energy.advection * Cv;

    /** record pack layer temperature **/
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
                     force_data_struct *force,
                     veg_con_struct    *veg_con,
                     double           **out_data,
                     timer_struct      *timer)
{
    // Calling put data will populate the save data storage terms
    put_data(all_vars, force, veg_con,
             out_data, timer);

    zero_output_list(out_data);
}