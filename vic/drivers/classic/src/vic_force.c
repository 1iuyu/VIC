/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes atmospheric variables for both the model time step,
 * and the time step used by the snow algorithm (if different).
 *****************************************************************************/

#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    Initialize atmospheric variables for the model and snow time steps.
 *****************************************************************************/
void
vic_force(force_data_struct *force,
          dmy_struct        *dmy,
          FILE             **infile,
          veg_con_struct    *veg_con,
          veg_hist_struct  **veg_hist,
          soil_con_struct   *soil_con)
{
    extern option_struct       options;
    extern param_set_struct    param_set;
    extern global_param_struct global_param;
    extern parameters_struct   param;
    extern size_t              NR, NF;

    size_t                     i;
    size_t                     j;
    size_t                     v;
    size_t                     rec;
    size_t                     uidx;
    double                     t_offset;
    double                   **forcing_data;
    double                  ***veg_hist_data;
    double                     avgJulyAirTemp;
    double                    *Tfactor;
    bool                      *AboveTreeLine;

    /*******************************
       Check that required inputs were supplied
    *******************************/

    if (!param_set.TYPE[AIR_TEMP].SUPPLIED) {
        log_err("Air temperature must be supplied as a forcing");
    }
    if (!param_set.TYPE[PREC].SUPPLIED && !param_set.TYPE[RAINF].SUPPLIED &&
                                          !param_set.TYPE[SNOWF].SUPPLIED) {
        log_err("PREC or RAINF and SNOWF must be supplied as a forcing");
    }
    if (!param_set.TYPE[SWDOWN].SUPPLIED) {
        log_err("Downward shortwave radiation must be supplied as a forcing");
    }
    if (!param_set.TYPE[LWDOWN].SUPPLIED) {
        log_err("Downward longwave radiation must be supplied as a forcing");
    }
    if (!param_set.TYPE[PRESSURE].SUPPLIED) {
        log_err("Atmospheric pressure must be supplied as a forcing");
    }
    if (!param_set.TYPE[VP].SUPPLIED) {
        log_err("Vapor ressure must be supplied as a forcing");
    }
    if (!param_set.TYPE[WIND].SUPPLIED) {
        log_err("Wind speed must be supplied as a forcing");
    }

    /*******************************
       Miscellaneous initialization
    *******************************/

    /* Assign local copies of some variables */
    avgJulyAirTemp = soil_con->avgJulyAirTemp;
    Tfactor = soil_con->Tfactor;
    AboveTreeLine = soil_con->AboveTreeLine;

    /* Assign N_ELEM for veg-dependent forcings */
    if (param_set.TYPE[ALBEDO].SUPPLIED) {
        param_set.TYPE[ALBEDO].N_ELEM = veg_con[0].vegetat_type_num;
    }
    if (param_set.TYPE[LAI].SUPPLIED) {
        param_set.TYPE[LAI].N_ELEM = veg_con[0].vegetat_type_num;
    }
    if (param_set.TYPE[FCANOPY].SUPPLIED) {
        param_set.TYPE[FCANOPY].N_ELEM = veg_con[0].vegetat_type_num;
    }
    if (param_set.TYPE[PREC].SUPPLIED && !param_set.TYPE[RAINF].SUPPLIED && 
                                         !param_set.TYPE[SNOWF].SUPPLIED) {
        hasRAINF_SNOWF = false;
    }
    else {
        hasRAINF_SNOWF = true;
    }

    /*******************************
       read in meteorological data
    *******************************/

    forcing_data = read_forcing_data(infile, global_param, &veg_hist_data);

    log_info("Read meteorological forcing file");

    /****************************************************
       Variables in the atmos_data structure
    ****************************************************/

    t_offset = Tfactor[0];
    for (i = 1; i < options.SNOW_BAND; i++) {
        if (Tfactor[i] < t_offset) {
            t_offset = Tfactor[i];
        }
    }

    for (rec = 0; rec < global_param.nrecs; rec++) {
        for (i = 0; i < NF; i++) {
            uidx = rec * NF + i;
            // temperature in Celsius
            force[rec].air_temp[i] = forcing_data[AIR_TEMP][uidx];
            if (hasRAINF_SNOWF) {
                // rain in mm/period
                force[rec].rainf[i] = forcing_data[RAINF][uidx];
                // snow in mm/period
                force[rec].snowf[i] = forcing_data[SNOWF][uidx];
                // snow flag
                force[rec].snowflag[i] = will_it_snow(&(force[rec].snowf[i]), 1);
            }
            else {
                /** Calculate Fraction of Precipitation that falls as Rain **/
                force[rec].prec[i] = forcing_data[PREC][uidx];
                rainonly = calc_rainonly(force[rec].air_temp[i], force[rec].prec[i],
                                         soil_con->MAX_SNOW_TEMP, soil_con->MIN_RAIN_TEMP);
                force[rec].rainf[i] = rainonly;
                force[rec].snowf[i] = force[rec].prec[i] - rainonly;
                // snow flag
                force[rec].snowflag[i] = will_it_snow(&(force[rec].snowf[i]), 1);
            }
            // downward shortwave in W/m2
            force[rec].shortwave[i] = forcing_data[SWDOWN][uidx];
            // downward longwave in W/m2
            force[rec].longwave[i] = forcing_data[LWDOWN][uidx];
            // pressure in Pa
            force[rec].pressure[i] = forcing_data[PRESSURE][uidx] * PA_PER_KPA;
            // vapor pressure in Pa
            force[rec].vp[i] = forcing_data[VP][uidx] * PA_PER_KPA;
            // vapor pressure deficit in Pa
            force[rec].vpd[i] = svp(force[rec].air_temp[i]) - force[rec].vp[i];
            if (force[rec].vpd[i] < 0) {
                force[rec].vpd[i] = 0;
                force[rec].vp[i] = svp(force[rec].air_temp[i]);
            }
            // air density in kg/m3
            force[rec].density[i] = air_density(force[rec].air_temp[i],
                                                force[rec].pressure[i]);
            // wind speed in m/s
            force[rec].wind[i] = forcing_data[WIND][uidx];
            // Optional inputs
            if (options.LAKES) {
                // Channel inflow from upstream (into lake)
                if (param_set.TYPE[CHANNEL_IN].SUPPLIED) {
                    force[rec].channel_in[i] = forcing_data[CHANNEL_IN][uidx];
                }
                else {
                    force[rec].channel_in[i] = 0;
                }
            }
            if (options.CARBON) {
                // Atmospheric CO2 concentration
                force[rec].Catm[i] = forcing_data[CATM][uidx];
                // Fraction of shortwave that is direct
                force[rec].fdir[i] = forcing_data[FDIR][uidx];
                // photosynthetically active radiation
                force[rec].par[i] = forcing_data[PAR][uidx];
                // Cosine of solar zenith angle
                force[rec].coszen[i] = compute_coszen(soil_con->lat,
                                                      soil_con->lng,
                                                      soil_con->time_zone_lng,
                                                      dmy[rec].day_in_year,
                                                      dmy[rec].dayseconds);
            }
        }
        if (NF > 1) {
            force[rec].air_temp[NR] = average(force[rec].air_temp, NF);
            // For precipitation put total
            force[rec].prec[NR] = average(force[rec].prec, NF) * NF;
            force[rec].rainf[NR] = average(force[rec].rainf, NF) * NF;
            force[rec].snowf[NR] = average(force[rec].snowf, NF) * NF;
            force[rec].shortwave[NR] = average(force[rec].shortwave, NF);
            force[rec].longwave[NR] = average(force[rec].longwave, NF);
            force[rec].pressure[NR] = average(force[rec].pressure, NF);
            force[rec].vp[NR] = average(force[rec].vp, NF);
            force[rec].vpd[NR] = average(force[rec].vpd, NF);
            force[rec].density[NR] = average(force[rec].density, NF);
            force[rec].wind[NR] = average(force[rec].wind, NF);
            force[rec].snowflag[NR] = false;
            for (i = 0; i < NF; i++) {
                if (force[rec].snowflag[i] == true) {
                    force[rec].snowflag[NR] = true;
                }
            }
            if (options.LAKES) {
                force[rec].channel_in[NR] =
                    average(force[rec].channel_in, NF) * NF;
            }
            if (options.CARBON) {
                force[rec].Catm[NR] = average(force[rec].Catm, NF);
                force[rec].fdir[NR] = average(force[rec].fdir, NF);
                force[rec].par[NR] = average(force[rec].par, NF);
                // for coszen, use value at noon
                force[rec].coszen[NR] = compute_coszen(soil_con->lat,
                                                       soil_con->lng,
                                                       soil_con->time_zone_lng,
                                                       dmy[rec].day_in_year,
                                                       SEC_PER_DAY / 2);
            }
        }
    }

    /****************************************************
       Variables in the veg_hist structure
    ****************************************************/

    /* First, assign default climatology */
    for (rec = 0; rec < global_param.nrecs; rec++) {
        for (v = 0; v <= veg_con[0].vegetat_type_num; v++) {
            for (i = 0; i < NF; i++) {
                veg_hist[rec][v].albedo[i] =
                    veg_con[v].albedo[dmy[rec].month - 1];
                veg_hist[rec][v].displacement[i] =
                    veg_con[v].displacement[dmy[rec].month - 1];
                veg_hist[rec][v].fcanopy[i] =
                    veg_con[v].fcanopy[dmy[rec].month - 1];
                veg_hist[rec][v].LAI[i] =
                    veg_con[v].LAI[dmy[rec].month - 1];
                veg_hist[rec][v].roughness[i] =
                    veg_con[v].roughness[dmy[rec].month - 1];
            }
        }
    }

    /* Next, overwrite with veg_hist values, validate, and average */
    for (rec = 0; rec < global_param.nrecs; rec++) {
        for (v = 0; v <= veg_con[0].vegetat_type_num; v++) {
            for (i = 0; i < NF; i++) {
                uidx = rec * NF + i;
                if (param_set.TYPE[ALBEDO].SUPPLIED &&
                    options.ALB_SRC == FROM_VEGHIST) {
                    if (veg_hist_data[ALBEDO][v][uidx] != NODATA_VH) {
                        veg_hist[rec][v].albedo[i] =
                            veg_hist_data[ALBEDO][v][uidx];
                    }
                }
                if (param_set.TYPE[LAI].SUPPLIED &&
                    options.LAI_SRC == FROM_VEGHIST) {
                    if (veg_hist_data[LAI][v][uidx] != NODATA_VH) {
                        veg_hist[rec][v].LAI[i] =
                            veg_hist_data[LAI][v][uidx];
                    }
                }
                if (param_set.TYPE[FCANOPY].SUPPLIED &&
                    options.FCAN_SRC == FROM_VEGHIST) {
                    if (veg_hist_data[FCANOPY][v][uidx] != NODATA_VH) {
                        veg_hist[rec][v].fcanopy[i] =
                            veg_hist_data[FCANOPY][v][uidx];
                    }
                }
                // Checks on fcanopy and LAI
                if (veg_hist[rec][v].fcanopy[i] < MIN_FCANOPY ||
                    veg_hist[rec][v].LAI[i] == 0) {
                    log_warn(
                        "rec %zu, veg %zu substep %zu fcanopy %f < "
                        "minimum of %f; setting = 0", rec, v, i,
                        veg_hist[rec][v].fcanopy[i], MIN_FCANOPY);
                    veg_hist[rec][v].fcanopy[i] = 0;
                    veg_hist[rec][v].LAI[i] = 0;
                }
            }
            if (NF > 1) {
                veg_hist[rec][v].albedo[NR] = average(veg_hist[rec][v].albedo,
                                                      NF);
                veg_hist[rec][v].displacement[NR] = average(
                    veg_hist[rec][v].displacement, NF);
                veg_hist[rec][v].fcanopy[NR] = average(
                    veg_hist[rec][v].fcanopy, NF);
                veg_hist[rec][v].LAI[NR] = average(veg_hist[rec][v].LAI, NF);
                veg_hist[rec][v].roughness[NR] = average(
                    veg_hist[rec][v].roughness, NF);
            }
        }
    }

    /****************************************************
       Free forcing_data and veg_hist_data structures
    ****************************************************/

    for (i = 0; i < N_FORCING_TYPES; i++) {
        if (param_set.TYPE[i].SUPPLIED) {
            if (i != ALBEDO && i != LAI && i != FCANOPY) {
                free(forcing_data[i]);
            }
            else {
                for (j = 0; j < param_set.TYPE[i].N_ELEM; j++) {
                    free(veg_hist_data[i][j]);
                }
                free(veg_hist_data[i]);
            }
        }
    }
    free(forcing_data);
    free(veg_hist_data);

    /****************************************************
       Compute treeline based on July average temperature
    ****************************************************/

    if (options.COMPUTE_TREELINE) {
        if (!(options.JULY_TAVG_SUPPLIED && avgJulyAirTemp == -999)) {
            compute_treeline(force, dmy, avgJulyAirTemp, Tfactor,
                             AboveTreeLine);
        }
    }
}
