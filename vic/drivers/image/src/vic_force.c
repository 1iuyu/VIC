/******************************************************************************
 * @section DESCRIPTION
 *
 * Read atmospheric forcing data.
 *****************************************************************************/

#include "vic_driver_image.h"

/******************************************************************************
 * @brief    Read atmospheric forcing data.
 *****************************************************************************/
void
vic_force(void)
{
    extern size_t              NF;
    extern size_t              NR;
    extern size_t              current;
    extern int                 mpi_rank;
    extern force_data_struct  *force;
    extern dmy_struct         *dmy;
    extern domain_struct       global_domain;
    extern domain_struct       local_domain;
    extern filenames_struct    filenames;
    extern global_param_struct global_param;
    extern option_struct       options;
    extern soil_con_struct    *soil_con;
    extern veg_con_struct    **veg_con;
    extern veg_hist_struct   **veg_hist;
    extern param_set_struct    param_set;

    double                    *dvar = NULL;
    size_t                     i, j, v;
    int                        status;
    size_t                     d3count[3];
    size_t                     d3start[3];
    size_t                     d4count[4];
    size_t                     d4start[4];

    // allocate memory for variables to be read
    dvar = malloc(local_domain.ncells_active * sizeof(*dvar));
    check_alloc_status(dvar, "Memory allocation error.");

    // only the time slice changes for the met file reads. The rest is constant
    d3start[1] = 0;
    d3start[2] = 0;
    d3count[0] = 1;
    d3count[1] = global_domain.n_ny;
    d3count[2] = global_domain.n_nx;

    /*******************************
      Miscellaneous initialization
    *******************************/
    int has_prec = param_set.TYPE[PREC].SUPPLIED;
    int has_rain = param_set.TYPE[RAINF].SUPPLIED;
    int has_snow = param_set.TYPE[SNOWF].SUPPLIED;

    log_info("Read meteorological forcing file");

    // Air temperature: air_temp
    for (j = 0; j < NF; j++) {
        d3start[0] = global_param.forceskip[0] + global_param.forceoffset[0] +
                     j;
        get_scatter_nc_field_double(&(filenames.forcing[0]),
                                    param_set.TYPE[AIR_TEMP].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            force[i].air_temp[j] = (double) dvar[i];
        }
    }
    if (has_prec && !has_rain && !has_snow) {
        // Precipitation: prec [mm/s]
        for (j = 0; j < NF; j++) {
            d3start[0] = global_param.forceskip[0] + global_param.forceoffset[0] +
                        j;
            get_scatter_nc_field_double(&(filenames.forcing[0]),
                                        param_set.TYPE[PREC].varname,
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                force[i].prec[j] = (double) dvar[i];
                force[i].rainf[j] = 0.;
                force[i].snowf[j] = 0.;
            }
        }
    }
    else {
        if (has_rain && has_snow) {
            // Snowfall: snowf [mm/s]
            for (j = 0; j < NF; j++) {
                d3start[0] = global_param.forceskip[0] + global_param.forceoffset[0] +
                            j;
                get_scatter_nc_field_double(&(filenames.forcing[0]),
                                            param_set.TYPE[SNOWF].varname,
                                            d3start, d3count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    force[i].snowf[j] = (double) dvar[i];
                }
            }
            // Rainfall: rainf [mm/s]
            for (j = 0; j < NF; j++) {
                d3start[0] = global_param.forceskip[0] + global_param.forceoffset[0] +
                            j;
                get_scatter_nc_field_double(&(filenames.forcing[0]),
                                            param_set.TYPE[RAINF].varname,
                                            d3start, d3count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    force[i].rainf[j] = (double) dvar[i];
                }
            }
        }
        else {
            log_err("Either prec must be supplied as a forcing, or both rain and snow must be supplied");
        }
    }

    // downward shortwave in W/m2
    for (j = 0; j < NF; j++) {
        d3start[0] = global_param.forceskip[0] + global_param.forceoffset[0] +
                     j;
        get_scatter_nc_field_double(&(filenames.forcing[0]),
                                    param_set.TYPE[SWDOWN].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            force[i].shortwave[j] = (double) dvar[i];
        }
    }

    // downward longwave in W/m2
    for (j = 0; j < NF; j++) {
        d3start[0] = global_param.forceskip[0] + global_param.forceoffset[0] +
                     j;
        get_scatter_nc_field_double(&(filenames.forcing[0]),
                                    param_set.TYPE[LWDOWN].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            force[i].longwave[j] = (double) dvar[i];
        }
    }

    // Wind speed: wind [m/s]
    for (j = 0; j < NF; j++) {
        d3start[0] = global_param.forceskip[0] + global_param.forceoffset[0] +
                     j;
        get_scatter_nc_field_double(&(filenames.forcing[0]),
                                    param_set.TYPE[WIND].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            force[i].wind[j] = (double) dvar[i];
        }
    }

    // rel_humid [%]
    for (j = 0; j < NF; j++) {
        d3start[0] = global_param.forceskip[0] + global_param.forceoffset[0] +
                     j;
        get_scatter_nc_field_double(&(filenames.forcing[0]),
                                    param_set.TYPE[REL_HUMID].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            force[i].rel_humid[j] = (double) dvar[i];
        }
    }

    // Pressure: pressure [Pa]
    for (j = 0; j < NF; j++) {
        d3start[0] = global_param.forceskip[0] + global_param.forceoffset[0] +
                     j;
        get_scatter_nc_field_double(&(filenames.forcing[0]),
                                    param_set.TYPE[PRESSURE].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            force[i].pressure[j] = (double) dvar[i];
        }
    }
    // specific humidity [kg/kg]: Qair
    for (j = 0; j < NF; j++) {
        d3start[0] = global_param.forceskip[0] + 
                     global_param.forceoffset[0] + j;
        get_scatter_nc_field_double(&(filenames.forcing[0]),
                                    param_set.TYPE[QAIR].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            force[i].Qair[j] = (double) dvar[i];
        }
    }

    if (options.CARBON) {
        // Atmospheric CO2 mixing ratio
        for (j = 0; j < NF; j++) {
            d3start[0] = global_param.forceskip[0] +
                         global_param.forceoffset[0] + j;
            get_scatter_nc_field_double(&(filenames.forcing[0]),
                                        param_set.TYPE[CATM].varname,
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                force[i].Catm[j] = (double) dvar[i];
            }
        }
        // Fraction of shortwave that is direct
        for (j = 0; j < NF; j++) {
            d3start[0] = global_param.forceskip[0] +
                         global_param.forceoffset[0] + j;
            get_scatter_nc_field_double(&(filenames.forcing[0]),
                                        param_set.TYPE[FDIR].varname,
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                force[i].fdir[j] = (double) dvar[i];
            }
        }
        // Photosynthetically active radiation
        for (j = 0; j < NF; j++) {
            d3start[0] = global_param.forceskip[0] +
                         global_param.forceoffset[0] + j;
            get_scatter_nc_field_double(&(filenames.forcing[0]),
                                        param_set.TYPE[PAR].varname,
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                force[i].par[j] = (double) dvar[i];
            }
        }
    }
    // Channel inflow for time step (m/s)
    if (options.ROUT) {
        for (j = 0; j < NF; j++) {
            d3start[0] = global_param.forceskip[0] +
                         global_param.forceoffset[0] + j;
            get_scatter_nc_field_double(&(filenames.forcing[0]),
                                        param_set.TYPE[CHANNEL_IN].varname,
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                force[i].channel_in[j] = (double) dvar[i];
            }
        }
    }

    if (mpi_rank == VIC_MPI_ROOT) {
        // Close forcing file if it is the last time step
        if (current == global_param.nrecs - 1) {
            status = nc_close(filenames.forcing[0].nc_id);
            check_nc_status(status, "Error closing %s",
                            filenames.forcing[0].nc_filename);
        }
    }

    // Update the offset counter
    global_param.forceoffset[0] += NF;

    // Initialize the veg_hist structure with the current climatological
    // vegetation parameters.  This may be overwritten with the historical
    // forcing time series.
    for (i = 0; i < local_domain.ncells_active; i++) {
        for (v = 0; v < local_domain.locations[i].nveg; v++) {
            for (j = 0; j < NF; j++) {
                veg_hist[i][v].fcanopy[j] =
                    veg_con[i][v].fcanopy[dmy[current].month - 1];
                veg_hist[i][v].LAI[j] =
                    veg_con[i][v].LAI[dmy[current].month - 1];
                veg_hist[i][v].SAI[j] =
                    veg_con[i][v].SAI[dmy[current].month - 1];
            }
        }
    }

    // Read veg_hist file
    if (options.LAI_SRC == FROM_VEGHIST ||
        options.FCAN_SRC == FROM_VEGHIST) {
        // only the time slice changes for the met file reads. The rest is constant
        d4start[2] = 0;
        d4start[3] = 0;
        d4count[0] = 1;
        d4count[1] = 1;
        d4count[2] = global_domain.n_ny;
        d4count[3] = global_domain.n_nx;

        // Leaf Area Index: LAI
        if (options.LAI_SRC == FROM_VEGHIST) {
            for (j = 0; j < NF; j++) {
                d4start[0] = global_param.forceskip[1] +
                             global_param.forceoffset[1] + j;
                for (v = 0; v < MAX_HRUS; v++) {
                    d4start[1] = v;
                    get_scatter_nc_field_double(&(filenames.forcing[1]), "LAI",
                                                d4start, d4count, dvar);
                    for (i = 0; i < local_domain.ncells_active; i++) {
                        if (v < local_domain.locations[i].nveg) {
                            veg_hist[i][v].LAI[j] = (double) dvar[i];
                        }
                    }
                }
            }
        }

        // Stem Area Index: SAI
        if (options.LAI_SRC == FROM_VEGHIST) {
            for (j = 0; j < NF; j++) {
                d4start[0] = global_param.forceskip[1] +
                             global_param.forceoffset[1] + j;
                for (v = 0; v < MAX_HRUS; v++) {
                    d4start[1] = v;
                    get_scatter_nc_field_double(&(filenames.forcing[1]), "SAI",
                                                d4start, d4count, dvar);
                    for (i = 0; i < local_domain.ncells_active; i++) {
                        if (v < local_domain.locations[i].nveg) {
                            veg_hist[i][v].SAI[j] = (double) dvar[i];
                        }
                    }
                }
            }
        }

        // Partial veg cover fraction: fcanopy
        if (options.FCAN_SRC == FROM_VEGHIST) {
            for (j = 0; j < NF; j++) {
                d4start[0] = global_param.forceskip[1] +
                             global_param.forceoffset[1] + j;
                for (v = 0; v < MAX_HRUS; v++) {
                    d4start[1] = v;
                    get_scatter_nc_field_double(&(filenames.forcing[1]),
                                                "fcanopy", d4start, d4count,
                                                dvar);
                    for (i = 0; i < local_domain.ncells_active; i++) {
                        if (v < local_domain.locations[i].nveg) {
                            veg_hist[i][v].fcanopy[j] = (double) dvar[i];
                        }
                    }
                }
            }
        }

        if (mpi_rank == VIC_MPI_ROOT) {
            // Close forcing file if it is the last time step
            if (current == global_param.nrecs - 1) {
                status = nc_close(filenames.forcing[1].nc_id);
                check_nc_status(status, "Error closing %s",
                                filenames.forcing[1].nc_filename);
            }
        }

        // Update the offset counter
        global_param.forceoffset[1] += NF;
    }

    // Convert forcings into what we need and calculate missing ones
    for (i = 0; i < local_domain.ncells_active; i++) {
        for (j = 0; j < NF; j++) {
            // pressure in Pa
            force[i].pressure[j] *= PA_PER_KPA;
            // vapor pressure in Pa
            force[i].vp[j] = q_to_vp(force[i].Qair[j], 
                                force[i].pressure[j]);

            // air density in kg/m3
            force[i].density[j] = air_density(force[i].air_temp[j],
                                              force[i].pressure[j],
                                              force[i].vp[j]);
            // Cosine of solar zenith angle: coszen
            compute_coszen(
                    local_domain.locations[i].latitude,
                    local_domain.locations[i].longitude,
                    soil_con[i].time_zone_lng, 
                    dmy[current].day_in_year,
                    dmy[current].dayseconds,
                    &force[i].coszen[j],
                    &force[i].daylen[j]);
            // atmospheric potential temperature (K)
            force[i].theta_pot[j] = compute_theta(force[i].air_temp[j],
                                                  force[i].pressure[j]);
            // atmospheric virtual potential temperature (K)
            force[i].theta_v[j] = compute_theta_v(force[i].Qair[j],
                                                  force[i].theta_pot[j]);
        }
    }

    // Checks on fcanopy and LAI
    for (i = 0; i < local_domain.ncells_active; i++) {
        for (v = 0; v < local_domain.locations[i].nveg; v++) {
            for (j = 0; j < NF; j++) {
                if (veg_hist[i][v].fcanopy[j] < MIN_FCANOPY ||
                    veg_hist[i][v].LAI[j] == 0) {
                    if (current == 0 || options.FCAN_SRC == FROM_VEGHIST) {
                        // Only issue this warning once
                        log_warn(
                            "cell %zu, veg %d substep %zu either fcanopy "
                            "%f < minimum of %f or LAI %f == 0; setting "
                            "both LAI and fcanopy to 0", i, v, j,
                            veg_hist[i][v].fcanopy[j], MIN_FCANOPY,
                            veg_hist[i][v].LAI[j]);
                    }
                    veg_hist[i][v].fcanopy[j] = 0;
                    veg_hist[i][v].LAI[j] = 0;
                }
            }
        }
    }

    // Put average value in NR field
    if (NF > 1) {
        for (i = 0; i < local_domain.ncells_active; i++) {
            force[i].air_temp[NR] = average(force[i].air_temp, NF);
            force[i].prec[NR] = average(force[i].prec, NF);
            force[i].rainf[NR] = average(force[i].rainf, NF);
            force[i].snowf[NR] = average(force[i].snowf, NF);
            force[i].shortwave[NR] = average(force[i].shortwave, NF);
            force[i].longwave[NR] = average(force[i].longwave, NF);
            force[i].pressure[NR] = average(force[i].pressure, NF);
            force[i].wind[NR] = average(force[i].wind, NF);
            force[i].vp[NR] = average(force[i].vp, NF);
            force[i].Qair[NR] = average(force[i].Qair, NF);
            force[i].rel_humid[NR] = average(force[i].rel_humid, NF);
            force[i].density[NR] = air_density(force[i].air_temp[NR],
                                            force[i].pressure[NR],
                                            force[i].vp[NR]);
            // for coszen, use value at noon
            compute_coszen(local_domain.locations[i].latitude,
                           local_domain.locations[i].longitude,
                           soil_con[i].time_zone_lng,
                           dmy[current].day_in_year,
                           SEC_PER_DAY / 2,
                           &force[i].coszen[NR],
                           &force[i].daylen[NR]);

            for (v = 0; v < local_domain.locations[i].nveg; v++) {
                    // not the correct way to calculate average albedo in general,
                    // but leave for now (it's correct if albedo is constant over
                    // the model step)
                    veg_hist[i][v].fcanopy[NR] = average(
                        veg_hist[i][v].fcanopy, NF);
                    veg_hist[i][v].LAI[NR] = average(veg_hist[i][v].LAI, NF);
                    veg_hist[i][v].SAI[NR] = average(veg_hist[i][v].SAI, NF);
            }

            // Optional inputs
            if (options.CARBON) {
                force[i].Catm[NR] = average(force[i].Catm, NF);
                force[i].fdir[NR] = average(force[i].fdir, NF);
                force[i].par[NR] = average(force[i].par, NF);
            }
            if (options.ROUT) {
                force[i].channel_in[NR] = average(force[i].channel_in, NF);
            }
        }
    }

    // cleanup
    free(dvar);
}

/******************************************************************************
 * @brief    Determine timestep and start year, month, day, and seconds of forcing files
 *****************************************************************************/
void
get_forcing_file_info(param_set_struct *param_set,
                      size_t            file_num)
{
    extern global_param_struct global_param;
    extern filenames_struct    filenames;

    double                     nc_times[2];
    double                     nc_time_origin;
    size_t                     start = 0;
    size_t                     count = 2;
    char                      *nc_unit_chars = NULL;
    char                      *calendar_char = NULL;
    unsigned short int         time_units;
    unsigned short int         calendar;
    dmy_struct                 nc_origin_dmy;
    dmy_struct                 nc_start_dmy;

    // read time info from netcdf file
    get_nc_field_double(&(filenames.forcing[file_num]), "time", &start, &count,
                        nc_times);
    get_nc_var_attr(&(filenames.forcing[file_num]), "time", "units",
                    &nc_unit_chars);
    get_nc_var_attr(&(filenames.forcing[file_num]), "time", "calendar",
                    &calendar_char);

    // parse the calendar string and check to make sure it matches the global clock
    calendar = str_to_calendar(calendar_char);

    // parse the time units
    parse_nc_time_units(nc_unit_chars, &time_units, &nc_origin_dmy);

    // Get date/time of the first entry in the forcing file.
    nc_time_origin =
        date2num(0., &nc_origin_dmy, 0., calendar, TIME_UNITS_DAYS);
    num2date(nc_time_origin, nc_times[0], 0., calendar, time_units,
             &nc_start_dmy);

    // Assign file start date/time
    global_param.forceyear[file_num] = nc_start_dmy.year;
    global_param.forcemonth[file_num] = nc_start_dmy.month;
    global_param.forceday[file_num] = nc_start_dmy.day;
    global_param.forcesec[file_num] = nc_start_dmy.dayseconds;

    // calculate timestep in forcing file
    if (time_units == TIME_UNITS_DAYS) {
        param_set->force_steps_per_day[file_num] =
            (size_t) nearbyint(1. / (nc_times[1] - nc_times[0]));
    }
    else if (time_units == TIME_UNITS_HOURS) {
        param_set->force_steps_per_day[file_num] =
            (size_t) nearbyint(HOURS_PER_DAY / (nc_times[1] - nc_times[0]));
    }
    else if (time_units == TIME_UNITS_MINUTES) {
        param_set->force_steps_per_day[file_num] =
            (size_t) nearbyint(MIN_PER_DAY / (nc_times[1] - nc_times[0]));
    }
    else if (time_units == TIME_UNITS_SECONDS) {
        param_set->force_steps_per_day[file_num] =
            (size_t) nearbyint(SEC_PER_DAY / (nc_times[1] - nc_times[0]));
    }

    // check that this forcing file will work
    if (param_set->force_steps_per_day[file_num] >
        global_param.model_steps_per_day) {
        log_err("Forcing file timestep must match the snow model timestep.  "
                "Snow model timesteps per day is set to %zu and the forcing "
                "file timestep is set to %zu",
                global_param.model_steps_per_day,
                param_set->force_steps_per_day[file_num])
    }
    if (calendar != global_param.calendar) {
        log_err("Calendar in forcing file (%s) does not match the calendar of "
                "VIC's clock", calendar_char);
    }

    // Free attribute character arrays
    free(nc_unit_chars);
    free(calendar_char);
}
