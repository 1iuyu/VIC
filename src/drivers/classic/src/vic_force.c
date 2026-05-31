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
vic_force(grid_cell_struct *grid,
          dmy_struct       *dmy,
          FILE            **infile)
{
    extern option_struct       options;
    extern param_set_struct    param_set;
    extern global_param_struct global_param;
    extern size_t              NR, NF;
    extern domain_struct       local_domain;
    extern filenames_struct    filenames;
    extern int                 mpi_rank;
    extern size_t              current;

    size_t             i, j, v;
    size_t             rec;
    size_t             uidx;
    int                status;
    force_data_struct *force;
    veg_con_struct    *veg_con;
    veg_hist_struct  **veg_hist;
    soil_con_struct   *soil_con;
    double           **forcing_data;
    double          ***veg_hist_data;
    double            *dvar = NULL;
    // 
    force = grid->force;
    veg_con = &grid->veg_con[0];
    veg_hist = grid->veg_hist;
    soil_con = grid->soil_con;

    /*******************************************
      Check that required inputs were supplied
    *******************************************/

    if (!param_set.TYPE[AIR_TEMP].SUPPLIED) {
        log_err("Air temperature must be supplied as a forcing");
    }
    if (!param_set.TYPE[PREC].SUPPLIED && 
        !(param_set.TYPE[RAINF].SUPPLIED && param_set.TYPE[SNOWF].SUPPLIED)) {
        log_err("Either precipitation must be supplied as a forcing, or both rain and snow must be supplied");
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
    if (!param_set.TYPE[QAIR].SUPPLIED) {
        log_err("Specific humidity must be supplied as a forcing");
    }
    if (!param_set.TYPE[WIND].SUPPLIED) {
        log_err("Wind speed must be supplied as a forcing");
    }
    if (!param_set.TYPE[REL_HUMID].SUPPLIED) {
        log_err("Wind speed must be supplied as a forcing");
    }

    /*******************************
       Miscellaneous initialization
    *******************************/
    int has_prec = param_set.TYPE[PREC].SUPPLIED;
    int has_rain = param_set.TYPE[RAINF].SUPPLIED;
    int has_snow = param_set.TYPE[SNOWF].SUPPLIED;

    /*******************************
       read in meteorological data
    *******************************/
    // allocate memory for variables to be read
    dvar = malloc(local_domain.ncells_active * sizeof(*dvar));
    check_alloc_status(dvar, "Memory allocation error.");

    // global_param.forceoffset[0] resets every year since the met file restarts
    // every year
    // global_param.forceskip[0] should also reset to 0 after the first year
    if (current > 0 && (dmy[current].year != dmy[current - 1].year)) {
        global_param.forceoffset[0] = 0;
        global_param.forceskip[0] = 0;
        // close the forcing file for the previous year and open the forcing
        // file for the current new year
        // (forcing file for the first year should already be open in
        // get_global_param)
        if (mpi_rank == VIC_MPI_ROOT) {
            // close previous forcing file
            status = nc_close(filenames.forcing[0].nc_id);
            check_nc_status(status, "Error closing %s",
                            filenames.forcing[0].nc_filename);
            // open new forcing file
            sprintf(filenames.forcing[0].nc_filename, "%s%4d.nc",
                    filenames.f_path_pfx[0], dmy[current].year);
            status = nc_open(filenames.forcing[0].nc_filename, NC_NOWRITE,
                             &(filenames.forcing[0].nc_id));
            check_nc_status(status, "Error opening %s",
                            filenames.forcing[0].nc_filename);
        }
    }
    // only the time slice changes for the met file reads. The rest is constant
    d3start[1] = 0;
    d3start[2] = 0;
    d3count[0] = 1;

    // Air temperature: tas
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

    // Precipitation: prcp
    for (j = 0; j < NF; j++) {
        d3start[0] = global_param.forceskip[0] + global_param.forceoffset[0] +
                     j;
        get_scatter_nc_field_double(&(filenames.forcing[0]),
                                    param_set.TYPE[PREC].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            force[i].prec[j] = (double) dvar[i];
        }
    }

    // Downward solar radiation: dswrf
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

    // Downward longwave radiation: dlwrf
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

    // Wind speed: wind
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

    // vapor pressure: vp
    for (j = 0; j < NF; j++) {
        d3start[0] = global_param.forceskip[0] + global_param.forceoffset[0] +
                     j;
        get_scatter_nc_field_double(&(filenames.forcing[0]),
                                    param_set.TYPE[VP].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            force[i].vp[j] = (double) dvar[i];
        }
    }

    // Pressure: pressure
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
    // Optional inputs
    if (options.LAKES) {
        // Channel inflow to lake
        // If channel_in not supplied, it defaults to 0 from calloc
        if (param_set.TYPE[CHANNEL_IN].SUPPLIED) {
            d3start[0] = global_param.forceskip[0] +
                         global_param.forceoffset[0] + j;
            get_scatter_nc_field_double(&(filenames.forcing[0]),
                                        param_set.TYPE[CHANNEL_IN].varname,
                                        d3start, d3count, dvar);
            for (j = 0; j < NF; j++) {
                for (i = 0; i < local_domain.ncells_active; i++) {
                    force[i].channel_in[j] = (double) dvar[i];
                }
            }
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
        // Cosine of solar zenith angle
        for (j = 0; j < NF; j++) {
            for (i = 0; i < local_domain.ncells_active; i++) {
                force[i].coszen[j] = compute_coszen(
                    local_domain.locations[i].latitude,
                    local_domain.locations[i].longitude,
                    soil_con[i].time_zone_lng, dmy[current].day_in_year,
                    dmy[current].dayseconds);
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
        for (v = 0; v < options.NVEGTYPES; v++) {
            vidx = veg_con_map[i].vidx[v];
            if (vidx != NODATA_VEG) {
                for (j = 0; j < NF; j++) {
                    veg_hist[i][vidx].albedo[j] =
                        veg_con[i][vidx].albedo[dmy[current].month - 1];
                    veg_hist[i][vidx].displacement[j] =
                        veg_con[i][vidx].displacement[dmy[current].month - 1];
                    veg_hist[i][vidx].fcanopy[j] =
                        veg_con[i][vidx].fcanopy[dmy[current].month - 1];
                    veg_hist[i][vidx].LAI[j] =
                        veg_con[i][vidx].LAI[dmy[current].month - 1];
                    veg_hist[i][vidx].roughness[j] =
                        veg_con[i][vidx].roughness[dmy[current].month - 1];
                }
            }
        }
    }

    // Read veg_hist file
    if (options.LAI_SRC == FROM_VEGHIST ||
        options.FCAN_SRC == FROM_VEGHIST ||
        options.ALB_SRC == FROM_VEGHIST) {
        // global_param.forceoffset[1] resets every year since the met file restarts
        // every year
        // global_param.forceskip[1] should also reset to 0 after the first year
        if (current > 0 && (dmy[current].year != dmy[current - 1].year)) {
            global_param.forceoffset[1] = 0;
            global_param.forceskip[1] = 0;
            // close the forcing file for the previous year and open the forcing
            // file for the current new year
            // (forcing file for the first year should already be open in
            // get_global_param)
            if (mpi_rank == VIC_MPI_ROOT) {
                // close previous forcing file
                status = nc_close(filenames.forcing[1].nc_id);
                check_nc_status(status, "Error closing %s",
                                filenames.forcing[1].nc_filename);
                // open new forcing file
                sprintf(filenames.forcing[1].nc_filename, "%s%4d.nc",
                        filenames.f_path_pfx[1],
                        dmy[current].year);
                status = nc_open(filenames.forcing[1].nc_filename, NC_NOWRITE,
                                 &(filenames.forcing[1].nc_id));
                check_nc_status(status, "Error opening %s",
                                filenames.forcing[1].nc_filename);
            }
        }

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
                for (v = 0; v < options.NVEGTYPES; v++) {
                    d4start[1] = v;
                    get_scatter_nc_field_double(&(filenames.forcing[1]), "LAI",
                                                d4start, d4count, dvar);
                    for (i = 0; i < local_domain.ncells_active; i++) {
                        vidx = veg_con_map[i].vidx[v];
                        if (vidx != NODATA_VEG) {
                            veg_hist[i][vidx].LAI[j] = (double) dvar[i];
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
                for (v = 0; v < options.NVEGTYPES; v++) {
                    d4start[1] = v;
                    get_scatter_nc_field_double(&(filenames.forcing[1]),
                                                "fcanopy", d4start, d4count,
                                                dvar);
                    for (i = 0; i < local_domain.ncells_active; i++) {
                        vidx = veg_con_map[i].vidx[v];
                        if (vidx != NODATA_VEG) {
                            veg_hist[i][vidx].fcanopy[j] = (double) dvar[i];
                        }
                    }
                }
            }
        }

        // Albedo: albedo
        if (options.ALB_SRC == FROM_VEGHIST) {
            for (j = 0; j < NF; j++) {
                d4start[0] = global_param.forceskip[1] +
                             global_param.forceoffset[1] + j;
                for (v = 0; v < options.NVEGTYPES; v++) {
                    d4start[1] = v;
                    get_scatter_nc_field_double(&(filenames.forcing[1]),
                                                "albedo", d4start, d4count,
                                                dvar);
                    for (i = 0; i < local_domain.ncells_active; i++) {
                        vidx = veg_con_map[i].vidx[v];
                        if (vidx != NODATA_VEG) {
                            veg_hist[i][vidx].albedo[j] = (double) dvar[i];
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


    // allocate memory for t_offset
    t_offset = malloc(local_domain.ncells_active * sizeof(*t_offset));
    check_alloc_status(t_offset, "Memory allocation error.");

    for (i = 0; i < local_domain.ncells_active; i++) {
        if (options.SNOW_BAND > 1) {
            Tfactor = soil_con[i].Tfactor;
            t_offset[i] = Tfactor[0];
            for (band = 1; band < options.SNOW_BAND; band++) {
                if (Tfactor[band] < t_offset[i]) {
                    t_offset[i] = Tfactor[band];
                }
            }
        }
        else {
            t_offset[i] = 0;
        }
    }
    // Convert forcings into what we need and calculate missing ones
    for (i = 0; i < local_domain.ncells_active; i++) {
        for (j = 0; j < NF; j++) {
            // pressure in Pa
            force[i].pressure[j] *= PA_PER_KPA;
            // vapor pressure in Pa
            force[i].vp[j] *= PA_PER_KPA;
            // air density in kg/m3
            force[i].density[j] = air_density(force[i].air_temp[j],
                                              force[i].pressure[j]);
        }
    }

    // Checks on fcanopy and LAI
    for (i = 0; i < local_domain.ncells_active; i++) {
        for (v = 0; v < options.NVEGTYPES; v++) {
            vidx = veg_con_map[i].vidx[v];
            if (vidx != NODATA_VEG) {
                for (j = 0; j < NF; j++) {
                    if (veg_hist[i][vidx].fcanopy[j] < MIN_FCANOPY ||
                        veg_hist[i][vidx].LAI[j] == 0) {
                        if (current == 0 || options.FCAN_SRC == FROM_VEGHIST) {
                            // Only issue this warning once
                            log_warn(
                                "cell %zu, veg %d substep %zu either fcanopy "
                                "%f < minimum of %f or LAI %f == 0; setting "
                                "both LAI and fcanopy to 0", i, vidx, j,
                                veg_hist[i][vidx].fcanopy[j], MIN_FCANOPY,
                                veg_hist[i][vidx].LAI[j]);
                        }
                        veg_hist[i][vidx].fcanopy[j] = 0;
                        veg_hist[i][vidx].LAI[j] = 0;
                    }
                }
            }
        }
    }


    // Put average value in NR field
    for (i = 0; i < local_domain.ncells_active; i++) {
        force[i].air_temp[NR] = average(force[i].air_temp, NF);
        // For precipitation put total
        force[i].prec[NR] = average(force[i].prec, NF) * NF;
        force[i].shortwave[NR] = average(force[i].shortwave, NF);
        force[i].longwave[NR] = average(force[i].longwave, NF);
        force[i].pressure[NR] = average(force[i].pressure, NF);
        force[i].wind[NR] = average(force[i].wind, NF);
        force[i].vp[NR] = average(force[i].vp, NF);
        force[i].density[NR] = air_density(force[i].air_temp[NR],
                                           force[i].pressure[NR]);

        for (v = 0; v < options.NVEGTYPES; v++) {
            vidx = veg_con_map[i].vidx[v];
            if (vidx != NODATA_VEG) {
                // not the correct way to calculate average albedo in general,
                // but leave for now (it's correct if albedo is constant over
                // the model step)
                veg_hist[i][vidx].albedo[NR] = average(veg_hist[i][vidx].albedo,
                                                       NF);
                veg_hist[i][vidx].displacement[NR] = average(
                    veg_hist[i][vidx].displacement, NF);
                veg_hist[i][vidx].fcanopy[NR] = average(
                    veg_hist[i][vidx].fcanopy, NF);
                veg_hist[i][vidx].LAI[NR] = average(veg_hist[i][vidx].LAI, NF);
                veg_hist[i][vidx].roughness[NR] = average(
                    veg_hist[i][vidx].roughness, NF);
            }
        }

    log_info("Read meteorological forcing file");

    /****************************************************
       Variables in the atmos_data structure
    ****************************************************/

    for (rec = 0; rec < global_param.nrecs; rec++) {
        for (i = 0; i < NF; i++) {
            uidx = rec * NF + i;
            // temperature in Celsius
            force[rec].air_temp[i] = forcing_data[AIR_TEMP][uidx];
            // precipitation in mm/period
            if (has_prec && !has_rain && !has_snow) {
                force[rec].prec[i] = forcing_data[PREC][uidx];
                force[rec].rainf[i] = 0.;
                force[rec].snowf[i] = 0.;
            }
            else {
                if (has_rain && has_snow) {
                    force[rec].prec[i] = forcing_data[RAINF][uidx] + forcing_data[SNOWF][uidx];
                    force[rec].rainf[i] = forcing_data[RAINF][uidx];
                    force[rec].snowf[i] = forcing_data[SNOWF][uidx];
                }
                else {
                    log_err("Either prec must be supplied as a forcing, or both rain and snow must be supplied");
                }
            }
            
            // downward shortwave in W/m2
            force[rec].shortwave[i] = forcing_data[SWDOWN][uidx];
       
            // downward longwave in W/m2
            force[rec].longwave[i] = forcing_data[LWDOWN][uidx];
            // pressure in Pa
            force[rec].pressure[i] = forcing_data[PRESSURE][uidx] * PA_PER_KPA;
            // vapor pressure in Pa
            force[rec].vp[i] = q_to_vp(forcing_data[QAIR][uidx],
                                       forcing_data[PRESSURE][uidx] * PA_PER_KPA);
            // specific humidity [kg/kg]
            force[rec].Qair[i] = forcing_data[QAIR][uidx];
            // rel_humid [%]
            force[rec].rel_humid[i] = forcing_data[REL_HUMID][uidx];

            // air density in kg/m3
            force[rec].density[i] = air_density(force[rec].air_temp[i],
                                                force[rec].pressure[i],
                                                force[rec].vp[i]);
            
            // wind speed in m/s
            force[rec].wind[i] = forcing_data[WIND][uidx];

            // Cosine of solar zenith angle
            force[rec].coszen[i] = compute_coszen(soil_con->lat,
                                                  soil_con->lng,
                                                  soil_con->time_zone_lng,
                                                  dmy[rec].day_in_year,
                                                  dmy[rec].dayseconds);
            
            if (options.CARBON) {
                // Atmospheric CO2 concentration
                force[rec].Catm[i] = forcing_data[CATM][uidx];
                // Fraction of shortwave that is direct
                force[rec].fdir[i] = forcing_data[FDIR][uidx];
                // photosynthetically active radiation
                force[rec].par[i] = forcing_data[PAR][uidx];
            }
        }

        if (NF > 1) {
            force[rec].air_temp[NR] = average(force[rec].air_temp, NF);
            // For precipitation put total
            if (has_prec && !has_rain && !has_snow) {
                force[rec].prec[NR] = average(force[rec].prec, NF) * NF;
            }
            else {
                if (has_rain && has_snow) {
                    force[rec].prec[NR] = average(force[rec].prec, NF);
                    force[rec].rainf[NR] = average(force[rec].rainf, NF);
                    force[rec].snowf[NR] = average(force[rec].snowf, NF);
                }
            }
            force[rec].shortwave[NR] = average(force[rec].shortwave, NF);
            force[rec].longwave[NR] = average(force[rec].longwave, NF);
            force[rec].pressure[NR] = average(force[rec].pressure, NF);
            force[rec].vp[NR] = average(force[rec].vp, NF);
            force[rec].rel_humid[NR] = average(force[rec].rel_humid, NF);
            force[rec].density[NR] = average(force[rec].density, NF);
            force[rec].wind[NR] = average(force[rec].wind, NF);
            // for coszen, use value at noon
            force[rec].coszen[NR] = compute_coszen(soil_con->lat,
                                                   soil_con->lng,
                                                   soil_con->time_zone_lng,
                                                   dmy[rec].day_in_year,
                                                   SEC_PER_DAY / 2);
            if (options.CARBON) {
                force[rec].Catm[NR] = average(force[rec].Catm, NF);
                force[rec].fdir[NR] = average(force[rec].fdir, NF);
                force[rec].par[NR] = average(force[rec].par, NF);
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
                veg_hist[rec][v].SAI[i] =
                    veg_con[v].SAI[dmy[rec].month - 1];
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
}
