/******************************************************************************
 * @section DESCRIPTION
 *
 * Read initial model state.
 *****************************************************************************/

#include "vic_driver_shared_image.h"
#include "rout.h"

/******************************************************************************
 * @brief    Read initial model state.
 *****************************************************************************/
void
vic_restore(void)
{
    extern int                 mpi_rank;
    extern all_vars_struct    *all_vars;
    extern domain_struct       global_domain;
    extern domain_struct       local_domain;
    extern filenames_struct    filenames;
    extern metadata_struct     state_metadata[N_STATE_VARS + N_STATE_VARS_EXT];

    size_t                     i,j,k,m;
    int                       *ivar = NULL;
    int                        status;
    double                    *dvar = NULL;
    size_t                     d3count[3];
    size_t                     d3start[3];
    size_t                     d4count[4];
    size_t                     d4start[4];

    if (mpi_rank == VIC_MPI_ROOT) {
        // open initial state file
        status = nc_open(filenames.init_state.nc_filename, NC_NOWRITE,
                         &(filenames.init_state.nc_id));
        check_nc_status(status, "Error opening %s",
                        filenames.init_state.nc_filename);
    }

    // validate state file dimensions and coordinate variables
    check_init_state_file();
    // read state variables

    // allocate memory for variables to be stored
    ivar = malloc(local_domain.ncells_active * sizeof(*ivar));
    check_alloc_status(ivar, "Memory allocation error");

    dvar = malloc(local_domain.ncells_active * sizeof(*dvar));
    check_alloc_status(dvar, "Memory allocation error");

    // initialize starts and counts
    d3start[0] = 0;
    d3start[1] = 0;
    d3start[2] = 0;
    d3count[0] = 1;
    d3count[1] = global_domain.n_ny;
    d3count[2] = global_domain.n_nx;

    d4start[0] = 0;
    d4start[1] = 0;
    d4start[2] = 0;
    d4start[3] = 0;
    d4count[0] = 1;
    d4count[1] = 1;
    d4count[2] = global_domain.n_ny;
    d4count[3] = global_domain.n_nx;

    // total soil moisture
    for (m = 0; m < MAX_HRUS; m++) {
        d4start[0] = m;
        for (j = 0; j < MAX_SOILS; j++) {
            d4start[1] = j;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_SOIL_MOISTURE].varname,
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg) {
                    all_vars[i].cell[m].moist[j] = dvar[i];
                }
            }
        }
    }

    // ice content
    for (m = 0; m < MAX_HRUS; m++) {
        d4start[0] = m;
        for (j = 0; j < MAX_SOILS; j++) {
            d4start[1] = j;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_SOIL_ICE].varname,
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg) {
                    all_vars[i].cell[m].ice[j] = dvar[i];
                }
            }
        }
    }

    // dew storage: tmpval = veg_var[veg].Wdew;
    for (m = 0; m < MAX_HRUS; m++) {
        d3start[0] = m;
        get_scatter_nc_field_double(&(filenames.init_state),
                                    state_metadata[STATE_CANOPY_WATER].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                all_vars[i].veg_var[m].Wdew = dvar[i];
            }
        }
    }

    // Tcanopy: tmpval = cell[veg].Tcanopy;
    for (m = 0; m < MAX_HRUS; m++) {
        d3start[0] = m;
        get_scatter_nc_field_double(&(filenames.init_state),
                                    state_metadata[STATE_TCANOPY].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                all_vars[i].energy[m].Tcanopy = dvar[i];
            }
        }
    }

    // Tfoliage: tmpval = cell[veg].Tfoliage;
    for (m = 0; m < MAX_HRUS; m++) {
        d3start[0] = m;
        get_scatter_nc_field_double(&(filenames.init_state),
                                    state_metadata[STATE_TCANOPY].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                all_vars[i].energy[m].Tfoliage = dvar[i];
            }
        }
    }

    // Tgrnd: tmpval = cell[veg].Tgrnd;
    for (m = 0; m < MAX_HRUS; m++) {
        d3start[0] = m;
        get_scatter_nc_field_double(&(filenames.init_state),
                                    state_metadata[STATE_TCANOPY].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                all_vars[i].energy[m].Tgrnd = dvar[i];
            }
        }
    }

    // snow age: snow[veg].SnowAge
    for (m = 0; m < MAX_HRUS; m++) {
        d3start[0] = m;
        get_scatter_nc_field_int(&(filenames.init_state),
                                 state_metadata[STATE_SNOW_AGE].varname,
                                 d3start, d3count, ivar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                all_vars[i].snow[m].snowage = ivar[i];
            }
        }
    }

    // snow covered fraction: snow[veg].coverage
    for (m = 0; m < MAX_HRUS; m++) {
        d3start[0] = m;
        get_scatter_nc_field_int(&(filenames.init_state),
                                 state_metadata[STATE_SNOW_COVERAGE].varname,
                                 d3start, d3count, ivar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                all_vars[i].snow[m].coverage = ivar[i];
            }
        }
    }

    // pack_ice: snow[veg].pack_ice
    for (m = 0; m < MAX_HRUS; m++) {
        d4start[0] = m;
        for (k = 0; k < MAX_SNOWS; k++) {
            d4start[1] = k;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_SNOW_PACK_ICE].varname,
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg && 
                                all_vars[i].snow[m].Nsnow > 0) {
                    all_vars[i].snow[m].pack_ice[k] = dvar[i];
                }
            }
        }
    }

    // snow water equivalent: snow[veg].swq
    for (m = 0; m < MAX_HRUS; m++) {
        d3start[0] = m;
        get_scatter_nc_field_double(&(filenames.init_state),
                                    state_metadata[
                                        STATE_SNOW_WATER_EQUIVALENT].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                all_vars[i].snow[m].swq = dvar[i];
            }
        }
    }

    // snow theta_ice: snow[veg].theta_ice
    for (m = 0; m < MAX_HRUS; m++) {
        d4start[0] = m;
        for (k = 0; k < MAX_SNOWS; k++) {
            d4start[1] = k;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_SNOW_THETA_ICE].varname,
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg &&
                                all_vars[i].snow[m].Nsnow > 0) {
                    all_vars[i].snow[m].theta_ice[k] = dvar[i];
                }
            }
        }
    }

    // snow surface water: snow[veg].theta_liq
    for (m = 0; m < MAX_HRUS; m++) {
        d4start[0] = m;
        for (k = 0; k < MAX_SNOWS; k++) {
            d4start[1] = k;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_SNOW_THETA_LIQ].varname,
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg &&
                                all_vars[i].snow[m].Nsnow > 0) {
                    all_vars[i].snow[m].theta_liq[k] = dvar[i];
                }
            }
        }
    }

    // snow pack temperature: snow[veg].pack_temp
    for (m = 0; m < MAX_HRUS; m++) {
        d4start[0] = m;
        for (k = 0; k < MAX_SNOWS; k++) {
            d4start[1] = k;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_SNOW_PACK_TEMP].varname,
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg &&
                                all_vars[i].snow[m].Nsnow > 0) {
                    all_vars[i].snow[m].pack_T[k] = dvar[i];
                }
            }
        }
    }

    // snow pack water: snow[veg].pack_liq
    for (m = 0; m < MAX_HRUS; m++) {
        d4start[0] = m;
        for (k = 0; k < MAX_SNOWS; k++) {
            d4start[1] = k;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_SNOW_PACK_LIQ].varname,
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg &&
                                all_vars[i].snow[m].Nsnow > 0) {
                    all_vars[i].snow[m].pack_liq[k] = dvar[i];
                }
            }
        }
    }

    // snow pack porosity: snow[veg].porosity
    for (m = 0; m < MAX_HRUS; m++) {
        d4start[0] = m;
        for (k = 0; k < MAX_SNOWS; k++) {
            d4start[1] = k;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_SNOW_POROSITY].varname,
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg &&
                                all_vars[i].snow[m].Nsnow > 0) {
                    all_vars[i].snow[m].porosity[k] = dvar[i];
                }
            }
        }
    }

    // snow density: snow[veg].density
    for (m = 0; m < MAX_HRUS; m++) {
        d4start[0] = m;
        for (k = 0; k < MAX_SNOWS; k++) {
            d4start[1] = k;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_SNOW_DENSITY].varname,
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg &&
                                all_vars[i].snow[m].Nsnow > 0) {
                    all_vars[i].snow[m].density[k] = dvar[i];
                }
            }
        }
    }

    // last_swq: snow[veg].last_swq
    for (m = 0; m < MAX_HRUS; m++) {
        d3start[0] = m;
        get_scatter_nc_field_double(&(filenames.init_state),
                                    state_metadata[STATE_SNOW_OLDSWQ].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                all_vars[i].snow[m].last_swq = dvar[i];
            }
        }
    }

    // canopy swe storage: veg_var[veg].canopy_swe
    for (m = 0; m < MAX_HRUS; m++) {
        d3start[0] = m;
        get_scatter_nc_field_double(&(filenames.init_state),
                                    state_metadata[STATE_SNOW_CANOPY].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                all_vars[i].veg_var[m].canopy_swq = dvar[i];
            }
        }
    }

    // thermal node temperatures: energy[veg].T[nidx]
    for (m = 0; m < MAX_HRUS; m++) {
        d4start[0] = m;
        for (j = 0; j < MAX_NODES; j++) {
            d4start[1] = j;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_NODE_TEMP].varname,
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg &&
                                k < all_vars[i].cell[m].Nnode) {
                    all_vars[i].energy[m].T[j] = dvar[i];
                }
            }
        }
    }

    // Foliage temperature: energy[veg].Tfoliage
    for (m = 0; m < MAX_HRUS; m++) {
        d3start[0] = m;
        get_scatter_nc_field_double(&(filenames.init_state),
                                    state_metadata[STATE_FOLIAGE_TEMPERATURE].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                all_vars[i].energy[m].Tfoliage = dvar[i];
            }
        }
    }

    // routing ring
    vic_restore_rout_extension(&(filenames.init_state), state_metadata);

    free(ivar);
    free(dvar);

    // close initial state file
    if (mpi_rank == VIC_MPI_ROOT) {
        status = nc_close(filenames.init_state.nc_id);
        check_nc_status(status, "Error closing %s",
                        filenames.init_state.nc_filename);
    }
}

/******************************************************************************
 * @brief    Check that the initial state file matches the global parameter
             settings
 *****************************************************************************/
void
check_init_state_file(void)
{
    extern filenames_struct filenames;
    extern domain_struct    global_domain;
    extern domain_struct    local_domain;
    extern option_struct    options;
    extern soil_con_struct *soil_con;
    extern int              mpi_rank;

    int                     status;
    size_t                  dimlen;
    size_t                  i;
    size_t                  j;
    size_t                  d1count[1];
    size_t                  d1start[1];
    size_t                  d2count[2];
    size_t                  d2start[2];
    size_t                  d3count[3];
    size_t                  d3start[3];
    int                     lon_var_id;
    int                     lat_var_id;
    double                 *dvar;
    double                  rtol = 0.0; // maybe move this to a .h file
    double                  abs_tol = 0.0001; // maybe move this to a .h file

    // read and validate dimension lengths
    if (mpi_rank == VIC_MPI_ROOT) {
        dimlen = get_nc_dimension(&(filenames.init_state),
                                  global_domain.info.x_dim);
        if (dimlen != global_domain.n_nx) {
            log_err("Number of grid columns in state file does not "
                    "match parameter file");
        }
        dimlen = get_nc_dimension(&(filenames.init_state),
                                  global_domain.info.y_dim);
        if (dimlen != global_domain.n_ny) {
            log_err("Number of grid rows in state file does not "
                    "match parameter file");
        }
        dimlen = get_nc_dimension(&(filenames.init_state), "veg_class");
        if (dimlen != options.NVEGTYPES) {
            log_err("Number of veg classes in state file does not "
                    "match parameter file");
        }
        dimlen = get_nc_dimension(&(filenames.init_state), "snow_band");
        if (dimlen != options.SNOW_BAND) {
            log_err("Number of snow bands in state file does not "
                    "match parameter file");
        }
        dimlen = get_nc_dimension(&(filenames.init_state), "nlayer");
        if (dimlen != options.Nlayer) {
            log_err("Number of soil layers in state file does not "
                    "match parameter file");
        }
        dimlen = get_nc_dimension(&(filenames.init_state), "frost_area");
        if (dimlen != options.Nfrost) {
            log_err("Number of frost areas in state file does not "
                    "match parameter file");
        }
        dimlen = get_nc_dimension(&(filenames.init_state), "soil_node");
        if (dimlen != MAX_SOILS) {
            log_err("Number of soil nodes in state file does not "
                    "match parameter file");
        }
    }

    // read dimension variables

    // lat/lon
    if (mpi_rank == VIC_MPI_ROOT) {
        status = nc_inq_varid(filenames.init_state.nc_id,
                              global_domain.info.lon_var, &lon_var_id);
        check_nc_status(status, "Unable to find variable \"%s\" in %s",
                        global_domain.info.lon_var,
                        filenames.init_state.nc_filename);
        status = nc_inq_varid(filenames.init_state.nc_id,
                              global_domain.info.lat_var, &lat_var_id);
        check_nc_status(status, "Unable to find variable \"%s\" in %s",
                        global_domain.info.lat_var,
                        filenames.init_state.nc_filename);
        if (global_domain.info.n_coord_dims == 1) {
            d1start[0] = 0;
            dvar = calloc(global_domain.n_nx, sizeof(*dvar));
            check_alloc_status(dvar, "Memory allocation error");

            d1count[0] = global_domain.n_nx;
            status = nc_get_vara_double(filenames.init_state.nc_id, lon_var_id,
                                        d1start, d1count, dvar);
            check_nc_status(status, "Error reading data from \"%s\" in %s",
                            global_domain.info.lon_var,
                            filenames.init_state.nc_filename);
            // implicitly nested loop over ni and nj with j set to 0
            for (i = 0; i < global_domain.n_nx; i++) {
                if (!assert_close_double(dvar[i],
                                         global_domain.locations[i].longitude,
                                         rtol,
                                         abs_tol)) {
                    log_err("Longitudes in initial state file do not "
                            "match parameter file");
                }
            }
            free(dvar);

            dvar = calloc(global_domain.n_ny, sizeof(*dvar));
            check_alloc_status(dvar, "Memory allocation error");

            d1count[0] = global_domain.n_ny;
            status = nc_get_vara_double(filenames.init_state.nc_id, lat_var_id,
                                        d1start, d1count, dvar);
            check_nc_status(status, "Error reading data from \"%s\" in %s",
                            global_domain.info.lat_var,
                            filenames.init_state.nc_filename);
            // implicitly nested loop over ni and nj with i set to 0;
            // j stride = n_nx
            for (j = 0; j < global_domain.n_ny; j++) {
                if (!assert_close_double(dvar[j],
                                         global_domain.locations[j *
                                                                 global_domain.
                                                                 n_nx]
                                         .latitude, rtol,
                                         abs_tol)) {
                    log_err("Latitudes in initial state file do not "
                            "match parameter file");
                }
            }
            free(dvar);
        }
        else if (global_domain.info.n_coord_dims == 2) {
            d2start[0] = 0;
            d2start[1] = 0;
            dvar =
                calloc(global_domain.n_ny * global_domain.n_nx, sizeof(*dvar));
            check_alloc_status(dvar, "Memory allocation error");

            d2count[0] = global_domain.n_ny;
            d2count[1] = global_domain.n_nx;
            status = nc_get_vara_double(filenames.init_state.nc_id, lon_var_id,
                                        d2start, d2count, dvar);
            check_nc_status(status, "Error reading data from \"%s\" in %s",
                            global_domain.info.lon_var,
                            filenames.init_state.nc_filename);
            for (i = 0; i < global_domain.n_ny * global_domain.n_nx; i++) {
                if (dvar[i] != (double) global_domain.locations[i].longitude) {
                    log_err("Longitudes in initial state file do not "
                            "match parameter file");
                }
            }
            status = nc_get_vara_double(filenames.init_state.nc_id, lat_var_id,
                                        d2start, d2count, dvar);
            check_nc_status(status, "Error reading data from \"%s\" in %s",
                            global_domain.info.lat_var,
                            filenames.init_state.nc_filename);
            for (i = 0; i < global_domain.n_ny * global_domain.n_nx; i++) {
                if (dvar[i] != (double) global_domain.locations[i].latitude) {
                    log_err("Latitudes in initial state file do not "
                            "match parameter file");
                }
            }
            free(dvar);
        }
        else {
            log_err("global_domain.info.n_coord_dims should be 1 or 2");
        }
    }

    // initialize dvar for soil thermal node deltas and depths
    dvar = malloc(local_domain.ncells_active * sizeof(*dvar));
    check_alloc_status(dvar, "Memory allocation error");

    // soil thermal node deltas (dimension: node, lat, lon)
    d3start[0] = 0;
    d3start[1] = 0;
    d3start[2] = 0;
    d3count[0] = 1;
    d3count[1] = global_domain.n_ny;
    d3count[2] = global_domain.n_nx;
    for (j = 0; j < MAX_SOILS; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(&(filenames.init_state),
                                    "dz_soil",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (dvar[i] != soil_con[i].dz_soil[j]) {
                log_err("Soil node intervals in state file do not match "
                        "those computed by VIC");
            }
        }
    }

    // soil thermal node depths
    d3start[0] = 0;
    d3start[1] = 0;
    d3start[2] = 0;
    d3count[0] = 1;
    d3count[1] = global_domain.n_ny;
    d3count[2] = global_domain.n_nx;
    for (j = 0; j < MAX_SOILS; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(&(filenames.init_state),
                                    "node_depth",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (dvar[i] != soil_con[i].Zsum_soil[j]) {
                log_err("Soil node depths in state file do not match "
                        "those computed by VIC");
            }
        }
    }
    free(dvar);
}
