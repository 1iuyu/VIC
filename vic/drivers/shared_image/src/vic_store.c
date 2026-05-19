/******************************************************************************
 * @section DESCRIPTION
 *
 * Save model state.
 *****************************************************************************/

#include "vic_driver_shared_image.h"
#include "rout.h"

/******************************************************************************
 * @brief    Save model state.
 *****************************************************************************/
void
vic_store(dmy_struct *dmy_state,
          char       *filename)
{
    extern filenames_struct    filenames;
    extern all_vars_struct    *all_vars;
    extern domain_struct       local_domain;
    extern option_struct       options;
    extern int                 mpi_rank;

    int                        status;
    size_t                     i,j,k,m;
    int                       *ivar = NULL;
    double                    *dvar = NULL;
    size_t                     d3start[3];
    size_t                     d4start[4];
    nc_file_struct             nc_state_file;
    nc_var_struct             *nc_var;

    set_nc_state_file_info(&nc_state_file);

    // create netcdf file for storing model state
    if (options.STATENAME_CESM) {
        sprintf(filename, "%s.vic.r.%04i-%02i-%02i-%05u.nc",
                filenames.statefile, dmy_state->year,
                dmy_state->month, dmy_state->day,
                dmy_state->dayseconds);
    }
    else {
        sprintf(filename, "%s.%04i%02i%02i_%05u.nc",
                filenames.statefile, dmy_state->year,
                dmy_state->month, dmy_state->day,
                dmy_state->dayseconds);
    }

    initialize_state_file(filename, &nc_state_file, dmy_state);

    if (mpi_rank == VIC_MPI_ROOT) {
        debug("writing state file: %s", filename);
    }

    // write state variables

    // allocate memory for variables to be stored
    ivar = malloc(local_domain.ncells_active * sizeof(*ivar));
    check_alloc_status(ivar, "Memory allocation error");


    dvar = malloc(local_domain.ncells_active * sizeof(*dvar));
    check_alloc_status(dvar, "Memory allocation error");


    // initialize starts and counts
    d3start[0] = 0;
    d3start[1] = 0;
    d3start[2] = 0;

    d4start[0] = 0;
    d4start[1] = 0;
    d4start[2] = 0;
    d4start[3] = 0;

    // set missing values
    for (i = 0; i < local_domain.ncells_active; i++) {
        ivar[i] = nc_state_file.i_fillvalue;
        dvar[i] = nc_state_file.d_fillvalue;
    }

    // total soil moisture
    nc_var = &(nc_state_file.nc_vars[STATE_SOIL_MOISTURE]);
    for (m = 0; m < options.MAX_HRU; m++) {
        d4start[0] = m;
        for (j = 0; j < MAX_SOILS; j++) {
            d4start[1] = j;
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg) {
                    dvar[i] =
                        (double) all_vars[i].cell[m].moist[j];
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.nc_id,
                                       nc_var->nc_varid,
                                       nc_state_file.d_fillvalue,
                                       d4start, nc_var->nc_counts, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }

    // ice content
    nc_var = &(nc_state_file.nc_vars[STATE_SOIL_ICE]);
    for (m = 0; m < options.MAX_HRU; m++) {
        d4start[0] = m;
        for (j = 0; j < MAX_SOILS; j++) {
            d4start[1] = j;
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg) {
                    dvar[i] =
                        (double) all_vars[i].cell[m].ice[j];
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.nc_id,
                                       nc_var->nc_varid,
                                       nc_state_file.d_fillvalue,
                                       d4start, nc_var->nc_counts,
                                       dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }


    // dew storage: tmpval = veg_var[veg].Wdew;
    nc_var = &(nc_state_file.nc_vars[STATE_CANOPY_WATER]);
    for (m = 0; m < options.MAX_HRU; m++) {
        d3start[0] = m;
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                dvar[i] = (double) all_vars[i].veg_var[m].Wdew;
            }
            else {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
        gather_put_nc_field_double(nc_state_file.nc_id,
                                   nc_var->nc_varid,
                                   nc_state_file.d_fillvalue,
                                   d3start, nc_var->nc_counts, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.d_fillvalue;
        }
    }

    // snow age: snow[veg].SnowAge
    nc_var = &(nc_state_file.nc_vars[STATE_SNOW_AGE]);
    for (m = 0; m < options.MAX_HRU; m++) {
        d3start[0] = m;
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                ivar[i] = (int) all_vars[i].snow[m].snowage;
            }
            else {
                ivar[i] = nc_state_file.i_fillvalue;
            }
        }
        gather_put_nc_field_int(nc_state_file.nc_id,
                                nc_var->nc_varid,
                                nc_state_file.i_fillvalue,
                                d3start, nc_var->nc_counts, ivar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            ivar[i] = nc_state_file.i_fillvalue;
        }
    }

    // snow covered fraction: snow[veg].coverage
    nc_var = &(nc_state_file.nc_vars[STATE_SNOW_COVERAGE]);
    for (m = 0; m < options.MAX_HRU; m++) {
        d3start[0] = m;
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                dvar[i] = (double) all_vars[i].snow[m].coverage;
            }
            else {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
        gather_put_nc_field_double(nc_state_file.nc_id,
                                   nc_var->nc_varid,
                                   nc_state_file.d_fillvalue,
                                   d3start, nc_var->nc_counts, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.d_fillvalue;
        }
    }


    // snow water equivalent: snow[veg].swq
    nc_var = &(nc_state_file.nc_vars[STATE_SNOW_WATER_EQUIVALENT]);
    for (m = 0; m < options.MAX_HRU; m++) {
        d3start[0] = m;
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                dvar[i] = (double) all_vars[i].snow[m].swq;
            }
            else {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
        gather_put_nc_field_double(nc_state_file.nc_id,
                                   nc_var->nc_varid,
                                   nc_state_file.d_fillvalue,
                                   d3start, nc_var->nc_counts, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.d_fillvalue;
        }
    }


    // snow theta_ice: snow[veg].theta_ice
    nc_var = &(nc_state_file.nc_vars[STATE_SNOW_THETA_ICE]);
    for (m = 0; m < options.MAX_HRU; m++) {
        d4start[0] = m;
        for (k = 0; k < MAX_SNOWS; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg && 
                                k < all_vars[i].snow[m].Nsnow) {
                    dvar[i] = (double) all_vars[i].snow[m].theta_ice[k];
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.nc_id,
                                       nc_var->nc_varid,
                                       nc_state_file.d_fillvalue,
                                       d4start, nc_var->nc_counts, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }


    // snow surface water: snow[veg].theta_liq
    nc_var = &(nc_state_file.nc_vars[STATE_SNOW_THETA_LIQ]);
    for (m = 0; m < options.MAX_HRU; m++) {
        d4start[0] = m;
        for (k = 0; k < MAX_SNOWS; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg &&
                            k < all_vars[i].snow[m].Nsnow) {
                    dvar[i] = (double) all_vars[i].snow[m].theta_liq[k];
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.nc_id,
                                       nc_var->nc_varid,
                                       nc_state_file.d_fillvalue,
                                       d4start, nc_var->nc_counts, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }


    // snow pack temperature: snow[veg].pack_temp
    nc_var = &(nc_state_file.nc_vars[STATE_SNOW_PACK_TEMP]);
    for (m = 0; m < options.MAX_HRU; m++) {
        d4start[0] = m;
        for (k = 0; k < MAX_SNOWS; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg &&
                            k < all_vars[i].snow[m].Nsnow) {
                    dvar[i] = (double) all_vars[i].snow[m].pack_T[k];
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.nc_id,
                                       nc_var->nc_varid,
                                       nc_state_file.d_fillvalue,
                                       d4start, nc_var->nc_counts, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }


    // snow pack water: snow[veg].pack_liq
    nc_var = &(nc_state_file.nc_vars[STATE_SNOW_PACK_LIQ]);
    for (m = 0; m < options.MAX_HRU; m++) {
        d4start[0] = m;
        for (k = 0; k < MAX_SNOWS; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg &&
                                k < all_vars[i].snow[m].Nsnow) {
                    dvar[i] = (double) all_vars[i].snow[m].pack_liq[k];
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.nc_id,
                                       nc_var->nc_varid,
                                       nc_state_file.d_fillvalue,
                                       d4start, nc_var->nc_counts, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }

    // snow pack porosity: snow[veg].porosity
    nc_var = &(nc_state_file.nc_vars[STATE_SNOW_POROSITY]);
    for (m = 0; m < options.MAX_HRU; m++) {
        d4start[0] = m;
        for (k = 0; k < MAX_SNOWS; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg &&
                                k < all_vars[i].snow[m].Nsnow) {
                    dvar[i] = (double) all_vars[i].snow[m].porosity[k];
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.nc_id,
                                       nc_var->nc_varid,
                                       nc_state_file.d_fillvalue,
                                       d4start, nc_var->nc_counts, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }

    // snow density: snow[veg].density
    nc_var = &(nc_state_file.nc_vars[STATE_SNOW_DENSITY]);
    for (m = 0; m < options.MAX_HRU; m++) {
        d4start[0] = m;
        for (k = 0; k < MAX_SNOWS; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg &&
                                    k < all_vars[i].snow[m].Nsnow) {
                    dvar[i] = (double) all_vars[i].snow[m].density[k];
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.nc_id,
                                    nc_var->nc_varid,
                                    nc_state_file.d_fillvalue,
                                    d3start, nc_var->nc_counts, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }


    // last_swq: snow[veg].last_swq
    nc_var = &(nc_state_file.nc_vars[STATE_SNOW_OLDSWQ]);
    for (m = 0; m < options.MAX_HRU; m++) {
        d3start[0] = m;
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                dvar[i] = (double) all_vars[i].snow[m].last_swq;
            }
            else {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
        gather_put_nc_field_double(nc_state_file.nc_id,
                                   nc_var->nc_varid,
                                   nc_state_file.d_fillvalue,
                                   d3start, nc_var->nc_counts, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.d_fillvalue;
        }
    }


    // canopy swe storage: veg_var[veg].canopy_swe
    nc_var = &(nc_state_file.nc_vars[STATE_SNOW_CANOPY]);
    for (m = 0; m < options.MAX_HRU; m++) {
        d3start[0] = m;
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                dvar[i] = (double) all_vars[i].veg_var[m].canopy_swq;
            }
            else {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
        gather_put_nc_field_double(nc_state_file.nc_id,
                                   nc_var->nc_varid,
                                   nc_state_file.d_fillvalue,
                                   d3start, nc_var->nc_counts, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.d_fillvalue;
        }
    }


    // soil node temperatures: energy[veg].T[nidx]
    nc_var = &(nc_state_file.nc_vars[STATE_NODE_TEMP]);
    for (m = 0; m < options.MAX_HRU; m++) {
        d4start[0] = m;
        for (k = 0; k < MAX_NODES; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg &&
                                k < all_vars[i].cell[m].Nnode) {
                    dvar[i] = (double) all_vars[i].energy[m].T[k];
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.nc_id,
                                       nc_var->nc_varid,
                                       nc_state_file.d_fillvalue,
                                       d4start, nc_var->nc_counts, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }


    // Foliage temperature: energy[veg].Tfoliage
    nc_var = &(nc_state_file.nc_vars[STATE_FOLIAGE_TEMPERATURE]);
    for (m = 0; m < options.MAX_HRU; m++) {
        d3start[0] = m;
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                dvar[i] = (double) all_vars[i].energy[m].Tfoliage;
            }
            else {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
        gather_put_nc_field_double(nc_state_file.nc_id,
                                   nc_var->nc_varid,
                                   nc_state_file.d_fillvalue,
                                   d3start, nc_var->nc_counts, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.d_fillvalue;
        }
    }

    // store extension variables
    vic_store_rout_extension(&nc_state_file);


    // close the netcdf file if it is still open
    if (mpi_rank == VIC_MPI_ROOT) {
        if (nc_state_file.open == true) {
            status = nc_close(nc_state_file.nc_id);
            check_nc_status(status, "Error closing %s", filename);
        }
    }

    free(ivar);
    free(dvar);
    free(nc_state_file.nc_vars);
}

/******************************************************************************
 * @brief   Setup state file netcdf structure
 *****************************************************************************/
void
set_nc_state_file_info(nc_file_struct *nc_state_file)
{
    extern option_struct options;
    extern domain_struct global_domain;

    // set fill values
    nc_state_file->c_fillvalue = NC_FILL_CHAR;
    nc_state_file->s_fillvalue = NC_FILL_SHORT;
    nc_state_file->i_fillvalue = NC_FILL_INT;
    nc_state_file->d_fillvalue = NC_FILL_DOUBLE;
    nc_state_file->f_fillvalue = NC_FILL_FLOAT;

    // set ids to MISSING
    nc_state_file->nc_id = MISSING;
    nc_state_file->band_dimid = MISSING;
    nc_state_file->front_dimid = MISSING;
    nc_state_file->layer_dimid = MISSING;
    nc_state_file->ni_dimid = MISSING;
    nc_state_file->nj_dimid = MISSING;
    nc_state_file->node_dimid = MISSING;
    nc_state_file->time_dimid = MISSING;
    nc_state_file->veg_dimid = MISSING;

    // set dimension sizes
    nc_state_file->band_size = options.SNOW_BAND;
    nc_state_file->front_size = MAX_FRONTS;
    nc_state_file->layer_size = options.Nlayer;
    nc_state_file->ni_size = global_domain.n_nx;
    nc_state_file->nj_size = global_domain.n_ny;
    nc_state_file->node_size = MAX_SOILS;
    nc_state_file->time_size = NC_UNLIMITED;
    nc_state_file->veg_size = options.MAX_HRU;

    // set ids and dimension sizes of the extension variables
    set_nc_state_file_info_rout_extension(nc_state_file);

    // allocate memory for nc_vars
    nc_state_file->nc_vars =
        calloc(N_STATE_VARS + N_STATE_VARS_EXT,
               sizeof(*(nc_state_file->nc_vars)));
    check_alloc_status(nc_state_file->nc_vars, "Memory allocation error");
}

/******************************************************************************
 * @brief   Setup state variable dimensions, types, etc.
 *****************************************************************************/
void
set_nc_state_var_info(nc_file_struct *nc)
{
    size_t i, j;

    for (i = 0; i < N_STATE_VARS; i++) {
        nc->nc_vars[i].nc_varid = i;
        for (j = 0; j < MAXDIMS; j++) {
            nc->nc_vars[i].nc_dimids[j] = -1;
            nc->nc_vars[i].nc_counts[j] = 0;
        }
        nc->nc_vars[i].nc_dims = 0;

        switch (i) {
        case STATE_SNOW_NSNOW:
            nc->nc_vars[i].nc_type = NC_INT;
            break;
        default:
            nc->nc_vars[i].nc_type = NC_DOUBLE;
        }

        // Set the number of dimensions and dimids for each state variable
        switch (i) {
        case STATE_SOIL_MOISTURE:
        case STATE_SOIL_ICE:
        case STATE_SNOW_PACK_ICE:
        case STATE_SNOW_THETA_ICE:
        case STATE_SNOW_THETA_LIQ:
        case STATE_SNOW_PACK_TEMP:
        case STATE_SNOW_PACK_LIQ:
        case STATE_SNOW_POROSITY:
            // 4d vars [veg, layer, j, i]
            nc->nc_vars[i].nc_dims = 4;
            nc->nc_vars[i].nc_dimids[0] = nc->veg_dimid;
            nc->nc_vars[i].nc_dimids[1] = nc->layer_dimid;
            nc->nc_vars[i].nc_dimids[2] = nc->nj_dimid;
            nc->nc_vars[i].nc_dimids[3] = nc->ni_dimid;
            nc->nc_vars[i].nc_counts[0] = 1;
            nc->nc_vars[i].nc_counts[1] = 1;
            nc->nc_vars[i].nc_counts[2] = nc->nj_size;
            nc->nc_vars[i].nc_counts[3] = nc->ni_size;
            break;
        case STATE_CANOPY_WATER:
        case STATE_TCANOPY:
        case STATE_ANNUALNPP:
        case STATE_ANNUALNPPPREV:
        case STATE_CLITTER:
        case STATE_CINTER:
        case STATE_CSLOW:
        case STATE_SNOW_AGE:
        case STATE_SNOW_COVERAGE:
        case STATE_SNOW_WATER_EQUIVALENT:
        case STATE_SNOW_DENSITY:
        case STATE_SNOW_OLDSWQ:
        case STATE_SNOW_CANOPY:
        case STATE_FOLIAGE_TEMPERATURE:
            // 3d vars [veg, j, i]
            nc->nc_vars[i].nc_dims = 3;
            nc->nc_vars[i].nc_dimids[0] = nc->veg_dimid;
            nc->nc_vars[i].nc_dimids[1] = nc->nj_dimid;
            nc->nc_vars[i].nc_dimids[2] = nc->ni_dimid;
            nc->nc_vars[i].nc_counts[0] = 1;
            nc->nc_vars[i].nc_counts[1] = nc->nj_size;
            nc->nc_vars[i].nc_counts[2] = nc->ni_size;
            break;
        case STATE_NODE_TEMP:
            // 4d vars [veg, node, j, i]
            nc->nc_vars[i].nc_dims = 4;
            nc->nc_vars[i].nc_dimids[0] = nc->veg_dimid;
            nc->nc_vars[i].nc_dimids[1] = nc->node_dimid;
            nc->nc_vars[i].nc_dimids[2] = nc->nj_dimid;
            nc->nc_vars[i].nc_dimids[3] = nc->ni_dimid;
            nc->nc_vars[i].nc_counts[0] = 1;
            nc->nc_vars[i].nc_counts[1] = 1;
            nc->nc_vars[i].nc_counts[2] = nc->nj_size;
            nc->nc_vars[i].nc_counts[3] = nc->ni_size;
            break;
        default:
            log_err("state variable %zu not found when setting dimensions", i);
        }

        if (nc->nc_vars[i].nc_dims > MAXDIMS) {
            log_err("Too many dimensions specified in variable %zu", i);
        }
    }
    set_nc_state_var_info_rout_extension(nc);
}

/******************************************************************************
 * @brief   Initialize state file by creating dimensions, variables,
            and adding metadata.
 *****************************************************************************/
void
initialize_state_file(char           *filename,
                      nc_file_struct *nc_state_file,
                      dmy_struct     *dmy_state)
{
    extern option_struct       options;
    extern domain_struct       global_domain;
    extern domain_struct       local_domain;
    extern global_param_struct global_param;
    extern metadata_struct     state_metadata[N_STATE_VARS + N_STATE_VARS_EXT];
    extern soil_con_struct    *soil_con;
    extern int                 mpi_rank;

    int                        status;
    int                        dimids[MAXDIMS];
    int                        old_fill_mode;
    size_t                     i;
    size_t                     j;
    size_t                     dcount[MAXDIMS];
    size_t                     dstart[MAXDIMS];
    int                        lon_var_id;
    int                        lat_var_id;
    int                        veg_var_id;
    int                        snow_band_var_id;
    int                        layer_var_id;
    int                        dz_node_var_id;
    int                        node_depth_var_id;
    char                       unit_str[MAXSTRING];
    char                       str[MAXSTRING];
    size_t                     ndims;
    double                     dtime;
    double                    *dvar = NULL;
    int                       *ivar = NULL;

    // open the netcdf file
    if (mpi_rank == VIC_MPI_ROOT) {
        status = nc_create(filename, get_nc_mode(options.STATE_FORMAT),
                           &(nc_state_file->nc_id));
        check_nc_status(status, "Error creating %s", filename);
        nc_state_file->open = true;
    }

    if (mpi_rank == VIC_MPI_ROOT) {
        // Set netcdf file global attributes
        set_global_nc_attributes(nc_state_file->nc_id, NC_STATE_FILE);

        // set the NC_FILL attribute
        status = nc_set_fill(nc_state_file->nc_id, NC_FILL, &old_fill_mode);
        check_nc_status(status, "Error setting fill value in %s", filename);

        // define the time dimension
        status = nc_def_dim(nc_state_file->nc_id, "time",
                            nc_state_file->time_size,
                            &(nc_state_file->time_dimid));
        check_nc_status(status, "Error defining time dimenension in %s",
                        filename);

        // define the variable time
        status = nc_def_var(nc_state_file->nc_id, "time", NC_DOUBLE, 1,
                            &(nc_state_file->time_dimid),
                            &(nc_state_file->time_varid));
        check_nc_status(status, "Error defining time variable in %s", filename);
        status = nc_put_att_text(nc_state_file->nc_id,
                                 nc_state_file->time_varid,
                                 "standard_name", strlen("time"), "time");
        check_nc_status(status, "Error adding attribute in %s", filename);

        // adding units attribute to time variable
        str_from_time_units(global_param.time_units, unit_str);

        snprintf(str, sizeof(str), "%s since %s", unit_str, global_param.time_origin_str);

        status = nc_put_att_text(nc_state_file->nc_id,
                                 nc_state_file->time_varid,
                                 "units", strlen(str), str);
        check_nc_status(status, "Error adding attribute in %s", filename);

        // adding calendar attribute to time variable
        str_from_calendar(global_param.calendar, str);

        status = nc_put_att_text(nc_state_file->nc_id,
                                 nc_state_file->time_varid,
                                 "calendar", strlen(str), str);
        check_nc_status(status, "Error adding calendar attribute in %s",
                        filename);

        // define netcdf dimensions
        status = nc_def_dim(nc_state_file->nc_id, global_domain.info.x_dim,
                            nc_state_file->ni_size, &(nc_state_file->ni_dimid));
        check_nc_status(status, "Error defining \"%s\" in %s",
                        global_domain.info.x_dim,
                        filename);

        status = nc_def_dim(nc_state_file->nc_id, global_domain.info.y_dim,
                            nc_state_file->nj_size, &(nc_state_file->nj_dimid));
        check_nc_status(status, "Error defining \"%s\" in %s",
                        global_domain.info.y_dim,
                        filename);

        status = nc_def_dim(nc_state_file->nc_id, "veg_class",
                            nc_state_file->veg_size,
                            &(nc_state_file->veg_dimid));
        check_nc_status(status, "Error defining veg_class in %s", filename);

        status = nc_def_dim(nc_state_file->nc_id, "snow_band",
                            nc_state_file->band_size,
                            &(nc_state_file->band_dimid));
        check_nc_status(status, "Error defining snow_band in %s", filename);

        status = nc_def_dim(nc_state_file->nc_id, "nlayer",
                            nc_state_file->layer_size,
                            &(nc_state_file->layer_dimid));
        check_nc_status(status, "Error defining nlayer in %s", filename);

        status = nc_def_dim(nc_state_file->nc_id, "soil_node",
                            nc_state_file->node_size,
                            &(nc_state_file->node_dimid));
        check_nc_status(status, "Error defining soil_node in %s", filename);


        // add extension dimensions
        initialize_state_file_rout_extension(filename, nc_state_file);

        set_nc_state_var_info(nc_state_file);
    }

    // initialize dimids to invalid values
    for (i = 0; i < MAXDIMS; i++) {
        dimids[i] = -1;
        dcount[i] = 0;
    }

    // write dimension variables

    // Coordinate variables
    if (mpi_rank == VIC_MPI_ROOT) {
        ndims = global_domain.info.n_coord_dims;
        dstart[0] = 0;
        dstart[1] = 0;

        if (global_domain.info.n_coord_dims == 1) {
            dimids[0] = nc_state_file->ni_dimid;
            dcount[0] = nc_state_file->ni_size;
        }
        else if (global_domain.info.n_coord_dims == 2) {
            dimids[0] = nc_state_file->nj_dimid;
            dcount[0] = nc_state_file->nj_size;

            dimids[1] = nc_state_file->ni_dimid;
            dcount[1] = nc_state_file->ni_size;
        }
        else {
            log_err("COORD_DIMS_OUT should be 1 or 2");
        }
    }

    // define the netcdf variable longitude
    if (mpi_rank == VIC_MPI_ROOT) {
        status = nc_def_var(nc_state_file->nc_id, global_domain.info.lon_var,
                            NC_DOUBLE, ndims, dimids, &(lon_var_id));
        check_nc_status(status, "Error defining lon variable in %s", filename);

        status = nc_put_att_text(nc_state_file->nc_id, lon_var_id, "long_name", strlen(
                                     "longitude"), "longitude");
        check_nc_status(status, "Error adding attribute in %s", filename);
        status = nc_put_att_text(nc_state_file->nc_id, lon_var_id, "units", strlen(
                                     "degrees_east"), "degrees_east");
        check_nc_status(status, "Error adding attribute in %s", filename);
        status = nc_put_att_text(nc_state_file->nc_id, lon_var_id,
                                 "standard_name", strlen(
                                     "longitude"), "longitude");
        check_nc_status(status, "Error adding attribute in %s", filename);

        if (global_domain.info.n_coord_dims == 1) {
            dimids[0] = nc_state_file->nj_dimid;
            dcount[0] = nc_state_file->nj_size;
        }

        // define the netcdf variable latitude
        status = nc_def_var(nc_state_file->nc_id, global_domain.info.lat_var,
                            NC_DOUBLE, ndims, dimids, &(lat_var_id));
        check_nc_status(status, "Error defining lat variable (%s) in %s",
                        global_domain.info.lat_var, filename);
        status = nc_put_att_text(nc_state_file->nc_id, lat_var_id, "long_name", strlen(
                                     "latitude"), "latitude");
        check_nc_status(status, "Error adding attribute in %s", filename);
        status = nc_put_att_text(nc_state_file->nc_id, lat_var_id, "units", strlen(
                                     "degrees_north"), "degrees_north");
        check_nc_status(status, "Error adding attribute in %s", filename);
        status = nc_put_att_text(nc_state_file->nc_id, lat_var_id,
                                 "standard_name", strlen("latitude"),
                                 "latitude");
        check_nc_status(status, "Error adding attribute in %s", filename);
        for (i = 0; i < MAXDIMS; i++) {
            dimids[i] = -1;
            dcount[i] = 0;
        }

        // veg_class
        dimids[0] = nc_state_file->veg_dimid;
        status = nc_def_var(nc_state_file->nc_id, "veg_class",
                            NC_INT, 1, dimids, &(veg_var_id));
        check_nc_status(status, "Error defining veg_class variable in %s",
                        filename);
        status = nc_put_att_text(nc_state_file->nc_id, veg_var_id, "long_name",
                                 strlen("veg_class"), "veg_class");
        check_nc_status(status, "Error adding attribute in %s", filename);
        status = nc_put_att_text(nc_state_file->nc_id, veg_var_id,
                                 "standard_name",
                                 strlen("vegetation_class_number"),
                                 "vegetation_class_number");
        check_nc_status(status, "Error adding attribute in %s", filename);
        dimids[0] = -1;

        // snow_band
        dimids[0] = nc_state_file->band_dimid;
        status = nc_def_var(nc_state_file->nc_id, "snow_band",
                            NC_INT, 1, dimids, &(snow_band_var_id));
        check_nc_status(status, "Error defining snow_band variable in %s",
                        filename);
        status = nc_put_att_text(nc_state_file->nc_id, snow_band_var_id,
                                 "long_name",
                                 strlen("snow_band"), "snow_band");
        check_nc_status(status, "Error adding attribute in %s", filename);
        status = nc_put_att_text(nc_state_file->nc_id, snow_band_var_id,
                                 "standard_name",
                                 strlen("snow_elevation_band_number"),
                                 "snow_elevation_band_number");
        check_nc_status(status, "Error adding attribute in %s", filename);
        dimids[0] = -1;

        // layer
        dimids[0] = nc_state_file->layer_dimid;
        status =
            nc_def_var(nc_state_file->nc_id, "layer", NC_INT, 1, dimids,
                       &(layer_var_id));
        check_nc_status(status, "Error defining layer variable in %s",
                        filename);
        status = nc_put_att_text(nc_state_file->nc_id, layer_var_id,
                                 "long_name",
                                 strlen("layer"), "layer");
        check_nc_status(status, "Error adding attribute in %s", filename);
        status = nc_put_att_text(nc_state_file->nc_id, layer_var_id,
                                 "standard_name", strlen("soil_layer_number"),
                                 "soil_layer_number");
        check_nc_status(status, "Error adding attribute in %s", filename);
        dimids[0] = -1;

        // dz_node (dimension: node, lat, lon)
        dimids[0] = nc_state_file->node_dimid;
        dimids[1] = nc_state_file->nj_dimid;
        dimids[2] = nc_state_file->ni_dimid;
        status = nc_def_var(nc_state_file->nc_id, "dz_node", NC_DOUBLE, 3,
                            dimids, &(dz_node_var_id));
        check_nc_status(status, "Error defining node variable in %s", filename);
        status = nc_put_att_text(nc_state_file->nc_id, dz_node_var_id,
                                 "long_name",
                                 strlen("dz_node"), "dz_node");
        check_nc_status(status, "Error adding attribute in %s", filename);
        status = nc_put_att_text(nc_state_file->nc_id, dz_node_var_id,
                                 "standard_name", strlen(
                                     "soil_thermal_node_spacing"),
                                 "soil_thermal_node_spacing");
        check_nc_status(status, "Error adding attribute in %s", filename);
        status = nc_put_att_text(nc_state_file->nc_id, dz_node_var_id, "units",
                                 strlen("m"), "m");
        check_nc_status(status, "Error adding attribute in %s", filename);
        dimids[0] = -1;
        dimids[1] = -1;
        dimids[2] = -1;

        // node_depth (dimension: node, lat, lon)
        dimids[0] = nc_state_file->node_dimid;
        dimids[1] = nc_state_file->nj_dimid;
        dimids[2] = nc_state_file->ni_dimid;
        status = nc_def_var(nc_state_file->nc_id, "node_depth", NC_DOUBLE, 3,
                            dimids, &(node_depth_var_id));
        check_nc_status(status, "Error defining node variable in %s", filename);
        status = nc_put_att_text(nc_state_file->nc_id, node_depth_var_id,
                                 "long_name",
                                 strlen("node_depth"), "node_depth");
        check_nc_status(status, "Error adding attribute in %s", filename);
        status = nc_put_att_text(nc_state_file->nc_id, node_depth_var_id,
                                 "standard_name",
                                 strlen("soil_thermal_node_depth"),
                                 "soil_thermal_node_depth");
        check_nc_status(status, "Error adding attribute in %s", filename);
        status = nc_put_att_text(nc_state_file->nc_id, node_depth_var_id,
                                 "units",
                                 strlen("m"), "m");
        check_nc_status(status, "Error adding attribute in %s", filename);

    }

    // Define state variables
    if (mpi_rank == VIC_MPI_ROOT) {
        for (i = 0; i < (N_STATE_VARS + N_STATE_VARS_EXT); i++) {
            if (strcasecmp(state_metadata[i].varname, MISSING_S) == 0) {
                // skip variables not set in set_state_meta_data_info
                continue;
            }

            // create the variable
            status = nc_def_var(nc_state_file->nc_id, state_metadata[i].varname,
                                nc_state_file->nc_vars[i].nc_type,
                                nc_state_file->nc_vars[i].nc_dims,
                                nc_state_file->nc_vars[i].nc_dimids,
                                &(nc_state_file->nc_vars[i].nc_varid));
            check_nc_status(status, "Error defining state variable %s in %s",
                            state_metadata[i].varname, filename);

            // set the fill value attribute
            if (nc_state_file->nc_vars[i].nc_type == NC_DOUBLE) {
                status = nc_put_att_double(nc_state_file->nc_id,
                                           nc_state_file->nc_vars[i].nc_varid,
                                           "_FillValue", NC_DOUBLE, 1,
                                           &(nc_state_file->d_fillvalue));
            }
            else if (nc_state_file->nc_vars[i].nc_type == NC_INT) {
                status = nc_put_att_int(nc_state_file->nc_id,
                                        nc_state_file->nc_vars[i].nc_varid,
                                        "_FillValue", NC_INT, 1,
                                        &(nc_state_file->i_fillvalue));
            }
            else {
                log_err("NC_TYPE %d not supported at this time",
                        nc_state_file->nc_vars[i].nc_type);
            }
            check_nc_status(status,
                            "Error putting _FillValue attribute to %s in %s",
                            state_metadata[i].varname, filename);

            // Set string attributes
            put_nc_attr(nc_state_file->nc_id,
                        nc_state_file->nc_vars[i].nc_varid,
                        "long_name", state_metadata[i].long_name);
            put_nc_attr(nc_state_file->nc_id,
                        nc_state_file->nc_vars[i].nc_varid,
                        "standard_name", state_metadata[i].standard_name);
            put_nc_attr(nc_state_file->nc_id,
                        nc_state_file->nc_vars[i].nc_varid,
                        "units", state_metadata[i].units);
            put_nc_attr(nc_state_file->nc_id,
                        nc_state_file->nc_vars[i].nc_varid,
                        "description", state_metadata[i].description);
        }

        // leave define mode
        status = nc_enddef(nc_state_file->nc_id);
        check_nc_status(status, "Error leaving define mode for %s", filename);
    }

    // time variable
    if (mpi_rank == VIC_MPI_ROOT) {
        dtime = date2num(global_param.time_origin_num, dmy_state, 0,
                         global_param.calendar, global_param.time_units);
        // put in netCDF file
        dstart[0] = 0;
        status = nc_put_var1_double(nc_state_file->nc_id,
                                    nc_state_file->time_varid,
                                    dstart, &dtime);
        check_nc_status(status, "Error writing time variable");
    }

    // populate lat/lon
    if (mpi_rank == VIC_MPI_ROOT) {
        if (global_domain.info.n_coord_dims == 1) {
            dvar = calloc(nc_state_file->ni_size, sizeof(*dvar));
            check_alloc_status(dvar, "Memory allocation error");

            dcount[0] = nc_state_file->ni_size;
            // implicitly nested loop over ni and nj with j set to 0
            for (i = 0; i < nc_state_file->ni_size; i++) {
                dvar[i] = (double) global_domain.locations[i].longitude;
            }
            status = nc_put_vara_double(nc_state_file->nc_id, lon_var_id,
                                        dstart,
                                        dcount, dvar);
            check_nc_status(status, "Error adding data to lon in %s", filename);
            free(dvar);

            dvar = calloc(nc_state_file->nj_size, sizeof(*dvar));
            check_alloc_status(dvar, "Memory allocation error");
            dcount[0] = nc_state_file->nj_size;
            // implicitly nested loop over ni and nj with i set to 0;
            // j stride = ni_size
            for (j = 0; j < nc_state_file->nj_size; j++) {
                dvar[j] =
                    (double) global_domain.locations[j *
                                                     nc_state_file->ni_size].
                    latitude;
            }

            status = nc_put_vara_double(nc_state_file->nc_id, lat_var_id,
                                        dstart,
                                        dcount, dvar);
            check_nc_status(status, "Error adding data to lon in %s", filename);
            free(dvar);
        }
        else if (global_domain.info.n_coord_dims == 2) {
            dvar = calloc(global_domain.ncells_total, sizeof(*dvar));
            check_alloc_status(dvar, "Memory allocation error");
            dcount[0] = nc_state_file->nj_size;
            dcount[1] = nc_state_file->ni_size;

            for (i = 0; i < global_domain.ncells_total; i++) {
                dvar[i] = (double) global_domain.locations[i].longitude;
            }
            status = nc_put_vara_double(nc_state_file->nc_id, lon_var_id,
                                        dstart,
                                        dcount, dvar);
            check_nc_status(status, "Error adding data to lon in %s", filename);

            for (i = 0; i < global_domain.ncells_total;
                 i++) {
                dvar[i] = (double) global_domain.locations[i].latitude;
            }
            status = nc_put_vara_double(nc_state_file->nc_id, lat_var_id,
                                        dstart,
                                        dcount, dvar);
            check_nc_status(status, "Error adding data to lat in %s", filename);

            free(dvar);
        }
        else {
            log_err("COORD_DIMS_OUT should be 1 or 2");
        }
    }

    // Variables for other dimensions (all 1-dimensional)
    if (mpi_rank == VIC_MPI_ROOT) {
        ndims = 1;

        // vegetation classes
        dimids[0] = nc_state_file->veg_dimid;
        dcount[0] = nc_state_file->veg_size;
        ivar = malloc(nc_state_file->veg_size * sizeof(*ivar));
        check_alloc_status(ivar, "Memory allocation error");

        for (j = 0; j < nc_state_file->veg_size; j++) {
            ivar[j] = (int) j + 1;
        }
        status = nc_put_vara_int(nc_state_file->nc_id, veg_var_id, dstart,
                                 dcount,
                                 ivar);
        check_nc_status(status, "Error writing veg var id");
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
            dcount[i] = 0;
        }
        free(ivar);

        // snow bands
        dimids[0] = nc_state_file->band_dimid;
        dcount[0] = nc_state_file->band_size;
        ivar = malloc(nc_state_file->band_size * sizeof(*ivar));
        check_alloc_status(ivar, "Memory allocation error");

        for (j = 0; j < nc_state_file->band_size; j++) {
            ivar[j] = (int) j;
        }
        status = nc_put_vara_int(nc_state_file->nc_id, snow_band_var_id, dstart,
                                 dcount, ivar);
        check_nc_status(status, "Error writing snow band id");
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
            dcount[i] = 0;
        }
        free(ivar);

        // soil layers
        dimids[0] = nc_state_file->layer_dimid;
        dcount[0] = nc_state_file->layer_size;
        ivar = malloc(nc_state_file->layer_size * sizeof(*ivar));
        check_alloc_status(ivar, "Memory allocation error");

        for (j = 0; j < nc_state_file->layer_size; j++) {
            ivar[j] = (int) j;
        }
        status = nc_put_vara_int(nc_state_file->nc_id, layer_var_id, dstart,
                                 dcount,
                                 ivar);
        check_nc_status(status, "Error writing layer id");
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
            dcount[i] = 0;
        }
        free(ivar);
    }

    // initialize dvar for soil thermal node deltas and depths
    dvar = malloc(local_domain.ncells_active * sizeof(*dvar));
    check_alloc_status(dvar, "Memory allocation error");
    // set missing values
    for (i = 0; i < local_domain.ncells_active; i++) {
        dvar[i] = nc_state_file->d_fillvalue;
    }

    // soil thermal node deltas (dimension: node, lat, lon)
    dstart[0] = 0;
    dstart[1] = 0;
    dstart[2] = 0;
    dcount[0] = 1;
    dcount[1] = nc_state_file->nj_size;
    dcount[2] = nc_state_file->ni_size;
    for (j = 0; j < nc_state_file->node_size; j++) {
        dstart[0] = j;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = soil_con[i].dz_soil[j];
        }
        gather_put_nc_field_double(nc_state_file->nc_id,
                                   dz_node_var_id,
                                   nc_state_file->d_fillvalue,
                                   dstart, dcount, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file->d_fillvalue;
        }
    }

    // soil thermal node depths (dimension: node, lat, lon)
    dstart[0] = 0;
    dstart[1] = 0;
    dstart[2] = 0;
    dcount[0] = 1;
    dcount[1] = nc_state_file->nj_size;
    dcount[2] = nc_state_file->ni_size;
    for (j = 0; j < nc_state_file->node_size; j++) {
        dstart[0] = j;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = soil_con[i].Zsum_soil[j];
        }
        gather_put_nc_field_double(nc_state_file->nc_id,
                                   node_depth_var_id,
                                   nc_state_file->d_fillvalue,
                                   dstart, dcount, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file->d_fillvalue;
        }
    }
    for (i = 0; i < MAXDIMS; i++) {
        dimids[i] = -1;
        dcount[i] = 0;
    }
    free(dvar);

}
