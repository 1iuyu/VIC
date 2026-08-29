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
    extern rout_struct        *rout;
    extern domain_struct       global_domain;
    extern domain_struct       local_domain;
    extern option_struct       options;
    extern filenames_struct    filenames;
    extern metadata_struct     state_metadata[N_STATE_VARS];

    size_t                     i, j, m;
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
                if (m < local_domain.locations[i].nveg &&
                            j < all_vars[i].cell[m].Nsoil) {
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
                if (m < local_domain.locations[i].nveg &&
                            j < all_vars[i].cell[m].Nsoil) {
                    all_vars[i].cell[m].ice[j] = dvar[i];
                }
            }
        }
    }

    // liq content
    for (m = 0; m < MAX_HRUS; m++) {
        d4start[0] = m;
        for (j = 0; j < MAX_SOILS; j++) {
            d4start[1] = j;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_SOIL_LIQ].varname,
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg &&
                            j < all_vars[i].cell[m].Nsoil) {
                    all_vars[i].cell[m].liq[j] = dvar[i];
                }
            }
        }
    }

    // last step ice content
    for (m = 0; m < MAX_HRUS; m++) {
        d4start[0] = m;
        for (j = 0; j < MAX_SOILS; j++) {
            d4start[1] = j;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_SOIL_LASTICE].varname,
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg &&
                            j < all_vars[i].cell[m].Nsoil) {
                    all_vars[i].cell[m].last_ice[j] = dvar[i];
                }
            }
        }
    }

    // last step water content
    for (m = 0; m < MAX_HRUS; m++) {
        d4start[0] = m;
        for (j = 0; j < MAX_SOILS; j++) {
            d4start[1] = j;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_SOIL_LASTLIQ].varname,
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg &&
                            j < all_vars[i].cell[m].Nsoil) {
                    all_vars[i].cell[m].last_liq[j] = dvar[i];
                }
            }
        }
    }

    // lastswq: tmpval = snow[veg].lastswq;
    for (m = 0; m < MAX_HRUS; m++) {
        d3start[0] = m;
        get_scatter_nc_field_double(&(filenames.init_state),
                                    state_metadata[STATE_SNOW_LASTSWQ].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                all_vars[i].snow[m].last_swq = dvar[i];
            }
        }
    }

    // int_snow: tmpval = veg_var[veg].int_snow;
    for (m = 0; m < MAX_HRUS; m++) {
        d3start[0] = m;
        get_scatter_nc_field_double(&(filenames.init_state),
                                    state_metadata[STATE_INT_SNOW].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                all_vars[i].veg_var[m].int_snow = dvar[i];
            }
        }
    }

    // int_rain: tmpval = veg_var[veg].int_rain;
    for (m = 0; m < MAX_HRUS; m++) {
        d3start[0] = m;
        get_scatter_nc_field_double(&(filenames.init_state),
                                    state_metadata[STATE_INT_RAIN].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                all_vars[i].veg_var[m].int_rain = dvar[i];
            }
        }
    }

    // snow age: snow[veg].snowage
    for (m = 0; m < MAX_HRUS; m++) {
        d3start[0] = m;
        get_scatter_nc_field_double(&(filenames.init_state),
                                    state_metadata[STATE_SNOW_AGE].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                all_vars[i].snow[m].snowage = dvar[i];
            }
        }
    }

    // snow covered fraction: snow[veg].coverage
    for (m = 0; m < MAX_HRUS; m++) {
        d3start[0] = m;
        get_scatter_nc_field_double(&(filenames.init_state),
                                    state_metadata[STATE_SNOW_COVERAGE].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                all_vars[i].snow[m].coverage = dvar[i];
            }
        }
    }

    // pack_ice: snow[veg].pack_ice
    for (m = 0; m < MAX_HRUS; m++) {
        d4start[0] = m;
        for (j = 0; j < MAX_SNOWS; j++) {
            d4start[1] = j;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_SNOW_PACK_ICE].varname,
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg && 
                                all_vars[i].snow[m].Nsnow > 0) {
                    all_vars[i].snow[m].pack_ice[j] = dvar[i];
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

    // snow last_thice: snow[veg].last_thice
    for (m = 0; m < MAX_HRUS; m++) {
        d4start[0] = m;
        for (j = 0; j < MAX_SNOWS; j++) {
            d4start[1] = j;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_SNOW_LASTICE].varname,
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg &&
                                all_vars[i].snow[m].Nsnow > 0) {
                    all_vars[i].snow[m].last_thice[j] = dvar[i];
                }
            }
        }
    }

    // snow last theta liq: snow[veg].last_thliq
    for (m = 0; m < MAX_HRUS; m++) {
        d4start[0] = m;
        for (j = 0; j < MAX_SNOWS; j++) {
            d4start[1] = j;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_SNOW_LASTLIQ].varname,
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg &&
                                all_vars[i].snow[m].Nsnow > 0) {
                    all_vars[i].snow[m].last_thliq[j] = dvar[i];
                }
            }
        }
    }

    // snow radius: snow[veg].radius
    for (m = 0; m < MAX_HRUS; m++) {
        d4start[0] = m;
        for (j = 0; j < MAX_SNOWS; j++) {
            d4start[1] = j;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_SNOW_RADIUS].varname,
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg &&
                                all_vars[i].snow[m].Nsnow > 0) {
                    all_vars[i].snow[m].radius[j] = dvar[i];
                }
            }
        }
    }

    // snow dz_snow snow[veg].dz_snow
    for (m = 0; m < MAX_HRUS; m++) {
        d4start[0] = m;
        for (j = 0; j < MAX_SNOWS; j++) {
            d4start[1] = j;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_SNOW_DZNODE].varname,
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg &&
                                all_vars[i].snow[m].Nsnow > 0) {
                    all_vars[i].snow[m].dz_snow[j] = dvar[i];
                }
            }
        }
    }

    // snow pack water: snow[veg].pack_liq
    for (m = 0; m < MAX_HRUS; m++) {
        d4start[0] = m;
        for (j = 0; j < MAX_SNOWS; j++) {
            d4start[1] = j;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_SNOW_PACK_LIQ].varname,
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg &&
                                all_vars[i].snow[m].Nsnow > 0) {
                    all_vars[i].snow[m].pack_liq[j] = dvar[i];
                }
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
                                j < all_vars[i].cell[m].Nnode) {
                    all_vars[i].energy[m].T[j] = dvar[i];
                }
            }
        }
    }

    // thermal node temperatures: energy[veg].last_T[nidx]
    for (m = 0; m < MAX_HRUS; m++) {
        d4start[0] = m;
        for (j = 0; j < MAX_NODES; j++) {
            d4start[1] = j;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_LAST_TEMP].varname,
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg &&
                                j < all_vars[i].cell[m].Nnode) {
                    all_vars[i].energy[m].last_T[j] = dvar[i];
                }
            }
        }
    }

    // cell[veg].matric[nidx]
    for (m = 0; m < MAX_HRUS; m++) {
        d4start[0] = m;
        for (j = 0; j < MAX_SOILS; j++) {
            d4start[1] = j;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_MATRIC].varname,
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg &&
                                j < all_vars[i].cell[m].Nsoil) {
                    all_vars[i].cell[m].matric[j] = dvar[i];
                }
            }
        }
    }

    // cell[veg].last_matric[nidx]
    for (m = 0; m < MAX_HRUS; m++) {
        d4start[0] = m;
        for (j = 0; j < MAX_SOILS; j++) {
            d4start[1] = j;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_LAST_MATRIC].varname,
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg &&
                                j < all_vars[i].cell[m].Nsoil) {
                    all_vars[i].cell[m].last_matric[j] = dvar[i];
                }
            }
        }
    }

    // veg_var[veg].mat_VEG[nidx]
    for (m = 0; m < MAX_HRUS; m++) {
        d4start[0] = m;
        for (j = 0; j < 4; j++) {
            d4start[1] = j;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_VEG_MATRIC].varname,
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg) {
                    all_vars[i].veg_var[m].mat_VEG[j] = dvar[i];
                }
            }
        }
    }

    // zwt
    for (m = 0; m < MAX_HRUS; m++) {
        d3start[0] = m;
        get_scatter_nc_field_double(&(filenames.init_state),
                                    state_metadata[STATE_ZWT].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                all_vars[i].cell[m].zwt = dvar[i];
            }
        }
    }

    // aqf_storage
    for (m = 0; m < MAX_HRUS; m++) {
        d3start[0] = m;
        get_scatter_nc_field_double(&(filenames.init_state),
                                    state_metadata[STATE_AQF_STORAGE].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                all_vars[i].cell[m].storage_aqf = dvar[i];
            }
        }
    }

    // nsnow
    for (m = 0; m < MAX_HRUS; m++) {
        d3start[0] = m;
        get_scatter_nc_field_int(&(filenames.init_state),
                                state_metadata[STATE_NSNOW].varname,
                                d3start, d3count, ivar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                all_vars[i].snow[m].Nsnow = ivar[i];
            }
        }
    }

    // nsoil
    for (m = 0; m < MAX_HRUS; m++) {
        d3start[0] = m;
        get_scatter_nc_field_int(&(filenames.init_state),
                                state_metadata[STATE_NSOIL].varname,
                                d3start, d3count, ivar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                all_vars[i].cell[m].Nsoil = ivar[i];
            }
        }
    }

    // nnode
    for (m = 0; m < MAX_HRUS; m++) {
        d3start[0] = m;
        get_scatter_nc_field_int(&(filenames.init_state),
                                state_metadata[STATE_NNODE].varname,
                                d3start, d3count, ivar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                all_vars[i].cell[m].Nnode = ivar[i];
            }
        }
    }

    // ncanopy
    for (m = 0; m < MAX_HRUS; m++) {
        d3start[0] = m;
        get_scatter_nc_field_int(&(filenames.init_state),
                                state_metadata[STATE_NCANOPY].varname,
                                d3start, d3count, ivar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                all_vars[i].veg_var[m].Ncanopy = ivar[i];
            }
        }
    }

    // nroot
    for (m = 0; m < MAX_HRUS; m++) {
        d3start[0] = m;
        get_scatter_nc_field_int(&(filenames.init_state),
                                state_metadata[STATE_NROOT].varname,
                                d3start, d3count, ivar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                all_vars[i].veg_var[m].Nroot = ivar[i];
            }
        }
    }

    // h2osfc
    for (m = 0; m < MAX_HRUS; m++) {
        d3start[0] = m;
        get_scatter_nc_field_double(&(filenames.init_state),
                                    state_metadata[STATE_H2OSFC].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                all_vars[i].cell[m].h2osfc = dvar[i];
            }
        }
    }

    // h2o_frac
    for (m = 0; m < MAX_HRUS; m++) {
        d3start[0] = m;
        get_scatter_nc_field_double(&(filenames.init_state),
                                    state_metadata[STATE_H2O_FRAC].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                all_vars[i].cell[m].frac_h2o = dvar[i];
            }
        }
    }

    // h2osfc_ice
    for (m = 0; m < MAX_HRUS; m++) {
        d3start[0] = m;
        get_scatter_nc_field_double(&(filenames.init_state),
                                    state_metadata[STATE_H2OSFC_ICE].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                all_vars[i].cell[m].h2osfc_ice = dvar[i];
            }
        }
    }

    // h2osfc_liq
    for (m = 0; m < MAX_HRUS; m++) {
        d3start[0] = m;
        get_scatter_nc_field_double(&(filenames.init_state),
                                    state_metadata[STATE_H2OSFC_LIQ].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (m < local_domain.locations[i].nveg) {
                all_vars[i].cell[m].h2osfc_liq = dvar[i];
            }
        }
    }

    // rout state
    if (options.ROUT) {
        // main_channel_storage
        for (m = 0; m < MAX_HRUS; m++) {
            d3start[0] = m;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_MAIN_CHANNEL_STORAGE].varname,
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg) {
                    rout[i].main_channel.wr = dvar[i];
                }
            }
        }

        // main_cross_section_area
        for (m = 0; m < MAX_HRUS; m++) {
            d3start[0] = m;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_MAIN_CROSS_SECTION_AREA].varname,
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg) {
                    rout[i].main_channel.mr = dvar[i];
                }
            }
        }

        // main_channel_depth
        for (m = 0; m < MAX_HRUS; m++) {
            d3start[0] = m;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_MAIN_CHANNEL_DEPTH].varname,
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg) {
                    rout[i].main_channel.yr = dvar[i];
                }
            }
        }

        // main_channel_manning_n
        for (m = 0; m < MAX_HRUS; m++) {
            d3start[0] = m;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_MAIN_CHANNEL_MANNING_N].varname,
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg) {
                    rout[i].main_channel.nr = dvar[i];
                }
            }
        }

        // main_wetted_perimeter
        for (m = 0; m < MAX_HRUS; m++) {
            d3start[0] = m;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_MAIN_WETTED_PERIMETER].varname,
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg) {
                    rout[i].main_channel.pr = dvar[i];
                }
            }
        }

        // main_hydraulic_radius
        for (m = 0; m < MAX_HRUS; m++) {
            d3start[0] = m;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_MAIN_HYDRAULIC_RADIUS].varname,
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg) {
                    rout[i].main_channel.rr = dvar[i];
                }
            }
        }

        // sub_channel_storage
        for (m = 0; m < MAX_HRUS; m++) {
            d3start[0] = m;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_SUB_CHANNEL_STORAGE].varname,
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg) {
                    rout[i].sub_channel.wt = dvar[i];
                }
            }
        }

        // sub_channel_manning_n
        for (m = 0; m < MAX_HRUS; m++) {
            d3start[0] = m;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_SUB_CHANNEL_MANNING_N].varname,
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg) {
                    rout[i].sub_channel.nt = dvar[i];
                }
            }
        }

        // sub_cross_section_area
        for (m = 0; m < MAX_HRUS; m++) {
            d3start[0] = m;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_SUB_CROSS_SECTION_AREA].varname,
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg) {
                    rout[i].sub_channel.mt = dvar[i];
                }
            }
        }

        // sub_channel_depth
        for (m = 0; m < MAX_HRUS; m++) {
            d3start[0] = m;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_SUB_CHANNEL_DEPTH].varname,
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg) {
                    rout[i].sub_channel.yt = dvar[i];
                }
            }
        }

        // sub_wetted_perimeter
        for (m = 0; m < MAX_HRUS; m++) {
            d3start[0] = m;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_SUB_WETTED_PERIMETER].varname,
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg) {
                    rout[i].sub_channel.pt = dvar[i];
                }
            }
        }

        // sub_hydraulic_radius
        for (m = 0; m < MAX_HRUS; m++) {
            d3start[0] = m;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_SUB_HYDRAULIC_RADIUS].varname,
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg) {
                    rout[i].sub_channel.rt = dvar[i];
                }
            }
        }

        // hillslope_depth
        for (m = 0; m < MAX_HRUS; m++) {
            d3start[0] = m;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_HILLSLOPE_DEPTH].varname,
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg) {
                    rout[i].hillslope.yh = dvar[i];
                }
            }
        }

        // hillslope_manning_n
        for (m = 0; m < MAX_HRUS; m++) {
            d3start[0] = m;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_HILLSLOPE_MANNING_N].varname,
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg) {
                    rout[i].hillslope.nh = dvar[i];
                }
            }
        }

        // hillslope_storage
        for (m = 0; m < MAX_HRUS; m++) {
            d3start[0] = m;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_HILLSLOPE_STORAGE].varname,
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg) {
                    rout[i].hillslope.wh = dvar[i];
                }
            }
        }

        // storage_prev
        for (m = 0; m < MAX_HRUS; m++) {
            d3start[0] = m;
            get_scatter_nc_field_double(&(filenames.init_state),
                                        state_metadata[STATE_STORAGE_PREV].varname,
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (m < local_domain.locations[i].nveg) {
                    rout[i].total_storage_prev = dvar[i];
                }
            }
        }
    }

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
        dimlen = get_nc_dimension(&(filenames.init_state), "nveg");
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
        dimlen = get_nc_dimension(&(filenames.init_state), "nsoil");
        if (dimlen != MAX_SOILS) {
            log_err("Number of soil layer in state file does not "
                    "match parameter file");
        }
        dimlen = get_nc_dimension(&(filenames.init_state), "nsnow");
        if (dimlen != MAX_SNOWS) {
            log_err("Number of snow layer in state file does not "
                    "match parameter file");
        }
        dimlen = get_nc_dimension(&(filenames.init_state), "nwave");
        if (dimlen != MAX_SWBANDS) {
            log_err("Number of solar wave in state file does not "
                    "match parameter file");
        }
        dimlen = get_nc_dimension(&(filenames.init_state), "ncanopy");
        if (dimlen != MAX_CANOPYS) {
            log_err("Number of canopy layer in state file does not "
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
