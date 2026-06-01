/******************************************************************************
 * @section DESCRIPTION
 *
 * Initialize model parameters
 *****************************************************************************/

#include "vic_driver_shared_image.h"

/******************************************************************************
 * @brief    Initialize model parameters
 *****************************************************************************/
void
vic_init(void)
{
    extern all_vars_struct    *all_vars;
    extern size_t              current;
    extern domain_struct       global_domain;
    extern domain_struct       local_domain;
    extern option_struct       options;
    extern filenames_struct    filenames;
    extern soil_con_struct    *soil_con;
    extern veg_con_struct    **veg_con;
    extern veg_lib_struct     *veg_lib;
    extern parameters_struct   param;
    extern int                 mpi_rank;

    bool                       found;
    char                       locstr[MAXSTRING];
    double                     mean;
    double                     sum;
    double                    *temp_array = NULL;
    double                    *Cv_sum = NULL;
    double                    *dvar = NULL;
    int                       *ivar = NULL;
    int                        status;
    size_t                     i,j,k,m;
    size_t                     nveg;
    size_t                     lidx;
    int                        vidx;
    size_t                     d2count[2];
    size_t                     d2start[2];
    size_t                     d3count[3];
    size_t                     d3start[3];
    size_t                     n_ny, n_nx;

    // allocate memory for Cv_sum
    Cv_sum = malloc(local_domain.ncells_active * sizeof(*Cv_sum));
    check_alloc_status(Cv_sum, "Memory allocation error.");

    // allocate memory for variables to be read
    dvar = malloc(local_domain.ncells_active * sizeof(*dvar));
    check_alloc_status(dvar, "Memory allocation error.");
    ivar = malloc(local_domain.ncells_active * sizeof(*ivar));
    check_alloc_status(ivar, "Memory allocation error.");
    int total_size = local_domain.ncells_active * options.Nlayer;
    temp_array = (double*) malloc(total_size * sizeof(*temp_array));
    check_alloc_status(temp_array, "Memory allocation error.");
    
    // The method used to convert the NetCDF fields to VIC structures for
    // individual grid cells is to read a 2D slice and then loop over the
    // domain cells to assign the values to the VIC structures
    n_ny = global_domain.n_ny;
    n_nx = global_domain.n_nx;

    d2start[0] = 0;
    d2start[1] = 0;
    d2count[0] = n_ny;
    d2count[1] = n_nx;

    d3start[0] = 0;
    d3start[1] = 0;
    d3start[2] = 0;
    d3count[0] = 1;
    d3count[1] = n_ny;
    d3count[2] = n_nx;

    current = 0;

    /*********************************************
       Reading the soil parameters for each grid
    *********************************************/

    // read_soilparam()
    
    // Validate Nlayer
    if (options.Nlayer < 1) {
        log_err("You must define at least 1 soil moisture layer to run "
                "the model.  Currently Nlayers is set to %zu.",
                options.Nlayer);
    }
    if (options.Nlayer > MAX_LAYERS) {
        log_err("Global file wants more soil moisture layers (%zu) than "
                "are defined by MAX_LAYERS (%d).  Edit vic_driver_shared.h and "
                "recompile.", options.Nlayer, MAX_LAYERS);
    }

    // latitude and longitude
    for (i = 0; i < local_domain.ncells_active; i++) {
        soil_con[i].lat = local_domain.locations[i].latitude;
        soil_con[i].lng = local_domain.locations[i].longitude;
    }

    // b_infilt
    get_scatter_nc_field_double(&(filenames.params), "infilt",
                                d2start, d2count, dvar);
    for (i = 0; i < local_domain.ncells_active; i++) {
        soil_con[i].b_infilt = (double) dvar[i];
    }

    // b_dynamic: heterogeniety parameter for infiltration
    get_scatter_nc_field_double(&(filenames.params), "b_dynamic",
                                d2start, d2count, dvar);
    for (i = 0; i < local_domain.ncells_active; i++) {
        soil_con[i].b_dynamic = (double) dvar[i];
    }

    // cap_drive: mean capilary drive (m) for dynamic VIC runoff
    get_scatter_nc_field_double(&(filenames.params), "capil_drive",
                                d2start, d2count, dvar);
    for (i = 0; i < local_domain.ncells_active; i++) {
        soil_con[i].capil_drive = (double) dvar[i];
    }

    // elevation: mean grid cell elevation
    get_scatter_nc_field_double(&(filenames.params), "elev",
                                d2start, d2count, dvar);
    for (i = 0; i < local_domain.ncells_active; i++) {
        soil_con[i].elevation = (double) dvar[i];
    }

    // depth: thickness for each soil layer
    for (j = 0; j < options.Nlayer; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(&(filenames.params), "depth",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            soil_con[i].depth[j] = (double) dvar[i];
        }
    }

    // off_gmt: cell gmt offset
    get_scatter_nc_field_double(&(filenames.params), "off_gmt",
                                d2start, d2count, dvar);
    for (i = 0; i < local_domain.ncells_active; i++) {
        soil_con[i].off_gmt = (double) dvar[i];
    }

    // slope: slope of the cell []
    get_scatter_nc_field_double(&(filenames.params), "slope",
                                d2start, d2count, dvar);
    for (i = 0; i < local_domain.ncells_active; i++) {
        soil_con[i].slope = (double) dvar[i];
    }

    // z_bedrock: Depth to bedrock [m]
    get_scatter_nc_field_double(&(filenames.params), "z_bedrock",
                                d2start, d2count, dvar);
    for (i = 0; i < local_domain.ncells_active; i++) {
        soil_con[i].z_bedrock = (double) dvar[i];
    }

    // Additional processing of the soil variables
    // (compute derived parameters)
    for (i = 0; i < local_domain.ncells_active; i++) {
        for (j = 0; j < options.Nlayer; j++) {
            // check layer thicknesses
            if (soil_con[i].depth[j] < MINSOILDEPTH) {
                sprint_location(locstr, &(local_domain.locations[i]));
                log_err("Model will not function with layer %zd depth %f < %f "
                        "m.\n%s", j, soil_con[i].depth[j], MINSOILDEPTH,
                        locstr);
            }
        }
        // check relative thickness of top two layers
        if (soil_con[i].depth[0] > soil_con[i].depth[1]) {
            sprint_location(locstr, &(local_domain.locations[i]));
            log_err("Model will not function with layer 0 thicker than layer "
                    "(%f m > %f m).\n%s", soil_con[i].depth[0],
                    soil_con[i].depth[1], locstr);
        }
        double z_bedrock = soil_con[i].z_bedrock;
        // Soil thermal node thicknesses and positions
        // 表层2m土层厚度固定
        soil_con[i].dz_soil[0] = 0.02;
        soil_con[i].dz_soil[1] = 0.04;
        soil_con[i].dz_soil[2] = 0.04;
        soil_con[i].dz_soil[3] = 0.1;
        soil_con[i].dz_soil[4] = 0.1;
        soil_con[i].dz_soil[5] = 0.1;
        soil_con[i].dz_soil[6] = 0.1;
        soil_con[i].dz_soil[7] = 0.3;
        soil_con[i].dz_soil[8] = 0.3;
        soil_con[i].dz_soil[9] = 0.3;
        soil_con[i].dz_soil[10] = 0.3;
        soil_con[i].dz_soil[11] = 0.3;
        // 剩余深度以等比数列累加，直到总厚度达到z_bedrock为止
        k = 12;
        double ratio = 1.30;
        double accum_depth = 2.0;
        double current_thickness = 0.5;
        while (k < MAX_SOILS && accum_depth < z_bedrock) {
            /* 检查加上当前层是否会超过总深度 */
            if (accum_depth + current_thickness < z_bedrock) {
                /* 正常分配 */
                soil_con[i].dz_soil[k] = current_thickness;
                accum_depth += current_thickness;
                current_thickness *= ratio;  /* 下一层厚度增长 */
            } else {
                /* 最后一层：设为基岩层 */
                soil_con[i].dz_soil[k] = current_thickness;
                accum_depth = z_bedrock;
                break;
            }
            k++;
            /* 防止厚度增长过大 */
            if (current_thickness > 5.0) {  /* 最大层厚限制 */
                current_thickness = 5.0;
            }
        }
        size_t Nbedrock = k + 1;
        soil_con[i].Nbedrock = Nbedrock;
        /* Compute soil node depths */
        double sum_dz = 0.;
        for (k = 0; k < Nbedrock; k++) {
            sum_dz += soil_con[i].dz_soil[k];
            soil_con[i].Zsum_soil[k] = sum_dz;
        }
        for (k = 0; k < Nbedrock; k++) {
            soil_con[i].zc_soil[k] = soil_con[i].Zsum_soil[k] - 
                                    soil_con[i].dz_soil[k] / 2.;
        }

        // Calculate grid cell area.
        soil_con[i].cell_area = local_domain.locations[i].area;

        // set soil albedo parameters
        soil_con[i].AlbedoSat[0] = 0.09;  // saturated soil albedo at visible band
        soil_con[i].AlbedoSat[1] = 0.18;  // saturated soil albedo at NIR band
        soil_con[i].AlbedoDry[0] = 0.18;  // dry soil albedo at visible band
        soil_con[i].AlbedoDry[1] = 0.36;  // dry soil albedo at NIR band

        /* Central Longitude of Current Time Zone */ 
        soil_con[i].time_zone_lng = soil_con[i].off_gmt * 360. / HOURS_PER_DAY;

    }
    size_t Nlayer = options.Nlayer;
    // clay: clay content for each soil layer
    for (j = 0; j < Nlayer; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(&(filenames.params), "clay",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            lidx = i * Nlayer + j;
            temp_array[lidx] = (double) dvar[i];
        }
    }
    for (i = 0; i < local_domain.ncells_active; i++) {
        set_node_parameters(soil_con[i].Nbedrock,
                            soil_con[i].depth,
                            soil_con[i].Zsum_soil,
                            &soil_con[i].clay_node[0],
                            &temp_array[i * Nlayer]);
    }

    // sand: sand content for each soil layer
    for (j = 0; j < Nlayer; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(&(filenames.params), "sand",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            lidx = i * Nlayer + j;
            temp_array[lidx] = (double) dvar[i];
        }
    }
    for (i = 0; i < local_domain.ncells_active; i++) {
        set_node_parameters(soil_con[i].Nbedrock,
                            soil_con[i].depth,
                            soil_con[i].Zsum_soil,
                            &soil_con[i].sand_node[0],
                            &temp_array[i * Nlayer]);
    }

    // silt: silt content for each soil layer
    for (j = 0; j < Nlayer; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(&(filenames.params), "silt",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            lidx = i * Nlayer + j;
            temp_array[lidx] = (double) dvar[i];
        }
    }
    for (i = 0; i < local_domain.ncells_active; i++) {
        set_node_parameters(soil_con[i].Nbedrock,
                            soil_con[i].depth,
                            soil_con[i].Zsum_soil,
                            &soil_con[i].silt_node[0],
                            &temp_array[i * Nlayer]);
    }

    // gravel: gravel content for each soil layer
    for (j = 0; j < Nlayer; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(&(filenames.params), "gravel",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            lidx = i * Nlayer + j;
            temp_array[lidx] = (double) dvar[i];
        }
    }
    for (i = 0; i < local_domain.ncells_active; i++) {
        set_node_parameters(soil_con[i].Nbedrock,
                            soil_con[i].depth,
                            soil_con[i].Zsum_soil,
                            &soil_con[i].gravel_node[0],
                            &temp_array[i * Nlayer]);
    }

    // Wsat: saturated point for each layer
    for (j = 0; j < Nlayer; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(&(filenames.params), "organic",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            lidx = i * Nlayer + j;
            temp_array[lidx] = (double) dvar[i];
        }
    }
    for (i = 0; i < local_domain.ncells_active; i++) {
        set_node_parameters(soil_con[i].Nbedrock,
                            soil_con[i].depth,
                            soil_con[i].Zsum_soil,
                            &soil_con[i].organic_node[0],
                            &temp_array[i * Nlayer]);
    }

    for (j = 0; j < Nlayer; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(&(filenames.params), "bulk_density",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            lidx = i * Nlayer + j;
            temp_array[lidx] = (double) dvar[i];
        }
    }
    for (i = 0; i < local_domain.ncells_active; i++) {
        set_node_parameters(soil_con[i].Nbedrock,
                            soil_con[i].depth,
                            soil_con[i].Zsum_soil,
                            &soil_con[i].bulk_dens_node[0],
                            &temp_array[i * Nlayer]);
    }
    
    /*  bulk density and soil density (particle density) read from 
        soil parameter file; otherwise campute from PedoTransfer.c */
    if (options.PARAM_FROM_SOIL) {
        // bulk_density: mineral bulk density for each soil layer

        // alpha: retention shape parameter in van Genuchten equation [1/m]
        for (j = 0; j < Nlayer; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(&(filenames.params), "alpha",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                lidx = i * Nlayer + j;
                temp_array[lidx] = (double) dvar[i];
            }
        }
        for (i = 0; i < local_domain.ncells_active; i++) {
            set_node_parameters(soil_con[i].Nbedrock,
                                soil_con[i].depth,
                                soil_con[i].Zsum_soil,
                                &soil_con[i].alpha_node[0],
                                &temp_array[i * Nlayer]);
        }
        // bexp: layer-specific exponent n (=3+2/lambda) in Campbell or van Genuchten eqn
        for (j = 0; j < Nlayer; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(&(filenames.params), "bexp",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                lidx = i * Nlayer + j;
                temp_array[lidx] = (double) dvar[i];
            }
        }
        for (i = 0; i < local_domain.ncells_active; i++) {
            set_node_parameters(soil_con[i].Nbedrock,
                                soil_con[i].depth,
                                soil_con[i].Zsum_soil,
                                &soil_con[i].bexp_node[0],
                                &temp_array[i * Nlayer]);
        }
        // Ksat: saturated hydraulic conductivity [m/s]
        for (j = 0; j < Nlayer; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(&(filenames.params), "Ksat",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                lidx = i * Nlayer + j;
                temp_array[lidx] = (double) dvar[i];
            }
        }
        for (i = 0; i < local_domain.ncells_active; i++) {
            set_node_parameters(soil_con[i].Nbedrock,
                                soil_con[i].depth,
                                soil_con[i].Zsum_soil,
                                &soil_con[i].Ksat_node[0],
                                &temp_array[i * Nlayer]);
        }
        // psisat: soil matric potential at saturation [m]
        for (j = 0; j < Nlayer; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(&(filenames.params), "psisat",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                lidx = i * Nlayer + j;
                temp_array[lidx] = (double) dvar[i];
            }
        }
        for (i = 0; i < local_domain.ncells_active; i++) {
            set_node_parameters(soil_con[i].Nbedrock,
                                soil_con[i].depth,
                                soil_con[i].Zsum_soil,
                                &soil_con[i].psisat_node[0],
                                &temp_array[i * Nlayer]);
        }
        // Wpwp: soil moisture content at permanent wilting point [m3/m3]
        for (j = 0; j < Nlayer; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(&(filenames.params), "Wpwp",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                lidx = i * Nlayer + j;
                temp_array[lidx] = (double) dvar[i];
            }
        }
        for (i = 0; i < local_domain.ncells_active; i++) {
            set_node_parameters(soil_con[i].Nbedrock,
                                soil_con[i].depth,
                                soil_con[i].Zsum_soil,
                                &soil_con[i].Wpwp_node[0],
                                &temp_array[i * Nlayer]);
        }
        // Wsat: soil moisture content at saturation [m3/m3]
        for (j = 0; j < Nlayer; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(&(filenames.params), "Wsat",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                lidx = i * Nlayer + j;
                temp_array[lidx] = (double) dvar[i];
            }
        }
        for (i = 0; i < local_domain.ncells_active; i++) {
            set_node_parameters(soil_con[i].Nbedrock,
                                soil_con[i].depth,
                                soil_con[i].Zsum_soil,
                                &soil_con[i].Wsat_node[0],
                                &temp_array[i * Nlayer]);
        }
        // lpar: unsaturated hydraulic conductivity exponent.
        for (j = 0; j < Nlayer; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(&(filenames.params), "lpar",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                lidx = i * Nlayer + j;
                temp_array[lidx] = (double) dvar[i];
            }
        }
        for (i = 0; i < local_domain.ncells_active; i++) {
            set_node_parameters(soil_con[i].Nbedrock,
                                soil_con[i].depth,
                                soil_con[i].Zsum_soil,
                                &soil_con[i].lpar_node[0],
                                &temp_array[i * Nlayer]);
        }
    }

    /******************************************
       Compute snow band elevaton properties
    ******************************************/
    // read_snowband()
    if (options.SNOW_BAND == 1) {
        for (i = 0; i < local_domain.ncells_active; i++) {
            soil_con[i].AreaFract[0] = 1.;
            soil_con[i].BandElev[0] = soil_con[i].elevation;
            soil_con[i].Pfactor[0] = 1.;
            soil_con[i].Tfactor[0] = 0.;
        }
    }
    else {
        // AreaFract: fraction of grid cell in each snow band
        for (j = 0; j < options.SNOW_BAND; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(&(filenames.params), "AreaFract",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                soil_con[i].AreaFract[j] = (double) dvar[i];
            }
        }
        // elevation: elevation of each snow band
        for (j = 0; j < options.SNOW_BAND; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(&(filenames.params), "elevation",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                soil_con[i].BandElev[j] = (double) dvar[i];
            }
        }
        // Run some checks and corrections for soil
        for (i = 0; i < local_domain.ncells_active; i++) {
            // Make sure area fractions are positive and add to 1
            sum = 0.;
            for (j = 0; j < options.SNOW_BAND; j++) {
                if (soil_con[i].AreaFract[j] < 0) {
                    sprint_location(locstr, &(local_domain.locations[i]));
                    log_err("Negative snow band area fraction (%f) read from "
                            "file\n%s", soil_con[i].AreaFract[j], locstr);
                }
                if (soil_con[i].AreaFract[j] > 0.0) {
                    sum += soil_con[i].AreaFract[j];
                }
            }
            if (!assert_close_double(sum, 1.0, 0., AREA_SUM_ERROR_THRESH)) {
                sprint_location(locstr, &(local_domain.locations[i]));
                if (sum > 0) {
                    log_warn("Sum of the snow band area fractions does not "
                             "equal 1 (%f), dividing each fraction by the "
                             "sum\n%s", sum, locstr);
                    for (j = 0; j < options.SNOW_BAND; j++) {
                        soil_con[i].AreaFract[j] /= sum;
                    }
                }
                else {
                    log_err("Sum of the snow band area fractions is 0\n%s",
                            locstr);
                }
            }
            // check that the mean elevation from the snow bands matches the
            // grid cell mean elevation. If not reset mean
            mean = 0.;
            for (j = 0; j < options.SNOW_BAND; j++) {
                if (soil_con[i].BandElev[j] > 0.0) {
                    mean += soil_con[i].BandElev[j] * soil_con[i].AreaFract[j];
                }
            }
            if (!assert_close_double(soil_con[i].elevation, mean, 0.,
                                     AREA_SUM_ERROR_THRESH)) {
                sprint_location(locstr, &(local_domain.locations[i]));
                log_warn("average band elevation %f not equal to grid_cell "
                         "average elevation %f; setting grid cell elevation "
                         "to average band elevation.\n%s",
                         mean, soil_con[i].elevation, locstr);
                soil_con[i].elevation = (double)mean;
            }
            // Tfactor: calculate the temperature factor
            for (j = 0; j < options.SNOW_BAND; j++) {
                if (soil_con[i].AreaFract[j] > 0) {
                    soil_con[i].Tfactor[j] = (soil_con[i].BandElev[j] -
                                            soil_con[i].elevation) *
                                            param.LAPSE_RATE;
                }
                else {
                    soil_con[i].Tfactor[j] = 0.;
                }
            }
            // Pfactor: calculate Pfactor from the precipitation fraction read
            // from file
            sum = 0.;
            // Calculate Precipitation Fraction
            for (j = 0; j < options.SNOW_BAND; j++) {
                if (soil_con[i].AreaFract[j] > 0) {
                    soil_con[i].Pfactor[j] = (1.0 + param.SNOW_PGRAD * (soil_con[i].BandElev[j] - soil_con[i].elevation)) * 
                                                soil_con[i].AreaFract[j];
                    if (soil_con[i].Pfactor[j] < 0.) {
                        sprint_location(locstr, &(local_domain.locations[i]));
                        log_err("Snow band precipitation fraction (%f) "
                                "must be between 0 and 1.\n%s",
                                soil_con[i].Pfactor[j], locstr);
                    }
                    sum += soil_con[i].Pfactor[j];
                }
                else {
                    soil_con[i].Pfactor[j] = 0.;
                }
            }
            if (!assert_close_double(sum, 1.0, 0., AREA_SUM_ERROR_THRESH)) {
                sprint_location(locstr, &(local_domain.locations[i]));
                log_warn("Sum of the snow band precipitation fractions does "
                         "not equal 1 (%f), dividing each fraction by the "
                         "sum\n%s", sum, locstr);
                for (j = 0; j < options.SNOW_BAND; j++) {
                    soil_con[i].Pfactor[j] /= sum;
                }
            }
            for (j = 0; j < options.SNOW_BAND; j++) {
                if (soil_con[i].AreaFract[j] > 0) {
                    soil_con[i].Pfactor[j] /= soil_con[i].AreaFract[j];
                }
                else {
                    soil_con[i].Pfactor[j] = 0.;
                }
            }
        }
    }

    /******************************************
       Reading the vegetation parameters 
    ******************************************/

    // reading the vegetation parameters is slightly more complicated because
    // VIC allocates memory for veg_con only if the vegetation type exists in
    // the grid cell. The veg_con_map_struct is used to provide some of this
    // mapping

    // number of vegetation types - in vic an extra veg tile is created
    // for above-treeline vegetation in some cases
    // TODO: handle above treeline vegetation tile
    for (i = 0; i < local_domain.ncells_active; i++) {
        nveg = local_domain.locations[i].nveg;
        for (j = 0; j < local_domain.locations[i].nveg; j++) {
            veg_con[i][j].vegetat_type_num = (size_t) nveg;
        }
    }

    // Cv: for each vegetation type, read the cover fraction into the mapping
    // structure. Then assign only the ones with a fraction greater than 0 to
    // the veg_con structure
    for (j = 0; j < options.MAX_HRU; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(&(filenames.params), "Cv",
                                d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (j < local_domain.locations[i].nveg) {
                veg_con[i][j].Cv = (double) dvar[i];
            }
        }
    }
    for (j = 0; j < options.MAX_HRU; j++) {
        d3start[0] = j;
        get_scatter_nc_field_int(&(filenames.params), "BandIndex",
                                d3start, d3count, ivar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (j < local_domain.locations[i].nveg) {
                veg_con[i][j].BandIndex = (int) ivar[i];
            }
        }
    }
    for (j = 0; j < options.MAX_HRU; j++) {
        d3start[0] = j;
        get_scatter_nc_field_int(&(filenames.params), "veg_class",
                                d3start, d3count, ivar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (j < local_domain.locations[i].nveg) {
                veg_con[i][j].veg_class = (int) ivar[i];
                if (veg_con[i][j].veg_class == options.GLACIER_ID) {
                    veg_con[i][j].IS_GLAC = true;
                }
                else {
                    veg_con[i][j].IS_GLAC = false;
                }
            }
        } 
    }

    // do the mapping
    for (i = 0; i < local_domain.ncells_active; i++) {
        for (j = 0; j < local_domain.locations[i].nveg; j++) {
            k = veg_con[i][j].veg_class;
            for (m = 0; m < MONTHS_PER_YEAR; m++) {
                if (options.FCAN_SRC == FROM_DEFAULT ||
                    options.FCAN_SRC == FROM_VEGLIB ||
                    options.FCAN_SRC == FROM_VEGPARAM) {
                    veg_con[i][j].fcanopy[m] = veg_lib[k].fcanopy[m];
                }
                if (options.LAI_SRC == FROM_VEGLIB ||
                    options.LAI_SRC == FROM_VEGPARAM) {
                    veg_con[i][j].LAI[m] = veg_lib[k].LAI[m];
                }
                if (options.LAI_SRC == FROM_VEGLIB ||
                    options.LAI_SRC == FROM_VEGPARAM) { 
                    veg_con[i][j].SAI[m] = veg_lib[k].SAI[m];
                }           
            }
        }
        // check the number of nonzero veg tiles
        if (j > local_domain.locations[i].nveg + 1) {
            sprint_location(locstr, &(local_domain.locations[i]));
            log_err("Number of veg tiles with nonzero area (%zu) > nveg + 1 "
                    "(%zu).\n%s", j, local_domain.locations[i].nveg,
                    locstr);
        }
        else if (j < local_domain.locations[i].nveg) {
            sprint_location(locstr, &(local_domain.locations[i]));
            log_err("Number of veg tiles with nonzero area (%zu) < nveg "
                    "(%zu).\n%s", j, local_domain.locations[i].nveg,
                    locstr);
        }
    }

    // calculate root fractions
    for (i = 0; i < local_domain.ncells_active; i++) {
        calc_root_fractions(veg_con[i], &(soil_con[i]));
    }

    // Run some checks and corrections for vegetation
    for (i = 0; i < local_domain.ncells_active; i++) {
        // Only run to options.NVEGTYPES - 1, assuming bare soil
        // is the last type
        for (j = 0; j < local_domain.locations[i].nveg; j++) {
            vidx = veg_con[i][j].veg_class;
            if (vidx != NODATA_VEG) {
                // check that the vegetation type is defined in the vegetation
                // library
                found = false;
                for (k = 0; k < options.NVEGTYPES; k++) {
                    if (veg_con[i][j].veg_class == veg_lib[k].veg_class) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    sprint_location(locstr, &(local_domain.locations[i]));
                    log_err("The vegetation class id %d in vegetation tile %zu "
                            "from cell %zd is not defined in the vegetation "
                            "library\n%s", veg_con[i][j].veg_class, j, i,
                            locstr);
                }
                Cv_sum[i] += veg_con[i][j].Cv;
            }
        }
        // handle the bare soil portion of the tile
        vidx = veg_con[i][j].Cv;
        if (vidx != NODATA_VEG) {
            Cv_sum[i] += veg_con[i][j].Cv;
        }

        // TODO: handle bare soil adjustment for compute treeline option

        // If the sum of the tile fractions is not within a tolerance,
        // readjust Cvs to sum to 1.0
        if (!assert_close_double(Cv_sum[i], 1., 0., AREA_SUM_ERROR_THRESH)) {
            sprint_location(locstr, &(local_domain.locations[i]));
            log_warn("Sum of veg tile area fractions !=  1.0 (%.16f) at grid "
                     "cell %zd. Adjusting fractions ...\n%s", Cv_sum[i], i,
                     locstr);
            for (j = 0; j < local_domain.locations[i].nveg; j++) {
                vidx = veg_con[i][j].Cv;
                if (vidx != NODATA_VEG) {
                    veg_con[i][j].Cv /= Cv_sum[i];
                }
            }
        }
    }

    // initialize state variables with default values
    for (i = 0; i < local_domain.ncells_active; i++) {
        nveg = veg_con[i][0].vegetat_type_num;
        initialize_snow(all_vars[i].snow, nveg);
        initialize_soil(all_vars[i].cell, nveg);
        initialize_veg(all_vars[i].veg_var, nveg);
        initialize_energy(all_vars[i].energy, nveg);
    }

    // set state metadata structure
    set_state_meta_data_info();

    // close parameter file
    if (mpi_rank == VIC_MPI_ROOT) {
        status = nc_close(filenames.params.nc_id);
        check_nc_status(status, "Error closing %s",
                        filenames.params.nc_filename);
    }

    // cleanup
    free(dvar);
    free(ivar);
    free(Cv_sum);
    free(temp_array);
}
