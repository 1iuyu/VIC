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
    extern optical_struct      optical;
    extern int                 mpi_rank;

    bool                       found;
    char                       locstr[MAXSTRING];
    double                     mean;
    double                     sum;
    double                    *array = NULL;
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
    array = (double*) malloc(total_size * sizeof(*array));
    check_alloc_status(array, "Memory allocation error.");
    
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

    // avg_temp: average soil temperature (K)
    get_scatter_nc_field_double(&(filenames.params), "avg_temp",
                                d2start, d2count, dvar);
    for (i = 0; i < local_domain.ncells_active; i++) {
        soil_con[i].avg_temp = (double) dvar[i];
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

    // topo_std: topographic std of the cell []
    get_scatter_nc_field_double(&(filenames.params), "topo_std",
                                d2start, d2count, dvar);
    for (i = 0; i < local_domain.ncells_active; i++) {
        soil_con[i].topo_std = (double) dvar[i];
    }

    // init_zwt: initial water table depth [m]
    get_scatter_nc_field_double(&(filenames.params), "init_zwt",
                                d2start, d2count, dvar);
    for (i = 0; i < local_domain.ncells_active; i++) {
        soil_con[i].init_zwt = (double) dvar[i];
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
        double current_thick = 0.5;
        while (k < MAX_SOILS && accum_depth < z_bedrock) {
            /* 检查加上当前层是否会超过总深度 */
            if (accum_depth + current_thick < z_bedrock) {
                /* 正常分配 */
                soil_con[i].dz_soil[k] = current_thick;
                accum_depth += current_thick;
                current_thick *= ratio;
            } else {
                /* 最后一层：设为基岩层 */
                soil_con[i].dz_soil[k] = current_thick;
                accum_depth = z_bedrock;
                break;
            }
            k++;
            /* 防止厚度增长过大 */
            if (current_thick > 5.0) {
                current_thick = 5.0;
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
                                    soil_con[i].dz_soil[k] / 2.0;
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
            array[lidx] = (double) dvar[i];
        }
    }
    for (i = 0; i < local_domain.ncells_active; i++) {
        set_node_parameters(soil_con[i].Nbedrock,
                            soil_con[i].depth,
                            soil_con[i].Zsum_soil,
                            &soil_con[i].clay_node[0],
                            &array[i * Nlayer]);
    }

    // sand: sand content for each soil layer
    for (j = 0; j < Nlayer; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(&(filenames.params), "sand",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            lidx = i * Nlayer + j;
            array[lidx] = (double) dvar[i];
        }
    }
    for (i = 0; i < local_domain.ncells_active; i++) {
        set_node_parameters(soil_con[i].Nbedrock,
                            soil_con[i].depth,
                            soil_con[i].Zsum_soil,
                            &soil_con[i].sand_node[0],
                            &array[i * Nlayer]);
    }

    // silt: silt content for each soil layer
    for (j = 0; j < Nlayer; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(&(filenames.params), "silt",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            lidx = i * Nlayer + j;
            array[lidx] = (double) dvar[i];
        }
    }
    for (i = 0; i < local_domain.ncells_active; i++) {
        set_node_parameters(soil_con[i].Nbedrock,
                            soil_con[i].depth,
                            soil_con[i].Zsum_soil,
                            &soil_con[i].silt_node[0],
                            &array[i * Nlayer]);
    }

    // gravel: gravel content for each soil layer
    for (j = 0; j < Nlayer; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(&(filenames.params), "gravel",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            lidx = i * Nlayer + j;
            array[lidx] = (double) dvar[i];
        }
    }
    for (i = 0; i < local_domain.ncells_active; i++) {
        set_node_parameters(soil_con[i].Nbedrock,
                            soil_con[i].depth,
                            soil_con[i].Zsum_soil,
                            &soil_con[i].gravel_node[0],
                            &array[i * Nlayer]);
    }

    // Wsat: saturated point for each layer
    for (j = 0; j < Nlayer; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(&(filenames.params), "organic",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            lidx = i * Nlayer + j;
            array[lidx] = (double) dvar[i];
        }
    }
    for (i = 0; i < local_domain.ncells_active; i++) {
        set_node_parameters(soil_con[i].Nbedrock,
                            soil_con[i].depth,
                            soil_con[i].Zsum_soil,
                            &soil_con[i].organic_node[0],
                            &array[i * Nlayer]);
    }

    for (j = 0; j < Nlayer; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(&(filenames.params), "bulk_density",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            lidx = i * Nlayer + j;
            array[lidx] = (double) dvar[i];
        }
    }
    for (i = 0; i < local_domain.ncells_active; i++) {
        set_node_parameters(soil_con[i].Nbedrock,
                            soil_con[i].depth,
                            soil_con[i].Zsum_soil,
                            &soil_con[i].bulk_dens_node[0],
                            &array[i * Nlayer]);
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
                array[lidx] = (double) dvar[i];
            }
        }
        for (i = 0; i < local_domain.ncells_active; i++) {
            set_node_parameters(soil_con[i].Nbedrock,
                                soil_con[i].depth,
                                soil_con[i].Zsum_soil,
                                &soil_con[i].alpha_node[0],
                                &array[i * Nlayer]);
        }
        // expt: layer-specific exponent n in van Genuchten eqn
        for (j = 0; j < Nlayer; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(&(filenames.params), "expt",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                lidx = i * Nlayer + j;
                array[lidx] = (double) dvar[i];
            }
        }
        for (i = 0; i < local_domain.ncells_active; i++) {
            set_node_parameters(soil_con[i].Nbedrock,
                                soil_con[i].depth,
                                soil_con[i].Zsum_soil,
                                &soil_con[i].expt_node[0],
                                &array[i * Nlayer]);
        }
        // Ksat: saturated hydraulic conductivity [m/s]
        for (j = 0; j < Nlayer; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(&(filenames.params), "Ksat",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                lidx = i * Nlayer + j;
                array[lidx] = (double) dvar[i];
            }
        }
        for (i = 0; i < local_domain.ncells_active; i++) {
            set_node_parameters(soil_con[i].Nbedrock,
                                soil_con[i].depth,
                                soil_con[i].Zsum_soil,
                                &soil_con[i].Ksat_node[0],
                                &array[i * Nlayer]);
        }
        // Wpwp: soil moisture content at permanent wilting point [m3/m3]
        for (j = 0; j < Nlayer; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(&(filenames.params), "Wpwp",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                lidx = i * Nlayer + j;
                array[lidx] = (double) dvar[i];
            }
        }
        for (i = 0; i < local_domain.ncells_active; i++) {
            set_node_parameters(soil_con[i].Nbedrock,
                                soil_con[i].depth,
                                soil_con[i].Zsum_soil,
                                &soil_con[i].Wpwp_node[0],
                                &array[i * Nlayer]);
        }
        // Wsat: soil moisture content at saturation [m3/m3]
        for (j = 0; j < Nlayer; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(&(filenames.params), "Wsat",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                lidx = i * Nlayer + j;
                array[lidx] = (double) dvar[i];
            }
        }
        for (i = 0; i < local_domain.ncells_active; i++) {
            set_node_parameters(soil_con[i].Nbedrock,
                                soil_con[i].depth,
                                soil_con[i].Zsum_soil,
                                &soil_con[i].Wsat_node[0],
                                &array[i * Nlayer]);
        }
        // lpar: unsaturated hydraulic conductivity exponent.
        for (j = 0; j < Nlayer; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(&(filenames.params), "lpar",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                lidx = i * Nlayer + j;
                array[lidx] = (double) dvar[i];
            }
        }
        for (i = 0; i < local_domain.ncells_active; i++) {
            set_node_parameters(soil_con[i].Nbedrock,
                                soil_con[i].depth,
                                soil_con[i].Zsum_soil,
                                &soil_con[i].lpar_node[0],
                                &array[i * Nlayer]);
        }
        for (i = 0; i < local_domain.ncells_active; i++) {
            for (j = 0; j < soil_con[i].Nbedrock - 1; j++) {
                soil_con[i].mpar_node[j] = 1.0 - 1.0 / soil_con[i].expt_node[i];
                soil_con[i].bulk_dens_node[j] = (soil_con[i].bulk_dens_node[j] * 
                    (1.0 - soil_con[i].gravel_node[j]) + soil_con[i].gravel_node[j] * 2650);
            }
        }     

    }
    else {
        // 土壤参数从PedoTransfer函数中计算得到
        PedoTransfer(soil_con); // not used in current version.       
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
                double value  = (double) dvar[i];
                if (!isnan(value) && !isinf(value) && value > 0) {
                    soil_con[i].AreaFract[j] = value;
                }
            }
        }
        // elevation: elevation of each snow band
        for (j = 0; j < options.SNOW_BAND; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(&(filenames.params), "elevation",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                double value  = (double) dvar[i];
                if (!isnan(value) && !isinf(value) && value > 0) {
                    soil_con[i].BandElev[j] = (double) dvar[i];
                }
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
                    soil_con[i].Tfactor[j] = 0.0;
                }
            }
            // Pfactor: calculate Pfactor from the precipitation fraction read
            // from file
            sum = 0.;
            // Calculate Precipitation Fraction
            for (j = 0; j < options.SNOW_BAND; j++) {
                if (soil_con[i].AreaFract[j] > 0) {
                    soil_con[i].Pfactor[j] = (1.0 + param.SNOW_PGRAD * (soil_con[i].BandElev[j] 
                        - soil_con[i].elevation)) * soil_con[i].AreaFract[j];
                    if (soil_con[i].Pfactor[j] < 0.) {
                        sprint_location(locstr, &(local_domain.locations[i]));
                        log_err("Snow band precipitation fraction (%f) "
                                "must be between 0 and 1.\n%s",
                                soil_con[i].Pfactor[j], locstr);
                    }
                    sum += soil_con[i].Pfactor[j];
                }
                else {
                    soil_con[i].Pfactor[j] = 0.0;
                }
            }
            if (!assert_close_double(sum, 1.0, 0.0, AREA_SUM_ERROR_THRESH)) {
                sprint_location(locstr, &(local_domain.locations[i]));
                log_warn("Sum of the snow band precipitation fractions does "
                         "not equal 1 (%f), dividing each fraction by the "
                         "sum\n%s", sum, locstr);
                for (j = 0; j < options.SNOW_BAND; j++) {
                    soil_con[i].Pfactor[j] /= sum;
                }
            }
            for (j = 0; j < options.SNOW_BAND; j++) {
                if (soil_con[i].AreaFract[j] > 0.0) {
                    soil_con[i].Pfactor[j] /= soil_con[i].AreaFract[j];
                }
                else {
                    soil_con[i].Pfactor[j] = 0.0;
                }
            }
        }
    }

    /******************************************
      Reading library of vegetation parameters
    ******************************************/
    // Initialize veg_lib parameters
    size_t n1dims = 1, n2dims = 2;
    size_t start1[1] = {0};
    size_t count1[1] = {(size_t) options.NVEGTYPES};
    size_t start2[2] = {0, 0};
    size_t count2[2] = {
            (size_t) options.NVEGTYPES,
            (size_t) MONTHS_PER_YEAR};
    double local_d1[options.NVEGTYPES];
    int local_int[options.NVEGTYPES];
    double d2_size = options.NVEGTYPES * MONTHS_PER_YEAR;
    double *local_d2 = malloc(d2_size * sizeof(double));
    check_alloc_status(local_d2, "Memory allocation error.");

    // 读取 Canopy_Upper
    get_scatter_nc_table_int(&(filenames.params), "IGBP_class", 
                             n1dims, start1, count1, local_int);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].veg_class = local_int[i];
    }

    // 读取 Canopy_Upper
    get_scatter_nc_table_double(&(filenames.params), "Canopy_Upper", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].Canopy_Upper = local_d1[i];
    }

    // 读取 Canopy_Lower
    get_scatter_nc_table_double(&(filenames.params), "Canopy_Lower", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].Canopy_Lower = local_d1[i];
    }

    // 读取 Canopy_Radius
    get_scatter_nc_table_double(&(filenames.params), "Canopy_Radius", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].Canopy_Radius = local_d1[i];
    }

    // 读取 COI
    get_scatter_nc_table_double(&(filenames.params), "COI", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].COI = local_d1[i];
    }

    // 读取 c_biomass
    get_scatter_nc_table_double(&(filenames.params), "c_biomass", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].c_biomass = local_d1[i];
    }

    // 读取 d_leaf
    get_scatter_nc_table_double(&(filenames.params), "d_leaf", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].d_leaf = local_d1[i];
    }

    // 读取 root_a
    get_scatter_nc_table_double(&(filenames.params), "root_a", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].root_a = local_d1[i];
    }

    // 读取 root_b
    get_scatter_nc_table_double(&(filenames.params), "root_b", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].root_b = local_d1[i];
    }

    // 读取 root_d
    get_scatter_nc_table_double(&(filenames.params), "root_d", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].root_d = local_d1[i];
    }

    // 读取 liq_bioms
    get_scatter_nc_table_double(&(filenames.params), "liq_bioms", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].liq_bioms = local_d1[i];
    }

    // 读取 slatop
    get_scatter_nc_table_double(&(filenames.params), "slatop", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].slatop = local_d1[i];
    }

    // 读取 stem_num
    get_scatter_nc_table_double(&(filenames.params), "stem_num", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].stem_num = local_d1[i];
    }

    // 读取 Z0sub_LAImax
    get_scatter_nc_table_double(&(filenames.params), "Z0sub_LAImax", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].Z0sub_LAImax = local_d1[i];
    }

    // 读取 Z0sub_Cs
    get_scatter_nc_table_double(&(filenames.params), "Z0sub_Cs", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].Z0sub_Cs = local_d1[i];
    }

    // 读取 Z0sub_Cr
    get_scatter_nc_table_double(&(filenames.params), "Z0sub_Cr", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].Z0sub_Cr = local_d1[i];
    }

    // 读取 Z0sub_c
    get_scatter_nc_table_double(&(filenames.params), "Z0sub_c", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].Z0sub_c = local_d1[i];
    }

    // 读取 Z0sub_cw
    get_scatter_nc_table_double(&(filenames.params), "Z0sub_cw", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].Z0sub_cw = local_d1[i];
    }

    // 读取 smpsc
    get_scatter_nc_table_double(&(filenames.params), "smpsc", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].smpsc = local_d1[i];
    }

    // 读取 smpso
    get_scatter_nc_table_double(&(filenames.params), "smpso", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].smpso = local_d1[i];
    }

    // 读取 trunk_dia
    get_scatter_nc_table_double(&(filenames.params), "trunk_dia", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].trunk_dia = local_d1[i];
    }

    // 读取 froot_leaf
    get_scatter_nc_table_double(&(filenames.params), "froot_leaf", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].froot_leaf = local_d1[i];
    }

    // 读取 theta_cj
    get_scatter_nc_table_double(&(filenames.params), "theta_cj", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].theta_cj = local_d1[i];
    }

    // 读取 kcano_max
    get_scatter_nc_table_double(&(filenames.params), "kcano_max", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].kcano_max = local_d1[i];
    }

    // 读取 kroot_max
    get_scatter_nc_table_double(&(filenames.params), "kroot_max", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].kroot_max = local_d1[i];
    }

    // 读取 matric50
    get_scatter_nc_table_double(&(filenames.params), "matric50", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].matric50 = local_d1[i];
    }

    // 读取 leaf_CN
    get_scatter_nc_table_double(&(filenames.params), "leaf_CN", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].leaf_CN = local_d1[i];
    }

    // 读取 SLA_top
    get_scatter_nc_table_double(&(filenames.params), "SLA_top", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].SLA_top = local_d1[i];
    }

    // 读取 fN_rub
    get_scatter_nc_table_double(&(filenames.params), "fN_rub", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].fN_rub = local_d1[i];
    }

    // 读取 medlynslope
    get_scatter_nc_table_double(&(filenames.params), "medlynslope", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].medlynslope = local_d1[i];
    }

    // 读取 medlynint
    get_scatter_nc_table_double(&(filenames.params), "medlynint", 
                                n1dims, start1, count1, local_d1);
    for (i = 0; i < options.NVEGTYPES; i++) {
        veg_lib[i].medlynint = local_d1[i];
    }

    // 读取 landtype
    get_scatter_nc_table_int(&(filenames.params), "landtype",
                             n1dims, start1, count1, local_int);

    for (i = 0; i < options.NVEGTYPES; i++) {
        if (local_int[i] == 0) {
            veg_lib[i].landtype = LAND_SOIL;
        }
        else if (local_int[i] == 1) {
            veg_lib[i].landtype = LAND_GLAC;
        }
        else if (local_int[i] == 2) {
            veg_lib[i].landtype = LAND_WET;
        }
        else if (local_int[i] == 3) {
            veg_lib[i].landtype = LAND_URBAN;
        }
        else {
            log_warn("Warning: Unknown landtype %d for veg %d, setting to default LAND_SOIL", 
                    local_int[i], i);
            veg_lib[i].landtype = LAND_SOIL;  // 默认值
        }
    }

    // 读取 Ctype
    get_scatter_nc_table_int(&(filenames.params), "Ctype",
                             n1dims, start1, count1, local_int);

    for (i = 0; i < options.NVEGTYPES; i++) {
        if (local_int[i] == 0) {
            veg_lib[i].Ctype = PHOTO_C3;
        }
        else if (local_int[i] == 1) {
            veg_lib[i].Ctype = PHOTO_C4;
        }
        else {
            log_warn("Warning: Unknown Ctype %d for veg %d, setting to default C3", 
                    local_int[i], i);
            veg_lib[i].Ctype = PHOTO_C3;  // 默认值
        }
    }

    // 读取 LAI
    get_scatter_nc_table_double(&(filenames.params), "LAI", 
                                n2dims, start2, count2, local_d2);
    for (i = 0; i < options.NVEGTYPES; i++) {
        for (m = 0; m < MONTHS_PER_YEAR; m++) {
            veg_lib[i].LAI[m] = local_d2[i * MONTHS_PER_YEAR + m];
        }
    }

    // 读取 SAI
    get_scatter_nc_table_double(&(filenames.params), "SAI", 
                                n2dims, start2, count2, local_d2);
    for (i = 0; i < options.NVEGTYPES; i++) {
        for (m = 0; m < MONTHS_PER_YEAR; m++) {
            veg_lib[i].SAI[m] = local_d2[i * MONTHS_PER_YEAR + m];
        }
    }

    // 读取 fcanopy
    if (options.FCAN_SRC != FROM_DEFAULT) {
        get_scatter_nc_table_double(&(filenames.params), "fcanopy", 
                                    n2dims, start2, count2, local_d2);
        for (i = 0; i < options.NVEGTYPES; i++) {
            for (m = 0; m < MONTHS_PER_YEAR; m++) {
                veg_lib[i].fcanopy[m] = local_d2[i * MONTHS_PER_YEAR + m];
                if (veg_lib[i].fcanopy[m] < 0 ||
                    veg_lib[i].fcanopy[m] > 1) {
                    log_err(
                        "Veg cover fraction must be between 0 and 1 " "(%f)",
                        veg_lib[i].fcanopy[m]);
                }
            }
        }
    }
    // 设置植被反射和透射率count2
    count2[1] = (size_t) MAX_SWBANDS;
    d2_size = options.NVEGTYPES * MAX_SWBANDS;
    double *local_sw = malloc(d2_size * sizeof(double));
    check_alloc_status(local_sw, "Memory allocation error.");

    // 读取 reflleaf
    get_scatter_nc_table_double(&(filenames.params), "reflleaf", 
                                n2dims, start2, count2, local_sw);
    for (i = 0; i < options.NVEGTYPES; i++) {
        for (m = 0; m < MAX_SWBANDS; m++) {
            veg_lib[i].reflleaf[m] = local_sw[i * MAX_SWBANDS + m];
        }
    }

    // 读取 reflstem
    get_scatter_nc_table_double(&(filenames.params), "reflstem", 
                                n2dims, start2, count2, local_sw);
    for (i = 0; i < options.NVEGTYPES; i++) {
        for (m = 0; m < MAX_SWBANDS; m++) {
            veg_lib[i].reflstem[m] = local_sw[i * MAX_SWBANDS + m];
        }
    }

    // 读取 transleaf
    get_scatter_nc_table_double(&(filenames.params), "transleaf", 
                                n2dims, start2, count2, local_sw);
    for (i = 0; i < options.NVEGTYPES; i++) {
        for (m = 0; m < MAX_SWBANDS; m++) {
            veg_lib[i].transleaf[m] = local_sw[i * MAX_SWBANDS + m];
        }
    }

    // 读取 transstem
    get_scatter_nc_table_double(&(filenames.params), "transstem", 
                                n2dims, start2, count2, local_sw);
    for (i = 0; i < options.NVEGTYPES; i++) {
        for (m = 0; m < MAX_SWBANDS; m++) {
            veg_lib[i].transstem[m] = local_sw[i * MAX_SWBANDS + m];
        }
    }

    free(local_d2);
    free(local_sw);

    /******************************************
       Reading the vegetation parameters 
    ******************************************/

    // number of vegetation types - in vic an extra veg tile is created
    // for above-treeline vegetation in some cases
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

    /******************************************
       Reading the SNICAR parameters 
    ******************************************/
    if (options.SNOW_ALBEDO == SNICAR) {

        // Initialize optical parameters
        memset(&optical, 0, sizeof(optical_struct));

        size_t n1dims = 1, n2dims = 2, n3dims = 3;
        size_t start1[1] = {0};
        size_t count1[1] = {(size_t) SNICAR_BANDS};
        size_t start[2] = {0, 0};
        size_t count[2] = {
                (size_t) SNICAR_BANDS,
                (size_t) SNICAR_RADII};
        size_t start3[3] = {0, 0, 0};
        size_t count3[3] = {
                (size_t) LOOKUP_TEMP,
                (size_t) LOOKUP_DTDZ,
                (size_t) LOOKUP_DENS};
        double d1_size = SNICAR_BANDS;
        double local_d1[SNICAR_BANDS];
        double d2_size = SNICAR_BANDS * SNICAR_RADII;
        double *local_d2 = malloc(d2_size * sizeof(double));
        check_alloc_status(local_d2, "Memory allocation error.");
        double d3_size = LOOKUP_TEMP * LOOKUP_DTDZ * LOOKUP_DENS;
        double *local_d3 = malloc(d3_size * sizeof(double));
        check_alloc_status(local_d3, "Memory allocation error.");
        
        // Mie single scatter albedos for direct-beam ice
        get_scatter_nc_table_double(&(filenames.params), "ss_alb_dir", n2dims,
                                    start, count, local_d2);
        memcpy(optical.ss_alb_dir, local_d2, d2_size * sizeof(double));

        // Mie single scatter albedos for diffuse ice
        get_scatter_nc_table_double(&(filenames.params), "ss_alb_dfs", n2dims,
                                    start, count, local_d2);
        memcpy(optical.ss_alb_dfs, local_d2, d2_size * sizeof(double));

        // asymmetry parameter of direct-beam ice  
        get_scatter_nc_table_double(&(filenames.params), "asym_snow_dir", n2dims,
                                    start, count, local_d2);
        memcpy(optical.asym_snow_dir, local_d2, d2_size * sizeof(double));

        // asymmetry parameter of diffuse ice
        get_scatter_nc_table_double(&(filenames.params), "asym_snow_dfs", n2dims,
                                    start, count, local_d2);
        memcpy(optical.asym_snow_dfs, local_d2, d2_size * sizeof(double));

        // mass extinction coefficient for direct-beam ice [m2/kg]
        get_scatter_nc_table_double(&(filenames.params), "mass_ext_dir", n2dims,
                                    start, count, local_d2);
        memcpy(optical.mass_ext_dir, local_d2, d2_size * sizeof(double));

        // mass extinction coefficient for diffuse ice [m2/kg]
        get_scatter_nc_table_double(&(filenames.params), "mass_ext_dfs", n2dims,
                                    start, count, local_d2);
        memcpy(optical.mass_ext_dfs, local_d2, d2_size * sizeof(double));

        // snow aging parameter retrieved from lookup table [hour]
        get_scatter_nc_table_double(&(filenames.params), "tau_table", n3dims,
                              start3, count3, local_d3);
        memcpy(optical.tau_table, local_d3, d3_size * sizeof(double));

        // snow aging parameter retrieved from lookup table [unitless]
        get_scatter_nc_table_double(&(filenames.params), "kappa_table", n3dims,
                              start3, count3, local_d3);
        memcpy(optical.kappa_table, local_d3, d3_size * sizeof(double));

        // snow aging parameter retrieved from lookup table [um hr-1]
        get_scatter_nc_table_double(&(filenames.params), "drdt_table", n3dims,
                              start3, count3, local_d3);
        memcpy(optical.drdt_table, local_d3, d3_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ss_alb_bcphil_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ss_alb_bcphil_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ss_alb_bcphil_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ss_alb_bcphil_dfs, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ss_alb_bcphob_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ss_alb_bcphob_dfs, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ss_alb_bcphob_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ss_alb_bcphob_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ss_alb_dust1_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ss_alb_dust1_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ss_alb_dust1_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ss_alb_dust1_dfs, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ss_alb_dust2_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ss_alb_dust2_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ss_alb_dust2_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ss_alb_dust2_dfs, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ss_alb_dust3_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ss_alb_dust3_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ss_alb_dust3_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ss_alb_dust3_dfs, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ss_alb_dust4_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ss_alb_dust4_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ss_alb_dust4_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ss_alb_dust4_dfs, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ss_alb_dust5_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ss_alb_dust5_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ss_alb_dust5_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ss_alb_dust5_dfs, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "asm_prm_bcphob_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.asm_prm_bcphob_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "asm_prm_bcphob_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.asm_prm_bcphob_dfs, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ss_alb_ocphil_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ss_alb_ocphil_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ss_alb_ocphil_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ss_alb_ocphil_dfs, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "asm_prm_ocphil_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.asm_prm_ocphil_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "asm_prm_ocphil_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.asm_prm_ocphil_dfs, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ext_cff_mss_ocphil_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ext_cff_mss_ocphil_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ext_cff_mss_ocphil_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ext_cff_mss_ocphil_dfs, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ss_alb_ocphob_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ss_alb_ocphob_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ss_alb_ocphob_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ss_alb_ocphob_dfs, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "asm_prm_ocphob_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.asm_prm_ocphob_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "asm_prm_ocphob_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.asm_prm_ocphob_dfs, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ext_cff_mss_ocphob_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ext_cff_mss_ocphob_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ext_cff_mss_ocphob_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ext_cff_mss_ocphob_dfs, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "asm_prm_bcphil_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.asm_prm_bcphil_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "asm_prm_bcphil_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.asm_prm_bcphil_dfs, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ext_cff_mss_bcphob_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ext_cff_mss_bcphob_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ext_cff_mss_bcphob_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ext_cff_mss_bcphob_dfs, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ext_cff_mss_bcphil_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ext_cff_mss_bcphil_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ext_cff_mss_bcphil_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ext_cff_mss_bcphil_dfs, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "asm_prm_dust1_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.asm_prm_dust1_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "asm_prm_dust1_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.asm_prm_dust1_dfs, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "asm_prm_dust2_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.asm_prm_dust2_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "asm_prm_dust2_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.asm_prm_dust2_dfs, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "asm_prm_dust3_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.asm_prm_dust3_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "asm_prm_dust3_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.asm_prm_dust3_dfs, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "asm_prm_dust4_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.asm_prm_dust4_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "asm_prm_dust4_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.asm_prm_dust4_dfs, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "asm_prm_dust5_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.asm_prm_dust5_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "asm_prm_dust5_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.asm_prm_dust5_dfs, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ext_cff_mss_dust1_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ext_cff_mss_dust1_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ext_cff_mss_dust1_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ext_cff_mss_dust1_dfs, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ext_cff_mss_dust2_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ext_cff_mss_dust2_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ext_cff_mss_dust2_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ext_cff_mss_dust2_dfs, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ext_cff_mss_dust3_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ext_cff_mss_dust3_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ext_cff_mss_dust3_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ext_cff_mss_dust3_dfs, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ext_cff_mss_dust4_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ext_cff_mss_dust4_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ext_cff_mss_dust4_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ext_cff_mss_dust4_dfs, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ext_cff_mss_dust5_dir", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ext_cff_mss_dust5_dir, local_d1, d1_size * sizeof(double));

        get_scatter_nc_table_double(&(filenames.params), "ext_cff_mss_dust5_dfs", 
                                    n1dims, start1, count1, local_d1);
        memcpy(optical.ext_cff_mss_dust5_dfs, local_d1, d1_size * sizeof(double));

        free(local_d2);
        free(local_d3);
    }
    
    // initialize state variables with default values
    for (i = 0; i < local_domain.ncells_active; i++) {
        nveg = veg_con[i][0].vegetat_type_num;
        initialize_snow(all_vars[i].snow, nveg);
        initialize_soil(all_vars[i].cell, nveg);
        initialize_veg(all_vars[i].veg_var, nveg);
        initialize_energy(all_vars[i].energy, nveg);
    }

    size_t veg_class = 0;
    // Initialize landunit types based on vegetation class
    for (i = 0; i < local_domain.ncells_active; i++) {
        for (j = 0; j <= local_domain.locations[i].nveg; j++) {
            if (veg_con[i][j].Cv > 0) {
                veg_class = veg_con[i][j].veg_class;
                if (veg_lib[veg_class].landtype == 0) {
                    all_vars[i].cell[j].IS_VEG = true;
                }
                else if (veg_lib[veg_class].landtype == 1) {
                    all_vars[i].cell[j].IS_GLAC = true;
                }
                else if (veg_lib[veg_class].landtype == 2) {
                    all_vars[i].cell[j].IS_WET = true;
                }
                else if (veg_lib[veg_class].landtype == 3) {
                    all_vars[i].cell[j].IS_URBAN = true;
                }
                else {
                    log_err("Unknown Landtype option");
                }
            }
        }
    }

    // Initialize layer roots fraction
    for (i = 0; i < local_domain.ncells_active; i++) {
        for (j = 0; j <= local_domain.locations[i].nveg; j++) {
            if (veg_con[i][j].Cv > 0) {
                veg_class = veg_con[i][j].veg_class;
                calc_root_fractions(veg_class, 
                                    &all_vars[i].cell[j],
                                    &all_vars[i].veg_var[j],
                                    soil_con, veg_lib);
            }
        }
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
    free(array);
}