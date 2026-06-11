/******************************************************************************
 * @section DESCRIPTION
 *
 * Initilization library.
 *****************************************************************************/

#include "vic_driver_shared_image.h"

/******************************************************************************
 * @brief    Initialize soil con sructure.
 *****************************************************************************/
void
initialize_soil_con(soil_con_struct *soil_con)
{
    extern option_struct options;
    size_t               i;
    
    soil_con->gridcel = MISSING_USI;

    soil_con->AlbedoPar = 0.0;
    soil_con->avg_temp = 0.0;
    soil_con->elevation = 0.0;
    soil_con->lat = 0.0;
    soil_con->lng = 0.0;
    soil_con->time_zone_lng = 0.0;
    soil_con->b_infilt = 0.0;
    soil_con->b_dynamic = 0.0;
    soil_con->capil_drive = 0.0;
    soil_con->cell_area = 0.0;
    soil_con->slope = 0.0;
    soil_con->init_zwt = 0.0;

    for (i = 0; i < MAX_LAYERS; i++) {
        soil_con->bulk_dens_min[i] = 0.0;
        soil_con->bulk_dens_org[i] = 0.0;
        soil_con->depth[i] = 0.0;
    }

    for (i = 0; i < MAX_SOILS; i++) {
        soil_con->expt_node[i] = 0.0;
        soil_con->mpar_node[i] = 0.0;
        soil_con->bulk_dens_node[i] = 0.0;
        soil_con->bubble_node[i] = 0.0;
        soil_con->clay_node[i] = 0.0;
        soil_con->sand_node[i] = 0.0;
        soil_con->silt_node[i] = 0.0;
        soil_con->gravel_node[i] = 0.0;
        soil_con->dz_soil[i] = 0.0;
        soil_con->Zsum_soil[i] = 0.0;
        soil_con->zc_soil[i] = 0.0;
        soil_con->organic_node[i] = 0.0;
        soil_con->soil_dens_node[i] = 0.0;
        soil_con->Ksat_node[i] = 0.0;
        soil_con->Wpwp_node[i] = 0.0;
        soil_con->Wsat_node[i] = 0.0;
        soil_con->soil_dens_min[i] = 2650.0;
        soil_con->soil_dens_org[i] = 1300.0;
    }

    for (i = 0; i < options.SNOW_BAND; i++) {
        soil_con->BandElev[i] = 0.;
        soil_con->AreaFract[i] = 1.;
        soil_con->Pfactor[i] = 0.;
        soil_con->Tfactor[i] = 0.;
    }
}

/******************************************************************************
 * @brief    Initialize veg con sructure.
 *****************************************************************************/
void
initialize_veg_con(veg_con_struct *veg_con)
{
    size_t   i;

    veg_con->Cv = 0.;
    veg_con->veg_class = NODATA_VEG; // -1 to force a crash if inappropriate
    veg_con->vegetat_type_num = 0;
    veg_con->BandIndex = 0;
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        veg_con->LAI[i] = 0.0;
        veg_con->SAI[i] = 0.0;
        veg_con->fcanopy[i] = 0.0;
    }
}

/******************************************************************************
 * @brief    Initialize domain info stucture
 *****************************************************************************/
void
initialize_domain_info(domain_info_struct *info)
{
    strcpy(info->lat_var, "MISSING");
    strcpy(info->lon_var, "MISSING");
    strcpy(info->mask_var, "MISSING");
    strcpy(info->area_var, "MISSING");
    strcpy(info->frac_var, "MISSING");
    strcpy(info->y_dim, "MISSING");
    strcpy(info->x_dim, "MISSING");
    info->n_coord_dims = 0;
}

/******************************************************************************
 * @brief    Initialize global structures
 *****************************************************************************/
void
initialize_global_structures(void)
{
    extern domain_struct global_domain;
    extern domain_struct local_domain;
    extern int           mpi_rank;

    initialize_domain_info(&local_domain.info);
    if (mpi_rank == VIC_MPI_ROOT) {
        initialize_options();
        initialize_global();
        initialize_parameters();
        initialize_filenames();
        initialize_domain_info(&global_domain.info);
        initialize_domain(&global_domain);
    }
}
