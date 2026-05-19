/******************************************************************************
 * @section DESCRIPTION
 *
 * Allocate memory for VIC structures.
 *****************************************************************************/

#include "vic_driver_shared_image.h"

/******************************************************************************
 * @brief    Allocate memory for VIC structures.
 *****************************************************************************/
void
vic_alloc(void)
{
    extern all_vars_struct    *all_vars;
    extern force_data_struct  *force;
    extern domain_struct       local_domain;
    extern option_struct       options;
    extern double           ***out_data;
    extern save_data_struct   *save_data;
    extern soil_con_struct    *soil_con;
    extern veg_con_struct    **veg_con;
    extern veg_hist_struct   **veg_hist;
    size_t                     i, j;

    // allocate memory for force structure
    force = malloc(local_domain.ncells_active * sizeof(*force));
    check_alloc_status(force, "Memory allocation error.");

    // allocate memory for veg_hist structure
    veg_hist = malloc(local_domain.ncells_active * sizeof(*veg_hist));
    check_alloc_status(veg_hist, "Memory allocation error.");

    // allocate memory for soil structure
    soil_con = malloc(local_domain.ncells_active * sizeof(*soil_con));
    check_alloc_status(soil_con, "Memory allocation error.");

    // allocate memory for vegetation structure
    veg_con = malloc(local_domain.ncells_active * sizeof(*veg_con));
    check_alloc_status(veg_con, "Memory allocation error.");

    // all_vars allocation
    all_vars = malloc(local_domain.ncells_active * sizeof(*all_vars));
    check_alloc_status(all_vars, "Memory allocation error.");

    // out_data allocation
    out_data = malloc(local_domain.ncells_active * sizeof(*out_data));
    check_alloc_status(out_data, "Memory allocation error.");

    // save_data allocation
    save_data = calloc(local_domain.ncells_active, sizeof(*save_data));
    check_alloc_status(save_data, "Memory allocation error.");

    // allocate memory for individual grid cells
    for (i = 0; i < local_domain.ncells_active; i++) {
        // force allocation - allocate enough memory for NR+1 steps
        alloc_force(&(force[i]));

        // snow band allocation
        soil_con[i].AreaFract = calloc(options.SNOW_BAND,
                                       sizeof(*(soil_con[i].AreaFract)));
        check_alloc_status(soil_con[i].AreaFract, "Memory allocation error.");
        soil_con[i].BandElev = calloc(options.SNOW_BAND,
                                      sizeof(*(soil_con[i].BandElev)));
        check_alloc_status(soil_con[i].BandElev, "Memory allocation error.");
        soil_con[i].Tfactor = calloc(options.SNOW_BAND,
                                     sizeof(*(soil_con[i].Tfactor)));
        check_alloc_status(soil_con[i].Tfactor, "Memory allocation error.");
        soil_con[i].Pfactor = calloc(options.SNOW_BAND,
                                     sizeof(*(soil_con[i].Pfactor)));
        check_alloc_status(soil_con[i].Pfactor, "Memory allocation error.");

        initialize_soil_con(&(soil_con[i]));

        // vegetation tile allocation
        size_t nv_active = local_domain.locations[i].nveg + 1;
    
        veg_con[i] = malloc((nv_active) * sizeof(*(veg_con[i])));
        check_alloc_status(veg_con[i], "Memory allocation error.");

        for (j = 0; j < nv_active; j++) {
            if (options.CARBON) {
                veg_con[i][j].CanopLayerBnd = calloc(options.Ncanopy,
                                                     sizeof(*(veg_con[i][j].
                                                              CanopLayerBnd)));
                check_alloc_status(veg_con[i][j].CanopLayerBnd,
                                   "Memory allocation error.");
            }
            initialize_veg_con(&(veg_con[i][j]));
        }

        all_vars[i] = make_all_vars(nv_active - 1);

        // allocate memory for veg_hist
        veg_hist[i] = calloc(nv_active, sizeof(*(veg_hist[i])));
        for (j = 0; j < nv_active; j++) {
            alloc_veg_hist(&(veg_hist[i][j]));
        }
    }
}
