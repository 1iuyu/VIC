/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes the snow variable arrays for each new grid cell.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    Initialize the snow variable arrays for each new grid cell.
 *****************************************************************************/
void
initialize_snow(snow_data_struct *snow,
                size_t            veg_num)
{
    size_t               i;

    for (i = 0; i <= veg_num; i++) {
            // Prognostic states
            snow[i].albedo = 0.0;
            snow[i].canopy_albedo = 0.0;
            snow[i].coldcontent = 0.0;
            snow[i].coverage = 0.0;
            snow[i].density = 0.0;
            snow[i].depth = 0.0;
            snow[i].last_snow = 0;
            snow[i].max_snow_depth = 0.0;
            snow[i].MELTING = false;
            snow[i].pack_temp = 0.0;
            snow[i].pack_water = 0.0;
            snow[i].snow = false;
            snow[i].snow_canopy = 0.0;
            snow[i].store_coverage = 0.0;
            snow[i].store_snow = false;
            snow[i].store_swq = 0.0;
            snow[i].surf_temp = 0.0;
            snow[i].surf_temp_fbcount = 0;
            snow[i].surf_temp_fbflag = false;
            snow[i].surf_water = 0.0;
            snow[i].swq = 0.0;
            snow[i].snow_distrib_slope = 0.0;
            snow[i].tmp_int_storage = 0.0;
            // Fluxes
            snow[i].blowing_flux = 0.0;
            snow[i].canopy_vapor_flux = 0.0;
            snow[i].mass_error = 0.0;
            snow[i].melt = 0.0;
            snow[i].Qnet = 0.0;
            snow[i].surface_flux = 0.0;
            snow[i].transport = 0.0;
            snow[i].vapor_flux = 0.0;
    }
}
