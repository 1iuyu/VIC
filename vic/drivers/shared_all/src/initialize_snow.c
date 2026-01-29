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
    size_t      lidx, i;

    static const double dz_snow_min[MAX_SNOWS] = {
        0.025, 0.025, 0.1, 0., 0.
    };

    for (i = 0; i <= veg_num; i++) {
            // Prognostic states
            snow[i].albedo = 0.0;
            snow[i].coverage = 0.0;
            snow[i].density = 0.0;
            snow[i].SnowAge = 0.0;
            snow[i].old_swq = 0.0;
            snow[i].snow_depth = 0.0;
            snow[i].new_snow_density = 0.0;
            for (lidx = 0; lidx < MAX_SNOWS; lidx++) {
                snow[i].dz_snow[lidx] = 0.0;
                snow[i].Zsum_snow[lidx] = 0.0;
                snow[i].zc_snow[lidx] = 0.0;
                snow[i].pack_liq[lidx] = 0.0;
                snow[i].pack_ice[lidx] = 0.0;
                snow[i].theta_ice[lidx] = 0.0;
                snow[i].theta_liq[lidx] = 0.0;
                snow[i].pack_T[lidx] = 0.0;
                snow[i].pack_outflow[lidx] = 0.0;
                snow[i].porosity[lidx] = 0.0;
                snow[i].snow_thresholds[lidx] = dz_snow_min[lidx];
            }
            snow[i].glac_excess = 0.0;
            snow[i].pack_melt = 0.0;
            snow[i].pack_comb = 0.0;
            snow[i].pack_transp = 0.0;
            snow[i].delta_depth = 0.0;
            snow[i].snow_canopy = 0.0;
            snow[i].swq = 0.0;
            snow[i].Nsnow = 0;
            
            // Fluxes
    }
}
