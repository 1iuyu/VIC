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
state_metadata_rout_extension()
{
    extern metadata_struct state_metadata[N_STATE_VARS + N_STATE_VARS_EXT];

    // STATE_MAIN_CHANNEL_STORAGE
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_CHANNEL_STORAGE].varname,
           "STATE_MAIN_CHANNEL_STORAGE");
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_CHANNEL_STORAGE].long_name,
           "main_channel_storage");
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_CHANNEL_STORAGE].standard_name,
           "main_channel_storage");
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_CHANNEL_STORAGE].units, "m^3");
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_CHANNEL_STORAGE].description,
           "water storage in main channel");

    // STATE_MAIN_MANNING_N
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_CHANNEL_MANNING_N].varname,
           "STATE_MAIN_CHANNEL_MANNING_N");
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_CHANNEL_MANNING_N].long_name,
           "main_channel_mannings_n");
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_CHANNEL_MANNING_N].standard_name,
           "main_channel_mannings_roughness_coefficient");
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_CHANNEL_MANNING_N].units, "-");
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_CHANNEL_MANNING_N].description,
           "Manning's roughness coefficient for main channel");

    // STATE_MAIN_CROSS_SECTION_AREA
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_CROSS_SECTION_AREA].varname,
           "STATE_MAIN_CROSS_SECTION_AREA");
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_CROSS_SECTION_AREA].long_name,
           "main_channel_cross_section_area");
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_CROSS_SECTION_AREA].standard_name,
           "main_channel_cross_section_area");
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_CROSS_SECTION_AREA].units, "m^2");
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_CROSS_SECTION_AREA].description,
           "main channel cross-sectional area");

    // STATE_MAIN_CHANNEL_DEPTH
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_CHANNEL_DEPTH].varname,
           "STATE_MAIN_CHANNEL_DEPTH");
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_CHANNEL_DEPTH].long_name,
           "main_channel_depth");
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_CHANNEL_DEPTH].standard_name,
           "main_channel_water_depth");
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_CHANNEL_DEPTH].units, "m");
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_CHANNEL_DEPTH].description,
           "water depth in main channel");

    // STATE_MAIN_WETTED_PERIMETER
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_WETTED_PERIMETER].varname,
           "STATE_MAIN_WETTED_PERIMETER");
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_WETTED_PERIMETER].long_name,
           "wetted_perimeter");
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_WETTED_PERIMETER].standard_name,
           "main_channel_wetted_perimeter");
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_WETTED_PERIMETER].units, "m");
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_WETTED_PERIMETER].description,
           "wetted perimeter of main channel cross-section");

    // STATE_MAIN_HYDRAULIC_RADIUS
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_HYDRAULIC_RADIUS].varname,
           "STATE_MAIN_HYDRAULIC_RADIUS");
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_HYDRAULIC_RADIUS].long_name,
           "hydraulic_radius");
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_HYDRAULIC_RADIUS].standard_name,
           "main_channel_hydraulic_radius");
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_HYDRAULIC_RADIUS].units, "m");
    strcpy(state_metadata[N_STATE_VARS + STATE_MAIN_HYDRAULIC_RADIUS].description,
           "hydraulic radius of main channel");

        // STATE_ROUT_HILLSLOPE_DEPTH
    strcpy(state_metadata[N_STATE_VARS + STATE_HILLSLOPE_DEPTH].varname,
           "STATE_ROUT_HILLSLOPE_DEPTH");
    strcpy(state_metadata[N_STATE_VARS + STATE_HILLSLOPE_DEPTH].long_name,
           "hillslope_flow_depth");
    strcpy(state_metadata[N_STATE_VARS + STATE_HILLSLOPE_DEPTH].standard_name,
           "hillslope_flow_depth");
    strcpy(state_metadata[N_STATE_VARS + STATE_HILLSLOPE_DEPTH].units, "m");
    strcpy(state_metadata[N_STATE_VARS + STATE_HILLSLOPE_DEPTH].description,
           "flow depth on hillslope");

    // STATE_ROUT_HILLSLOPE_MANNING_N
    strcpy(state_metadata[N_STATE_VARS + STATE_HILLSLOPE_MANNING_N].varname,
           "STATE_ROUT_HILLSLOPE_MANNING_N");
    strcpy(state_metadata[N_STATE_VARS + STATE_HILLSLOPE_MANNING_N].long_name,
           "hillslope_mannings_n");
    strcpy(state_metadata[N_STATE_VARS + STATE_HILLSLOPE_MANNING_N].standard_name,
           "hillslope_mannings_roughness_coefficient");
    strcpy(state_metadata[N_STATE_VARS + STATE_HILLSLOPE_MANNING_N].units, "-");
    strcpy(state_metadata[N_STATE_VARS + STATE_HILLSLOPE_MANNING_N].description,
           "Manning's roughness coefficient for hillslope");

    // STATE_ROUT_HILLSLOPE_STORAGE
    strcpy(state_metadata[N_STATE_VARS + STATE_HILLSLOPE_STORAGE].varname,
           "STATE_ROUT_HILLSLOPE_STORAGE");
    strcpy(state_metadata[N_STATE_VARS + STATE_HILLSLOPE_STORAGE].long_name,
           "hillslope_storage");
    strcpy(state_metadata[N_STATE_VARS + STATE_HILLSLOPE_STORAGE].standard_name,
           "hillslope_water_storage");
    strcpy(state_metadata[N_STATE_VARS + STATE_HILLSLOPE_STORAGE].units, "mm");
    strcpy(state_metadata[N_STATE_VARS + STATE_HILLSLOPE_STORAGE].description,
           "water storage on hillslope");
       // STATE_ROUT_SUB_CHANNEL_STORAGE
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_CHANNEL_STORAGE].varname,
           "STATE_SUB_CHANNEL_STORAGE");
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_CHANNEL_STORAGE].long_name,
           "sub_channel_storage");
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_CHANNEL_STORAGE].standard_name,
           "sub_channel_water_storage");
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_CHANNEL_STORAGE].units, "mm");
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_CHANNEL_STORAGE].description,
           "water storage in sub channel");

    // STATE_ROUT_SUB_CHANNEL_MANNING_N
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_CHANNEL_MANNING_N].varname,
           "STATE_SUB_CHANNEL_MANNING_N");
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_CHANNEL_MANNING_N].long_name,
           "sub_channel_mannings_n");
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_CHANNEL_MANNING_N].standard_name,
           "sub_channel_mannings_roughness_coefficient");
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_CHANNEL_MANNING_N].units, "-");
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_CHANNEL_MANNING_N].description,
           "Manning's roughness coefficient for sub channel");

    // STATE_ROUT_CROSS_SECTION_AREA
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_CROSS_SECTION_AREA].varname,
           "STATE_SUB_CROSS_SECTION_AREA");
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_CROSS_SECTION_AREA].long_name,
           "sub_channel_cross_section_area");
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_CROSS_SECTION_AREA].standard_name,
           "sub_channel_cross_section_area");
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_CROSS_SECTION_AREA].units, "m^2");
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_CROSS_SECTION_AREA].description,
           "sub channel cross-sectional area");

    // STATE_ROUT_CHANNEL_DEPTH
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_CHANNEL_DEPTH].varname,
           "STATE_SUB_CHANNEL_DEPTH");
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_CHANNEL_DEPTH].long_name,
           "sub_channel_depth");
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_CHANNEL_DEPTH].standard_name,
           "sub_channel_water_depth");
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_CHANNEL_DEPTH].units, "m");
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_CHANNEL_DEPTH].description,
           "water depth in sub channel");

    // STATE_ROUT_WETTED_PERIMETER
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_WETTED_PERIMETER].varname,
           "STATE_SUB_WETTED_PERIMETER");
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_WETTED_PERIMETER].long_name,
           "wetted_perimeter");
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_WETTED_PERIMETER].standard_name,
           "sub_channel_wetted_perimeter");
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_WETTED_PERIMETER].units, "m");
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_WETTED_PERIMETER].description,
           "wetted perimeter of sub channel cross-section");

    // STATE_ROUT_HYDRAULIC_RADIUS
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_HYDRAULIC_RADIUS].varname,
           "STATE_SUB_HYDRAULIC_RADIUS");
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_HYDRAULIC_RADIUS].long_name,
           "hydraulic_radius");
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_HYDRAULIC_RADIUS].standard_name,
           "sub_channel_hydraulic_radius");
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_HYDRAULIC_RADIUS].units, "m");
    strcpy(state_metadata[N_STATE_VARS + STATE_SUB_HYDRAULIC_RADIUS].description,
           "hydraulic radius of sub channel");
}

