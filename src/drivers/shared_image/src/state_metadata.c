/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine sets the metadata structure for VIC state variables
 *****************************************************************************/

#include "vic_driver_shared_image.h"
#include "rout.h"

/******************************************************************************
 * @brief    Set output met data information
 *****************************************************************************/
void
set_state_meta_data_info()
{
    size_t                 v;

    extern option_struct   options;
    extern metadata_struct state_metadata[N_STATE_VARS];

    // Build the list of state variables

    // Set missing and/or default values
    for (v = 0; v < N_STATE_VARS; v++) {
        // Set default string values
        strcpy(state_metadata[v].varname, MISSING_S);
        strcpy(state_metadata[v].long_name, MISSING_S);
        strcpy(state_metadata[v].standard_name, MISSING_S);
        strcpy(state_metadata[v].units, MISSING_S);
        strcpy(state_metadata[v].description, MISSING_S);
        // Set default number of elements
        state_metadata[v].nelem = 1;
    }

    // STATE_SOIL_MOISTURE
    strcpy(state_metadata[STATE_SOIL_MOISTURE].varname, "STATE_SOIL_MOISTURE");
    strcpy(state_metadata[STATE_SOIL_MOISTURE].long_name, "soil_moisture");
    strcpy(state_metadata[STATE_SOIL_MOISTURE].standard_name,
           "soil_layer_moisture");
    strcpy(state_metadata[STATE_SOIL_MOISTURE].units, "m3/m3");
    strcpy(state_metadata[STATE_SOIL_MOISTURE].description,
           "soil total moisture contents including ice for each soil layer");

    // STATE_SOIL_ICE
    strcpy(state_metadata[STATE_SOIL_ICE].varname, "STATE_SOIL_ICE");
    strcpy(state_metadata[STATE_SOIL_ICE].long_name, "soil_ice");
    strcpy(state_metadata[STATE_SOIL_ICE].standard_name,
           "soil_ice_moisture");
    strcpy(state_metadata[STATE_SOIL_ICE].units, "m3/m3");
    strcpy(state_metadata[STATE_SOIL_ICE].description,
           "soil ice content for each soil layer");

    // STATE_SOIL_LIQ
    strcpy(state_metadata[STATE_SOIL_LIQ].varname, "STATE_SOIL_LIQ");
    strcpy(state_metadata[STATE_SOIL_LIQ].long_name, "soil_ice");
    strcpy(state_metadata[STATE_SOIL_LIQ].standard_name,
           "soil_water_moisture");
    strcpy(state_metadata[STATE_SOIL_LIQ].units, "m3/m3");
    strcpy(state_metadata[STATE_SOIL_LIQ].description,
           "soil water content for each soil layer");
    
    // STATE_SNOW_AGE
    strcpy(state_metadata[STATE_SNOW_AGE].varname, "STATE_SNOW_AGE");
    strcpy(state_metadata[STATE_SNOW_AGE].long_name, "snow_age");
    strcpy(state_metadata[STATE_SNOW_AGE].standard_name,
           "age_since_last_new_snow");
    strcpy(state_metadata[STATE_SNOW_AGE].units, "model_time_step");
    strcpy(state_metadata[STATE_SNOW_AGE].description,
           "number of model time steps since the last new snow");

    // STATE_SNOW_LASTSWQ
    strcpy(state_metadata[STATE_SNOW_LASTSWQ].varname,
           "STATE_SNOW_LASTSWQ");
    strcpy(state_metadata[STATE_SNOW_LASTSWQ].long_name, "last_step_swq");
    strcpy(state_metadata[STATE_SNOW_LASTSWQ].standard_name,
           "last_step_swq");
    strcpy(state_metadata[STATE_SNOW_LASTSWQ].units,
           "mm");
    strcpy(state_metadata[STATE_SNOW_LASTSWQ].description,
           "snow water equivalent at previous time step");

    // STATE_SNOW_COVERAGE
    strcpy(state_metadata[STATE_SNOW_COVERAGE].varname, "STATE_SNOW_COVERAGE");
    strcpy(state_metadata[STATE_SNOW_COVERAGE].long_name, "snow_coverage");
    strcpy(state_metadata[STATE_SNOW_COVERAGE].standard_name,
           "snow_coverage_fraction");
    strcpy(state_metadata[STATE_SNOW_COVERAGE].units, "1");
    strcpy(state_metadata[STATE_SNOW_COVERAGE].description,
           "fraction of grid cell area covered by snow");

    // STATE_SNOW_WATER_EQUIVALENT
    strcpy(state_metadata[STATE_SNOW_WATER_EQUIVALENT].varname,
           "STATE_SNOW_WATER_EQUIVALENT");
    strcpy(state_metadata[STATE_SNOW_WATER_EQUIVALENT].long_name,
           "snow_water_equivalent");
    strcpy(state_metadata[STATE_SNOW_WATER_EQUIVALENT].standard_name,
           "snow_water_equivalent");
    strcpy(state_metadata[STATE_SNOW_WATER_EQUIVALENT].units, "mm");
    strcpy(state_metadata[STATE_SNOW_WATER_EQUIVALENT].description,
           "snow water equivalent");

    // STATE_SNOW_RADIUS
    strcpy(state_metadata[STATE_SNOW_RADIUS].varname, "STATE_SNOW_RADIUS");
    strcpy(state_metadata[STATE_SNOW_RADIUS].long_name, "snow_radius");
    strcpy(state_metadata[STATE_SNOW_RADIUS].standard_name, "snow_grain_radius");
    strcpy(state_metadata[STATE_SNOW_RADIUS].units, "m-6");
    strcpy(state_metadata[STATE_SNOW_RADIUS].description, "effective grain radius");

    // STATE_SNOW_DZNODE
    strcpy(state_metadata[STATE_SNOW_DZNODE].varname, "STATE_SNOW_DZNODE");
    strcpy(state_metadata[STATE_SNOW_DZNODE].long_name, "snow_dz_node");
    strcpy(state_metadata[STATE_SNOW_DZNODE].standard_name, "snow_layer_thickness");
    strcpy(state_metadata[STATE_SNOW_DZNODE].units, "m");
    strcpy(state_metadata[STATE_SNOW_DZNODE].description, "each snow pack depth");

    // STATE_SNOW_LASTICE
    strcpy(state_metadata[STATE_SNOW_LASTICE].varname, "STATE_SNOW_LASTICE");
    strcpy(state_metadata[STATE_SNOW_LASTICE].long_name, "snow_last_ice");
    strcpy(state_metadata[STATE_SNOW_LASTICE].standard_name, "snow_ice_previous_step");
    strcpy(state_metadata[STATE_SNOW_LASTICE].units, "m3/m3");
    strcpy(state_metadata[STATE_SNOW_LASTICE].description, "partial volume of snow ice from previous time step");

    // STATE_SNOW_LASTLIQ
    strcpy(state_metadata[STATE_SNOW_LASTLIQ].varname, "STATE_SNOW_LASTLIQ");
    strcpy(state_metadata[STATE_SNOW_LASTLIQ].long_name, "snow_last_liq");
    strcpy(state_metadata[STATE_SNOW_LASTLIQ].standard_name, "snow_liquid_previous_step");
    strcpy(state_metadata[STATE_SNOW_LASTLIQ].units, "m3/m3");
    strcpy(state_metadata[STATE_SNOW_LASTLIQ].description, "partial volume of snow liquid water from previous time step");

    // STATE_INT_SNOW
    strcpy(state_metadata[STATE_INT_SNOW].varname, "STATE_INT_SNOW");
    strcpy(state_metadata[STATE_INT_SNOW].long_name, "intercepted_snow");
    strcpy(state_metadata[STATE_INT_SNOW].standard_name, "snow_intercepted");
    strcpy(state_metadata[STATE_INT_SNOW].units, "mm");
    strcpy(state_metadata[STATE_INT_SNOW].description, "snow intercepted on canopy");

    // STATE_INT_RAIN
    strcpy(state_metadata[STATE_INT_RAIN].varname, "STATE_INT_RAIN");
    strcpy(state_metadata[STATE_INT_RAIN].long_name, "intercepted_rain");
    strcpy(state_metadata[STATE_INT_RAIN].standard_name, "rain_intercepted");
    strcpy(state_metadata[STATE_INT_RAIN].units, "mm");
    strcpy(state_metadata[STATE_INT_RAIN].description, "rain intercepted on canopy");

    // STATE_VEG_MATRIC
    strcpy(state_metadata[STATE_VEG_MATRIC].varname, "STATE_VEG_MATRIC");
    strcpy(state_metadata[STATE_VEG_MATRIC].long_name, "vegetation_matric_potential");
    strcpy(state_metadata[STATE_VEG_MATRIC].standard_name, "vegetation_water_matric_potential");
    strcpy(state_metadata[STATE_VEG_MATRIC].units, "m");
    strcpy(state_metadata[STATE_VEG_MATRIC].description, "vegetation water matric potential [sun, shade, xylem, root]");

    // STATE_LAST_TEMP
    strcpy(state_metadata[STATE_LAST_TEMP].varname, "STATE_LAST_TEMP");
    strcpy(state_metadata[STATE_LAST_TEMP].long_name, "last_temperature");
    strcpy(state_metadata[STATE_LAST_TEMP].standard_name, "soil_temperature_previous_step");
    strcpy(state_metadata[STATE_LAST_TEMP].units, "k");
    strcpy(state_metadata[STATE_LAST_TEMP].description, "soil_temperature_previous_step");

    // STATE_MATRIC
    strcpy(state_metadata[STATE_MATRIC].varname, "STATE_MATRIC");
    strcpy(state_metadata[STATE_MATRIC].long_name, "matric_potential");
    strcpy(state_metadata[STATE_MATRIC].standard_name, "soil_matric_potential");
    strcpy(state_metadata[STATE_MATRIC].units, "m");
    strcpy(state_metadata[STATE_MATRIC].description, "soil matric potential");

    // STATE_ZWT
    strcpy(state_metadata[STATE_ZWT].varname, "STATE_ZWT");
    strcpy(state_metadata[STATE_ZWT].long_name, "zwt");
    strcpy(state_metadata[STATE_ZWT].standard_name, "water_table_depth");
    strcpy(state_metadata[STATE_ZWT].units, "m");
    strcpy(state_metadata[STATE_ZWT].description, "water table depth");

    // STATE_AQF_STORAGE
    strcpy(state_metadata[STATE_AQF_STORAGE].varname, "STATE_AQF_STORAGE");
    strcpy(state_metadata[STATE_AQF_STORAGE].long_name, "aqf_storage");
    strcpy(state_metadata[STATE_AQF_STORAGE].standard_name, "aquifer_storage");
    strcpy(state_metadata[STATE_AQF_STORAGE].units, "m");
    strcpy(state_metadata[STATE_AQF_STORAGE].description, "water storage in aquifer");

    // STATE_SNOW_PACK_ICE
    strcpy(state_metadata[STATE_SNOW_PACK_ICE].varname,
           "STATE_SNOW_PACK_ICE");
    strcpy(state_metadata[STATE_SNOW_PACK_ICE].long_name, "snow_pack_ice");
    strcpy(state_metadata[STATE_SNOW_PACK_ICE].standard_name,
           "snow_pack_ice_content");
    strcpy(state_metadata[STATE_SNOW_PACK_ICE].units, "mm");
    strcpy(state_metadata[STATE_SNOW_PACK_ICE].description,
           "snow pack ice content expressed as equivalent water depth");

    // STATE_SNOW_PACK_LIQ
    strcpy(state_metadata[STATE_SNOW_PACK_LIQ].varname,
           "STATE_SNOW_PACK_LIQ");
    strcpy(state_metadata[STATE_SNOW_PACK_LIQ].long_name, "snow_pack_liq");
    strcpy(state_metadata[STATE_SNOW_PACK_LIQ].standard_name,
           "snow_pack_liquid_water_content");
    strcpy(state_metadata[STATE_SNOW_PACK_LIQ].units, "mm");
    strcpy(state_metadata[STATE_SNOW_PACK_LIQ].description,
           "snow pack liquid water content expressed as equivalent water depth");

    // STATE_SNOW_DENSITY
    strcpy(state_metadata[STATE_SNOW_DENSITY].varname, "STATE_SNOW_DENSITY");
    strcpy(state_metadata[STATE_SNOW_DENSITY].long_name, "snow_density");
    strcpy(state_metadata[STATE_SNOW_DENSITY].standard_name,
           "snowpack_density");
    strcpy(state_metadata[STATE_SNOW_DENSITY].units, "kg m-3");
    strcpy(state_metadata[STATE_SNOW_DENSITY].description, "snowpack density");

    // STATE_SNOW_NSNOW
    strcpy(state_metadata[STATE_NSNOW].varname,
           "STATE_NSNOW");
    strcpy(state_metadata[STATE_NSNOW].long_name,
           "number_of_snow_layers");
    strcpy(state_metadata[STATE_NSNOW].standard_name,
           "number_of_snow_layers");
    strcpy(state_metadata[STATE_NSNOW].units, "1");
    strcpy(state_metadata[STATE_NSNOW].description,
           "Number of snow layers in the model");

    // STATE_NODE_TEMP
    strcpy(state_metadata[STATE_NODE_TEMP].varname,
           "STATE_NODE_TEMP");
    strcpy(state_metadata[STATE_NODE_TEMP].long_name, "node_temp");
    strcpy(state_metadata[STATE_NODE_TEMP].standard_name,
           "node_temperature");
    strcpy(state_metadata[STATE_NODE_TEMP].units, "k");
    strcpy(state_metadata[STATE_NODE_TEMP].description,
           "soil temperature of each soil thermal node");

    // STATE_NSOIL
    strcpy(state_metadata[STATE_NSOIL].varname, "STATE_NSOIL");
    strcpy(state_metadata[STATE_NSOIL].long_name, "nsoil");
    strcpy(state_metadata[STATE_NSOIL].standard_name, "number_soil_layers");
    strcpy(state_metadata[STATE_NSOIL].units, "1");
    strcpy(state_metadata[STATE_NSOIL].description, "Number of soil layers");

    // STATE_NNODE
    strcpy(state_metadata[STATE_NNODE].varname, "STATE_NNODE");
    strcpy(state_metadata[STATE_NNODE].long_name, "nnode");
    strcpy(state_metadata[STATE_NNODE].standard_name, "number_thermal_nodes");
    strcpy(state_metadata[STATE_NNODE].units, "1");
    strcpy(state_metadata[STATE_NNODE].description, "Number of thermal nodes");

    // STATE_NCANOPY
    strcpy(state_metadata[STATE_NCANOPY].varname, "STATE_NCANOPY");
    strcpy(state_metadata[STATE_NCANOPY].long_name, "ncanopy");
    strcpy(state_metadata[STATE_NCANOPY].standard_name, "number_canopy_layers");
    strcpy(state_metadata[STATE_NCANOPY].units, "1");
    strcpy(state_metadata[STATE_NCANOPY].description, "Number of canopy layers");

    // STATE_NROOT
    strcpy(state_metadata[STATE_NROOT].varname, "STATE_NROOT");
    strcpy(state_metadata[STATE_NROOT].long_name, "nroot");
    strcpy(state_metadata[STATE_NROOT].standard_name, "number_root_layers");
    strcpy(state_metadata[STATE_NROOT].units, "1");
    strcpy(state_metadata[STATE_NROOT].description, "Number of root layers");

    // STATE_SOIL_LASTICE
    strcpy(state_metadata[STATE_SOIL_LASTICE].varname, "STATE_SOIL_LASTICE");
    strcpy(state_metadata[STATE_SOIL_LASTICE].long_name, "soil_last_ice");
    strcpy(state_metadata[STATE_SOIL_LASTICE].standard_name, "soil_ice_previous_step");
    strcpy(state_metadata[STATE_SOIL_LASTICE].units, "m3/m3");
    strcpy(state_metadata[STATE_SOIL_LASTICE].description, "last step ice content of the soil sublayer");

    // STATE_SOIL_LASTLIQ
    strcpy(state_metadata[STATE_SOIL_LASTLIQ].varname, "STATE_SOIL_LASTLIQ");
    strcpy(state_metadata[STATE_SOIL_LASTLIQ].long_name, "soil_last_liq");
    strcpy(state_metadata[STATE_SOIL_LASTLIQ].standard_name, "soil_liquid_previous_step");
    strcpy(state_metadata[STATE_SOIL_LASTLIQ].units, "m3/m3");
    strcpy(state_metadata[STATE_SOIL_LASTLIQ].description, "last step liquid content of the soil sublayer");

    // STATE_LAST_MATRIC
    strcpy(state_metadata[STATE_LAST_MATRIC].varname, "STATE_LAST_MATRIC");
    strcpy(state_metadata[STATE_LAST_MATRIC].long_name, "last_matric");
    strcpy(state_metadata[STATE_LAST_MATRIC].standard_name, "matric_potential_previous_step");
    strcpy(state_metadata[STATE_LAST_MATRIC].units, "m");
    strcpy(state_metadata[STATE_LAST_MATRIC].description, "last step matric potential");

    // STATE_H2OSFC
    strcpy(state_metadata[STATE_H2OSFC].varname, "STATE_H2OSFC");
    strcpy(state_metadata[STATE_H2OSFC].long_name, "h2osfc");
    strcpy(state_metadata[STATE_H2OSFC].standard_name, "surface_water");
    strcpy(state_metadata[STATE_H2OSFC].units, "mm");
    strcpy(state_metadata[STATE_H2OSFC].description, "surface water or glacier");

    // STATE_H2O_FRAC
    strcpy(state_metadata[STATE_H2O_FRAC].varname, "STATE_H2O_FRAC");
    strcpy(state_metadata[STATE_H2O_FRAC].long_name, "h2o_frac");
    strcpy(state_metadata[STATE_H2O_FRAC].standard_name, "surface_water_fraction");
    strcpy(state_metadata[STATE_H2O_FRAC].units, "1");
    strcpy(state_metadata[STATE_H2O_FRAC].description, "fraction of ground covered by surface water or glacier");

    // STATE_H2OSFC_ICE
    strcpy(state_metadata[STATE_H2OSFC_ICE].varname, "STATE_H2OSFC_ICE");
    strcpy(state_metadata[STATE_H2OSFC_ICE].long_name, "h2osfc_ice");
    strcpy(state_metadata[STATE_H2OSFC_ICE].standard_name, "surface_water_ice");
    strcpy(state_metadata[STATE_H2OSFC_ICE].units, "mm");
    strcpy(state_metadata[STATE_H2OSFC_ICE].description, "surface water or glacier ice content");

    // STATE_H2OSFC_LIQ
    strcpy(state_metadata[STATE_H2OSFC_LIQ].varname, "STATE_H2OSFC_LIQ");
    strcpy(state_metadata[STATE_H2OSFC_LIQ].long_name, "h2osfc_liq");
    strcpy(state_metadata[STATE_H2OSFC_LIQ].standard_name, "surface_water_liquid");
    strcpy(state_metadata[STATE_H2OSFC_LIQ].units, "mm");
    strcpy(state_metadata[STATE_H2OSFC_LIQ].description, "surface water or glacier liquid content");

    // STATE_ROUT
    if (options.ROUT) {
        // STATE_MAIN_CHANNEL_STORAGE
        strcpy(state_metadata[STATE_MAIN_CHANNEL_STORAGE].varname, "STATE_MAIN_CHANNEL_STORAGE");
        strcpy(state_metadata[STATE_MAIN_CHANNEL_STORAGE].long_name, "main_channel_storage");
        strcpy(state_metadata[STATE_MAIN_CHANNEL_STORAGE].standard_name, "main_channel_water_storage");
        strcpy(state_metadata[STATE_MAIN_CHANNEL_STORAGE].units, "m3");
        strcpy(state_metadata[STATE_MAIN_CHANNEL_STORAGE].description, "main channel storage");

        // STATE_MAIN_CROSS_SECTION_AREA
        strcpy(state_metadata[STATE_MAIN_CROSS_SECTION_AREA].varname, "STATE_MAIN_CROSS_SECTION_AREA");
        strcpy(state_metadata[STATE_MAIN_CROSS_SECTION_AREA].long_name, "main_cross_section_area");
        strcpy(state_metadata[STATE_MAIN_CROSS_SECTION_AREA].standard_name, "main_channel_cross_section_area");
        strcpy(state_metadata[STATE_MAIN_CROSS_SECTION_AREA].units, "m2");
        strcpy(state_metadata[STATE_MAIN_CROSS_SECTION_AREA].description, "main channel cross section area");

        // STATE_MAIN_CHANNEL_DEPTH
        strcpy(state_metadata[STATE_MAIN_CHANNEL_DEPTH].varname, "STATE_MAIN_CHANNEL_DEPTH");
        strcpy(state_metadata[STATE_MAIN_CHANNEL_DEPTH].long_name, "main_channel_depth");
        strcpy(state_metadata[STATE_MAIN_CHANNEL_DEPTH].standard_name, "main_channel_depth");
        strcpy(state_metadata[STATE_MAIN_CHANNEL_DEPTH].units, "m");
        strcpy(state_metadata[STATE_MAIN_CHANNEL_DEPTH].description, "main channel depth");

        // STATE_MAIN_CHANNEL_MANNING_N
        strcpy(state_metadata[STATE_MAIN_CHANNEL_MANNING_N].varname, "STATE_MAIN_CHANNEL_MANNING_N");
        strcpy(state_metadata[STATE_MAIN_CHANNEL_MANNING_N].long_name, "main_channel_manning_n");
        strcpy(state_metadata[STATE_MAIN_CHANNEL_MANNING_N].standard_name, "main_channel_manning_coefficient");
        strcpy(state_metadata[STATE_MAIN_CHANNEL_MANNING_N].units, "s/m^(1/3)");
        strcpy(state_metadata[STATE_MAIN_CHANNEL_MANNING_N].description, "main channel manning n");

        // STATE_MAIN_WETTED_PERIMETER
        strcpy(state_metadata[STATE_MAIN_WETTED_PERIMETER].varname, "STATE_MAIN_WETTED_PERIMETER");
        strcpy(state_metadata[STATE_MAIN_WETTED_PERIMETER].long_name, "main_wetted_perimeter");
        strcpy(state_metadata[STATE_MAIN_WETTED_PERIMETER].standard_name, "main_channel_wetted_perimeter");
        strcpy(state_metadata[STATE_MAIN_WETTED_PERIMETER].units, "m");
        strcpy(state_metadata[STATE_MAIN_WETTED_PERIMETER].description, "main channel wetted perimeter");

        // STATE_MAIN_HYDRAULIC_RADIUS
        strcpy(state_metadata[STATE_MAIN_HYDRAULIC_RADIUS].varname, "STATE_MAIN_HYDRAULIC_RADIUS");
        strcpy(state_metadata[STATE_MAIN_HYDRAULIC_RADIUS].long_name, "main_hydraulic_radius");
        strcpy(state_metadata[STATE_MAIN_HYDRAULIC_RADIUS].standard_name, "main_channel_hydraulic_radius");
        strcpy(state_metadata[STATE_MAIN_HYDRAULIC_RADIUS].units, "m");
        strcpy(state_metadata[STATE_MAIN_HYDRAULIC_RADIUS].description, "main channel hydraulic radius");

        // STATE_SUB_CHANNEL_STORAGE
        strcpy(state_metadata[STATE_SUB_CHANNEL_STORAGE].varname, "STATE_SUB_CHANNEL_STORAGE");
        strcpy(state_metadata[STATE_SUB_CHANNEL_STORAGE].long_name, "sub_channel_storage");
        strcpy(state_metadata[STATE_SUB_CHANNEL_STORAGE].standard_name, "sub_channel_water_storage");
        strcpy(state_metadata[STATE_SUB_CHANNEL_STORAGE].units, "m3");
        strcpy(state_metadata[STATE_SUB_CHANNEL_STORAGE].description, "sub channel storage");

        // STATE_SUB_CHANNEL_MANNING_N
        strcpy(state_metadata[STATE_SUB_CHANNEL_MANNING_N].varname, "STATE_SUB_CHANNEL_MANNING_N");
        strcpy(state_metadata[STATE_SUB_CHANNEL_MANNING_N].long_name, "sub_channel_manning_n");
        strcpy(state_metadata[STATE_SUB_CHANNEL_MANNING_N].standard_name, "sub_channel_manning_coefficient");
        strcpy(state_metadata[STATE_SUB_CHANNEL_MANNING_N].units, "s/m^(1/3)");
        strcpy(state_metadata[STATE_SUB_CHANNEL_MANNING_N].description, "sub channel manning n");

        // STATE_SUB_CROSS_SECTION_AREA
        strcpy(state_metadata[STATE_SUB_CROSS_SECTION_AREA].varname, "STATE_SUB_CROSS_SECTION_AREA");
        strcpy(state_metadata[STATE_SUB_CROSS_SECTION_AREA].long_name, "sub_cross_section_area");
        strcpy(state_metadata[STATE_SUB_CROSS_SECTION_AREA].standard_name, "sub_channel_cross_section_area");
        strcpy(state_metadata[STATE_SUB_CROSS_SECTION_AREA].units, "m2");
        strcpy(state_metadata[STATE_SUB_CROSS_SECTION_AREA].description, "sub channel cross section area");

        // STATE_SUB_CHANNEL_DEPTH
        strcpy(state_metadata[STATE_SUB_CHANNEL_DEPTH].varname, "STATE_SUB_CHANNEL_DEPTH");
        strcpy(state_metadata[STATE_SUB_CHANNEL_DEPTH].long_name, "sub_channel_depth");
        strcpy(state_metadata[STATE_SUB_CHANNEL_DEPTH].standard_name, "sub_channel_depth");
        strcpy(state_metadata[STATE_SUB_CHANNEL_DEPTH].units, "m");
        strcpy(state_metadata[STATE_SUB_CHANNEL_DEPTH].description, "sub channel depth");

        // STATE_SUB_WETTED_PERIMETER
        strcpy(state_metadata[STATE_SUB_WETTED_PERIMETER].varname, "STATE_SUB_WETTED_PERIMETER");
        strcpy(state_metadata[STATE_SUB_WETTED_PERIMETER].long_name, "sub_wetted_perimeter");
        strcpy(state_metadata[STATE_SUB_WETTED_PERIMETER].standard_name, "sub_channel_wetted_perimeter");
        strcpy(state_metadata[STATE_SUB_WETTED_PERIMETER].units, "m");
        strcpy(state_metadata[STATE_SUB_WETTED_PERIMETER].description, "sub channel wetted perimeter");

        // STATE_SUB_HYDRAULIC_RADIUS
        strcpy(state_metadata[STATE_SUB_HYDRAULIC_RADIUS].varname, "STATE_SUB_HYDRAULIC_RADIUS");
        strcpy(state_metadata[STATE_SUB_HYDRAULIC_RADIUS].long_name, "sub_hydraulic_radius");
        strcpy(state_metadata[STATE_SUB_HYDRAULIC_RADIUS].standard_name, "sub_channel_hydraulic_radius");
        strcpy(state_metadata[STATE_SUB_HYDRAULIC_RADIUS].units, "m");
        strcpy(state_metadata[STATE_SUB_HYDRAULIC_RADIUS].description, "sub channel hydraulic radius");

        // STATE_HILLSLOPE_DEPTH
        strcpy(state_metadata[STATE_HILLSLOPE_DEPTH].varname, "STATE_HILLSLOPE_DEPTH");
        strcpy(state_metadata[STATE_HILLSLOPE_DEPTH].long_name, "hillslope_depth");
        strcpy(state_metadata[STATE_HILLSLOPE_DEPTH].standard_name, "hillslope_water_depth");
        strcpy(state_metadata[STATE_HILLSLOPE_DEPTH].units, "m");
        strcpy(state_metadata[STATE_HILLSLOPE_DEPTH].description, "hillslope depth");

        // STATE_HILLSLOPE_MANNING_N
        strcpy(state_metadata[STATE_HILLSLOPE_MANNING_N].varname, "STATE_HILLSLOPE_MANNING_N");
        strcpy(state_metadata[STATE_HILLSLOPE_MANNING_N].long_name, "hillslope_manning_n");
        strcpy(state_metadata[STATE_HILLSLOPE_MANNING_N].standard_name, "hillslope_manning_coefficient");
        strcpy(state_metadata[STATE_HILLSLOPE_MANNING_N].units, "s/m^(1/3)");
        strcpy(state_metadata[STATE_HILLSLOPE_MANNING_N].description, "hillslope manning n");

        // STATE_HILLSLOPE_STORAGE
        strcpy(state_metadata[STATE_HILLSLOPE_STORAGE].varname, "STATE_HILLSLOPE_STORAGE");
        strcpy(state_metadata[STATE_HILLSLOPE_STORAGE].long_name, "hillslope_storage");
        strcpy(state_metadata[STATE_HILLSLOPE_STORAGE].standard_name, "hillslope_water_storage");
        strcpy(state_metadata[STATE_HILLSLOPE_STORAGE].units, "mm");
        strcpy(state_metadata[STATE_HILLSLOPE_STORAGE].description, "hillslope storage");

        // STATE_STORAGE_PREV
        strcpy(state_metadata[STATE_STORAGE_PREV].varname, "STATE_STORAGE_PREV");
        strcpy(state_metadata[STATE_STORAGE_PREV].long_name, "storage_prev");
        strcpy(state_metadata[STATE_STORAGE_PREV].standard_name, "previous_storage");
        strcpy(state_metadata[STATE_STORAGE_PREV].units, "m3");
        strcpy(state_metadata[STATE_STORAGE_PREV].description, "previous storage");
    }
}
