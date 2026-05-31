/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine sets the metadata structure for VIC state variables
 *****************************************************************************/

#include <vic_driver_shared_image.h>
#include <rout.h>

/******************************************************************************
 * @brief    Set output met data information
 *****************************************************************************/
void
set_state_meta_data_info()
{
    size_t                 v;

    extern option_struct   options;
    extern metadata_struct state_metadata[N_STATE_VARS + N_STATE_VARS_EXT];

    // Build the list of state variables

    // Set missing and/or default values
    for (v = 0; v < (N_STATE_VARS + N_STATE_VARS_EXT); v++) {
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
    strcpy(state_metadata[STATE_SOIL_MOISTURE].units, "mm");
    strcpy(state_metadata[STATE_SOIL_MOISTURE].description,
           "soil total moisture contents including ice for each soil layer");

    // STATE_SOIL_ICE
    strcpy(state_metadata[STATE_SOIL_ICE].varname, "STATE_SOIL_ICE");
    strcpy(state_metadata[STATE_SOIL_ICE].long_name, "soil_ice");
    strcpy(state_metadata[STATE_SOIL_ICE].standard_name,
           "soil_moisture_ice_depth");
    strcpy(state_metadata[STATE_SOIL_ICE].units, "mm");
    strcpy(state_metadata[STATE_SOIL_ICE].description,
           "soil ice content for each soil layer");

    // STATE_CANOPY_WATER
    strcpy(state_metadata[STATE_CANOPY_WATER].varname, "STATE_CANOPY_WATER");
    strcpy(state_metadata[STATE_CANOPY_WATER].long_name, "canopy_water");
    strcpy(state_metadata[STATE_CANOPY_WATER].standard_name, "water_in_canopy");
    strcpy(state_metadata[STATE_CANOPY_WATER].units, "mm");
    strcpy(state_metadata[STATE_CANOPY_WATER].description,
           "amount of water stored in the vegetation canopy");

    // STATE_TCANOPY
    strcpy(state_metadata[STATE_TCANOPY].varname, "STATE_TCANOPY");
    strcpy(state_metadata[STATE_TCANOPY].long_name, "Tcanopy");
    strcpy(state_metadata[STATE_TCANOPY].standard_name, "Tcanopy");
    strcpy(state_metadata[STATE_TCANOPY].units, "k");
    strcpy(state_metadata[STATE_TCANOPY].description,
           "temperature of the canopy");

    if (options.CARBON) {
        // STATE_ANNUALNPP
        strcpy(state_metadata[STATE_ANNUALNPP].varname, "STATE_ANNUALNPP");
        strcpy(state_metadata[STATE_ANNUALNPP].long_name, "annualnpp");
        strcpy(state_metadata[STATE_ANNUALNPP].standard_name,
               "running_total_annual_NPP");
        strcpy(state_metadata[STATE_ANNUALNPP].units, "g m-2");
        strcpy(state_metadata[STATE_ANNUALNPP].description,
               "running total annual NPP");

        // STATE_ANNUALNPPPREV
        strcpy(state_metadata[STATE_ANNUALNPPPREV].varname,
               "STATE_ANNUALNPPPREV");
        strcpy(state_metadata[STATE_ANNUALNPPPREV].long_name, "annualnppprev");
        strcpy(state_metadata[STATE_ANNUALNPPPREV].standard_name,
               "previous_year_total_annual_NPP");
        strcpy(state_metadata[STATE_ANNUALNPPPREV].units, "g m-2");
        strcpy(state_metadata[STATE_ANNUALNPPPREV].description,
               "total annual NPP from previous year");

        // STATE_CLITTER
        strcpy(state_metadata[STATE_CLITTER].varname, "STATE_CLITTER");
        strcpy(state_metadata[STATE_CLITTER].long_name, "clitter");
        strcpy(state_metadata[STATE_CLITTER].standard_name,
               "carbon_in_litter_pool");
        strcpy(state_metadata[STATE_CLITTER].units, "g m-2");
        strcpy(state_metadata[STATE_CLITTER].description,
               "carbon storage in litter pool");

        // STATE_CINTER
        strcpy(state_metadata[STATE_CINTER].varname, "STATE_CINTER");
        strcpy(state_metadata[STATE_CINTER].long_name, "cinter");
        strcpy(state_metadata[STATE_CINTER].standard_name,
               "carbon_in_intermediate_pool");
        strcpy(state_metadata[STATE_CINTER].units, "g m-2");
        strcpy(state_metadata[STATE_CINTER].description,
               "carbon storage in intermediate pool");

        // STATE_CSLOW
        strcpy(state_metadata[STATE_CSLOW].varname, "STATE_CSLOW");
        strcpy(state_metadata[STATE_CSLOW].long_name, "cslow");
        strcpy(state_metadata[STATE_CSLOW].standard_name,
               "carbon_in_slow_pool");
        strcpy(state_metadata[STATE_CSLOW].units, "g m-2");
        strcpy(state_metadata[STATE_CSLOW].description,
               "carbon storage in slow pool");
    }
    
    // STATE_SNOW_AGE
    strcpy(state_metadata[STATE_SNOW_AGE].varname, "STATE_SNOW_AGE");
    strcpy(state_metadata[STATE_SNOW_AGE].long_name, "snow_age");
    strcpy(state_metadata[STATE_SNOW_AGE].standard_name,
           "age_since_last_new_snow");
    strcpy(state_metadata[STATE_SNOW_AGE].units, "model_time_step");
    strcpy(state_metadata[STATE_SNOW_AGE].description,
           "number of model time steps since the last new snow");

    // STATE_SNOW_OLDSWQ
    strcpy(state_metadata[STATE_SNOW_OLDSWQ].varname,
           "STATE_SNOW_OLDSWQ");
    strcpy(state_metadata[STATE_SNOW_OLDSWQ].long_name, "last_step_swq");
    strcpy(state_metadata[STATE_SNOW_OLDSWQ].standard_name,
           "last_step_swq");
    strcpy(state_metadata[STATE_SNOW_OLDSWQ].units,
           "mm");
    strcpy(state_metadata[STATE_SNOW_OLDSWQ].description,
           "last step swq: snow[veg].old_swq");

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
    strcpy(state_metadata[STATE_SNOW_WATER_EQUIVALENT].units, "m");
    strcpy(state_metadata[STATE_SNOW_WATER_EQUIVALENT].description,
           "snow water equivalent");

    // STATE_SNOW_PACK_TEMP
    strcpy(state_metadata[STATE_SNOW_PACK_TEMP].varname,
           "STATE_SNOW_PACK_TEMP");
    strcpy(state_metadata[STATE_SNOW_PACK_TEMP].long_name, "snow_pack_temp");
    strcpy(state_metadata[STATE_SNOW_PACK_TEMP].standard_name,
           "snow_pack_temperature");
    strcpy(state_metadata[STATE_SNOW_PACK_TEMP].units, "C");
    strcpy(state_metadata[STATE_SNOW_PACK_TEMP].description,
           "snow pack layer temperature");

    // STATE_SNOW_THETA_LIQ
    strcpy(state_metadata[STATE_SNOW_THETA_LIQ].varname,
           "STATE_SNOW_THETA_LIQ");
    strcpy(state_metadata[STATE_SNOW_THETA_LIQ].long_name, "snow_theta_liq");
    strcpy(state_metadata[STATE_SNOW_THETA_LIQ].standard_name,
           "liquid_water_content_of_surface_snow");
    strcpy(state_metadata[STATE_SNOW_THETA_LIQ].units, "-");
    strcpy(state_metadata[STATE_SNOW_THETA_LIQ].description,
           "volumetric liquid water content in snow pack");

    // STATE_SNOW_THETA_ICE
    strcpy(state_metadata[STATE_SNOW_THETA_ICE].varname,
           "STATE_SNOW_THETA_ICE");
    strcpy(state_metadata[STATE_SNOW_THETA_ICE].long_name, "snow_theta_ice");
    strcpy(state_metadata[STATE_SNOW_THETA_ICE].standard_name,
           "ice_content_of_surface_snow");
    strcpy(state_metadata[STATE_SNOW_THETA_ICE].units, "-");
    strcpy(state_metadata[STATE_SNOW_THETA_ICE].description,
           "volumetric ice content in snow pack");

    // STATE_SNOW_PACK_ICE
    strcpy(state_metadata[STATE_SNOW_PACK_ICE].varname,
           "STATE_SNOW_PACK_ICE");
    strcpy(state_metadata[STATE_SNOW_PACK_ICE].long_name, "snow_pack_ice");
    strcpy(state_metadata[STATE_SNOW_PACK_ICE].standard_name,
           "snow_pack_ice_content");
    strcpy(state_metadata[STATE_SNOW_PACK_ICE].units, "m");
    strcpy(state_metadata[STATE_SNOW_PACK_ICE].description,
           "snow pack ice content expressed as equivalent water depth");

    // STATE_SNOW_PACK_LIQ
    strcpy(state_metadata[STATE_SNOW_PACK_LIQ].varname,
           "STATE_SNOW_PACK_LIQ");
    strcpy(state_metadata[STATE_SNOW_PACK_LIQ].long_name, "snow_pack_liq");
    strcpy(state_metadata[STATE_SNOW_PACK_LIQ].standard_name,
           "snow_pack_liquid_water_content");
    strcpy(state_metadata[STATE_SNOW_PACK_LIQ].units, "m");
    strcpy(state_metadata[STATE_SNOW_PACK_LIQ].description,
           "snow pack liquid water content expressed as equivalent water depth");

    // STATE_SNOW_POROSITY
    strcpy(state_metadata[STATE_SNOW_POROSITY].varname,
           "STATE_SNOW_POROSITY");
    strcpy(state_metadata[STATE_SNOW_POROSITY].long_name, "snow_porosity");
    strcpy(state_metadata[STATE_SNOW_POROSITY].standard_name,
           "snow_pack_porosity");
    strcpy(state_metadata[STATE_SNOW_POROSITY].units, "-");
    strcpy(state_metadata[STATE_SNOW_POROSITY].description,
           "porosity of snow pack");

    // STATE_SNOW_DENSITY
    strcpy(state_metadata[STATE_SNOW_DENSITY].varname, "STATE_SNOW_DENSITY");
    strcpy(state_metadata[STATE_SNOW_DENSITY].long_name, "snow_density");
    strcpy(state_metadata[STATE_SNOW_DENSITY].standard_name,
           "snowpack_density");
    strcpy(state_metadata[STATE_SNOW_DENSITY].units, "kg m-3");
    strcpy(state_metadata[STATE_SNOW_DENSITY].description, "snowpack density");

    // STATE_SNOW_NSNOW
    strcpy(state_metadata[STATE_SNOW_NSNOW].varname,
           "STATE_SNOW_NSNOW");
    strcpy(state_metadata[STATE_SNOW_NSNOW].long_name,
           "number_of_snow_layers");
    strcpy(state_metadata[STATE_SNOW_NSNOW].standard_name,
           "number_of_snow_layers");
    strcpy(state_metadata[STATE_SNOW_NSNOW].units, "int[0,1,2,3]");
    strcpy(state_metadata[STATE_SNOW_NSNOW].description,
           "Number of snow layers in the model");

    // STATE_SNOW_CANOPY
    strcpy(state_metadata[STATE_SNOW_CANOPY].varname, "STATE_SNOW_CANOPY");
    strcpy(state_metadata[STATE_SNOW_CANOPY].long_name, "snow_canopy");
    strcpy(state_metadata[STATE_SNOW_CANOPY].standard_name,
           "snow_water_equivalent_intercepted_in_canopy");
    strcpy(state_metadata[STATE_SNOW_CANOPY].units, "m");
    strcpy(state_metadata[STATE_SNOW_CANOPY].description,
           "snow interception storage in canopy");

    // STATE_NODE_TEMP
    strcpy(state_metadata[STATE_NODE_TEMP].varname,
           "STATE_NODE_TEMP");
    strcpy(state_metadata[STATE_NODE_TEMP].long_name, "node_temp");
    strcpy(state_metadata[STATE_NODE_TEMP].standard_name,
           "node_temperature");
    strcpy(state_metadata[STATE_NODE_TEMP].units, "k");
    strcpy(state_metadata[STATE_NODE_TEMP].description,
           "soil temperature of each soil thermal node");

    // STATE_FOLIAGE_TEMPERATURE
    strcpy(state_metadata[STATE_FOLIAGE_TEMPERATURE].varname,
           "STATE_FOLIAGE_TEMPERATURE");
    strcpy(state_metadata[STATE_FOLIAGE_TEMPERATURE].long_name,
           "foliage_temperature");
    strcpy(state_metadata[STATE_FOLIAGE_TEMPERATURE].standard_name,
           "foliage_temperature");
    strcpy(state_metadata[STATE_FOLIAGE_TEMPERATURE].units, "k");
    strcpy(state_metadata[STATE_FOLIAGE_TEMPERATURE].description,
           "overstory vegetaion temperature");

    // STATE_ROUT_RING
    state_metadata_rout_extension();
}
