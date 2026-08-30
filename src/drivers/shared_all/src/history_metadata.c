/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine sets the metadata structure for VIC output variables
 *****************************************************************************/

#include "vic_driver_shared_all.h"

/******************************************************************************
 * @brief    Set output met data information
 *****************************************************************************/
void
set_output_met_data_info()
{
    size_t                 v;

    extern option_struct   options;
    extern metadata_struct out_metadata[N_OUTVAR_TYPES];
    
    // Build the list of supported output variables
    // Set missing and/or default values
    for (v = 0; v < N_OUTVAR_TYPES; v++) {
       // Set default string values
       strcpy(out_metadata[v].varname, MISSING_S);
       strcpy(out_metadata[v].long_name, MISSING_S);
       strcpy(out_metadata[v].standard_name, MISSING_S);
       strcpy(out_metadata[v].units, MISSING_S);
       strcpy(out_metadata[v].description, MISSING_S);
       // Set default number of elements
       out_metadata[v].nelem = 1;
       out_metadata[v].elem_type = OUT_ELEM_DEFAULT;
    }

    // Water Balance Terms - state variables
    /* saturated area fraction */
    strcpy(out_metadata[OUT_ASAT].varname, "OUT_ASAT");
    strcpy(out_metadata[OUT_ASAT].long_name, "asat");
    strcpy(out_metadata[OUT_ASAT].standard_name, "saturated_area_fraction");
    strcpy(out_metadata[OUT_ASAT].units, "1");
    strcpy(out_metadata[OUT_ASAT].description, "saturated area fraction");

    /* river discharge [m3 s-1] */
    strcpy(out_metadata[OUT_DISCHARGE].varname, "OUT_DISCHARGE");
    strcpy(out_metadata[OUT_DISCHARGE].long_name,
           "water_volume_transport_in_river_channel");
    strcpy(out_metadata[OUT_DISCHARGE].standard_name, "river_discharge");
    strcpy(out_metadata[OUT_DISCHARGE].units, "m3 s-1");
    strcpy(
           out_metadata[OUT_DISCHARGE].description,
           "The water flux or volume transport in rivers is the amount of water flowing in the river channel and flood plain. 'Water' means water in all phases");
    
    /* root zone soil moisture [mm] */
    strcpy(out_metadata[OUT_ROOTMOIST].varname, "OUT_ROOTMOIST");
    strcpy(out_metadata[OUT_ROOTMOIST].long_name, "rootmoist");
    strcpy(out_metadata[OUT_ROOTMOIST].standard_name,
           "soil_moisture_in_root_zone");
    strcpy(out_metadata[OUT_ROOTMOIST].units, "mm");
    strcpy(out_metadata[OUT_ROOTMOIST].description, "root zone soil moisture");

    /* precipitation interception storage in canopy [mm] */
    strcpy(out_metadata[OUT_CANOPY_SWQ].varname, "OUT_CANOPY_SWQ");
    strcpy(out_metadata[OUT_CANOPY_SWQ].long_name, "canopy_swe");
    strcpy(out_metadata[OUT_CANOPY_SWQ].standard_name,
           "lwe_thickness_of_canopy_intercepted_swe");
    strcpy(out_metadata[OUT_CANOPY_SWQ].units, "mm");
    strcpy(out_metadata[OUT_CANOPY_SWQ].description,
           "precipitation interception storage in canopy");

    /* surface water or glacier (mm) */
    strcpy(out_metadata[OUT_H2OSFC].varname, "OUT_H2OSFC");
    strcpy(out_metadata[OUT_H2OSFC].long_name, "h2osfc");
    strcpy(out_metadata[OUT_H2OSFC].standard_name, "surface_water_content");
    strcpy(out_metadata[OUT_H2OSFC].units, "mm");
    strcpy(out_metadata[OUT_H2OSFC].description, "surface water or glacier");

    /* fraction of ground covered by surface water or glacier */
    strcpy(out_metadata[OUT_H2O_FRAC].varname, "OUT_H2O_FRAC");
    strcpy(out_metadata[OUT_H2O_FRAC].long_name, "h2o_frac");
    strcpy(out_metadata[OUT_H2O_FRAC].standard_name, "surface_water_fraction");
    strcpy(out_metadata[OUT_H2O_FRAC].units, "1");
    strcpy(out_metadata[OUT_H2O_FRAC].description, "fraction of ground covered by surface water or glacier");

    /* surface water or glacier liquid content (mm) */
    strcpy(out_metadata[OUT_H2OSFC_LIQ].varname, "OUT_H2OSFC_LIQ");
    strcpy(out_metadata[OUT_H2OSFC_LIQ].long_name, "h2osfc_liq");
    strcpy(out_metadata[OUT_H2OSFC_LIQ].standard_name, "surface_water_liquid_content");
    strcpy(out_metadata[OUT_H2OSFC_LIQ].units, "mm");
    strcpy(out_metadata[OUT_H2OSFC_LIQ].description, "surface water or glacier liquid content");

    /* surface water or glacier ice content (mm) */
    strcpy(out_metadata[OUT_H2OSFC_ICE].varname, "OUT_H2OSFC_ICE");
    strcpy(out_metadata[OUT_H2OSFC_ICE].long_name, "h2osfc_ice");
    strcpy(out_metadata[OUT_H2OSFC_ICE].standard_name, "surface_water_ice_content");
    strcpy(out_metadata[OUT_H2OSFC_ICE].units, "mm");
    strcpy(out_metadata[OUT_H2OSFC_ICE].description, "surface water or glacier ice content");

    /* soil matric potential [m] */
    strcpy(out_metadata[OUT_MATRIC].varname, "OUT_MATRIC");
    strcpy(out_metadata[OUT_MATRIC].long_name, "matric");
    strcpy(out_metadata[OUT_MATRIC].standard_name, "soil_matric_potential");
    strcpy(out_metadata[OUT_MATRIC].units, "m");
    strcpy(out_metadata[OUT_MATRIC].description, "soil matric potential");

    /* bulk density of snowfall [kg/m3] */
    strcpy(out_metadata[OUT_NEW_DENSITY].varname, "OUT_NEW_DENSITY");
    strcpy(out_metadata[OUT_NEW_DENSITY].long_name, "new_density");
    strcpy(out_metadata[OUT_NEW_DENSITY].standard_name, "fresh_snow_bulk_density");
    strcpy(out_metadata[OUT_NEW_DENSITY].units, "kg m-3");
    strcpy(out_metadata[OUT_NEW_DENSITY].description, "bulk density of snowfall");

    /* snow age (s) */
    strcpy(out_metadata[OUT_SNOW_AGE].varname, "OUT_SNOW_AGE");
    strcpy(out_metadata[OUT_SNOW_AGE].long_name, "snow_age");
    strcpy(out_metadata[OUT_SNOW_AGE].standard_name, "snow_age");
    strcpy(out_metadata[OUT_SNOW_AGE].units, "1");
    strcpy(out_metadata[OUT_SNOW_AGE].description, "snow age");

    /* fractional area of snow cover [fraction] */
    strcpy(out_metadata[OUT_SNOW_COVER].varname, "OUT_SNOW_COVER");
    strcpy(out_metadata[OUT_SNOW_COVER].long_name, "snow_cover");
    strcpy(out_metadata[OUT_SNOW_COVER].standard_name,
           "snow_cover_area_fraction");
    strcpy(out_metadata[OUT_SNOW_COVER].units, "1");
    strcpy(out_metadata[OUT_SNOW_COVER].description,
           "fractional area of snow cover");

    /* snow density (kg/m^3) */
    strcpy(out_metadata[OUT_SNOW_DENSITY].varname, "OUT_SNOW_DENSITY");
    strcpy(out_metadata[OUT_SNOW_DENSITY].long_name, "snow_density");
    strcpy(out_metadata[OUT_SNOW_DENSITY].standard_name, "snow_bulk_density");
    strcpy(out_metadata[OUT_SNOW_DENSITY].units, "kg m-3");
    strcpy(out_metadata[OUT_SNOW_DENSITY].description, "snow density");

    /* excess liquid water when snow layers combine (mm) */
    strcpy(out_metadata[OUT_SNOW_COMB].varname, "OUT_SNOW_COMB");
    strcpy(out_metadata[OUT_SNOW_COMB].long_name, "snow_comb");
    strcpy(out_metadata[OUT_SNOW_COMB].standard_name, "snow_layer_combine_excess_water");
    strcpy(out_metadata[OUT_SNOW_COMB].units, "mm");
    strcpy(out_metadata[OUT_SNOW_COMB].description, "excess liquid water when snow layers combine");

    /* depth of snow pack [m] */
    strcpy(out_metadata[OUT_SNOW_DEPTH].varname, "OUT_SNOW_DEPTH");
    strcpy(out_metadata[OUT_SNOW_DEPTH].long_name, "snow_depth");
    strcpy(out_metadata[OUT_SNOW_DEPTH].standard_name, "thickness_of_snow");
    strcpy(out_metadata[OUT_SNOW_DEPTH].units, "m");
    strcpy(out_metadata[OUT_SNOW_DEPTH].description, "depth of snow pack");

    /* ice content of the snow pack (mm) */
    strcpy(out_metadata[OUT_SNOW_PACK_ICE].varname, "OUT_SNOW_PACK_ICE");
    strcpy(out_metadata[OUT_SNOW_PACK_ICE].long_name, "snow_pack_ice");
    strcpy(out_metadata[OUT_SNOW_PACK_ICE].standard_name, "snow_pack_ice_content");
    strcpy(out_metadata[OUT_SNOW_PACK_ICE].units, "mm");
    strcpy(out_metadata[OUT_SNOW_PACK_ICE].description, "ice content of the snow pack");

    /* liquid water content of the snow pack (mm) */
    strcpy(out_metadata[OUT_SNOW_PACK_LIQ].varname, "OUT_SNOW_PACK_LIQ");
    strcpy(out_metadata[OUT_SNOW_PACK_LIQ].long_name, "snow_pack_liq");
    strcpy(out_metadata[OUT_SNOW_PACK_LIQ].standard_name, "snow_pack_liquid_water_content");
    strcpy(out_metadata[OUT_SNOW_PACK_LIQ].units, "mm");
    strcpy(out_metadata[OUT_SNOW_PACK_LIQ].description, "liquid water content of the snow pack");

    /* partial volume of snow ice [m3/m3] */
    strcpy(out_metadata[OUT_SNOW_ICEFRAC].varname, "OUT_SNOW_ICEFRAC");
    strcpy(out_metadata[OUT_SNOW_ICEFRAC].long_name, "snow_icefrac");
    strcpy(out_metadata[OUT_SNOW_ICEFRAC].standard_name, "snow_ice_volume_fraction");
    strcpy(out_metadata[OUT_SNOW_ICEFRAC].units, "m3 m-3");
    strcpy(out_metadata[OUT_SNOW_ICEFRAC].description, "partial volume of snow ice");

    /* partial volume of snow liquid water [m3/m3] */
    strcpy(out_metadata[OUT_SNOW_LIQFRAC].varname, "OUT_SNOW_LIQFRAC");
    strcpy(out_metadata[OUT_SNOW_LIQFRAC].long_name, "snow_liqfrac");
    strcpy(out_metadata[OUT_SNOW_LIQFRAC].standard_name, "snow_liquid_water_volume_fraction");
    strcpy(out_metadata[OUT_SNOW_LIQFRAC].units, "m3 m-3");
    strcpy(out_metadata[OUT_SNOW_LIQFRAC].description, "partial volume of snow liquid water");

    /* partial volume of snow liquid water [m3/m3] */
    strcpy(out_metadata[OUT_SNOW_POROSITY].varname, "OUT_SNOW_POROSITY");
    strcpy(out_metadata[OUT_SNOW_POROSITY].long_name, "snow_porosity");
    strcpy(out_metadata[OUT_SNOW_POROSITY].standard_name, "snow_porosity");
    strcpy(out_metadata[OUT_SNOW_POROSITY].units, "m3 m-3");
    strcpy(out_metadata[OUT_SNOW_POROSITY].description, "partial volume of snow liquid water");

    /* effective grain radius [m-6] */
    strcpy(out_metadata[OUT_SNOW_RADIUS].varname, "OUT_SNOW_RADIUS");
    strcpy(out_metadata[OUT_SNOW_RADIUS].long_name, "snow_radius");
    strcpy(out_metadata[OUT_SNOW_RADIUS].standard_name, "snow_grain_effective_radius");
    strcpy(out_metadata[OUT_SNOW_RADIUS].units, "m-6");
    strcpy(out_metadata[OUT_SNOW_RADIUS].description, "effective grain radius");

    /* soil ice content [mm] for each soil layer */
    strcpy(out_metadata[OUT_SOIL_ICE].varname, "OUT_SOIL_ICE");
    strcpy(out_metadata[OUT_SOIL_ICE].long_name, "soil_ice");
    strcpy(out_metadata[OUT_SOIL_ICE].standard_name, "soil_moisture_ice_depth");
    strcpy(out_metadata[OUT_SOIL_ICE].units, "m3/m3");
    strcpy(out_metadata[OUT_SOIL_ICE].description,
           "soil ice content for each soil layer");

    /* soil liquid moisture content [m3/m3] for each soil layer */
    strcpy(out_metadata[OUT_SOIL_LIQ].varname, "OUT_SOIL_LIQ");
    strcpy(out_metadata[OUT_SOIL_LIQ].long_name, "soil_liq");
    strcpy(out_metadata[OUT_SOIL_LIQ].standard_name,
           "soil_moisture_liquid_depth");
    strcpy(out_metadata[OUT_SOIL_LIQ].units, "m3/m3");
    strcpy(out_metadata[OUT_SOIL_LIQ].description,
           "soil liquid moisture content for each soil layer");

    /* soil total moisture content [mm] for each soil layer */
    strcpy(out_metadata[OUT_SOIL_MOIST].varname, "OUT_SOIL_MOIST");
    strcpy(out_metadata[OUT_SOIL_MOIST].long_name, "soil_moist");
    strcpy(out_metadata[OUT_SOIL_MOIST].standard_name, "soil_moisture");
    strcpy(out_metadata[OUT_SOIL_MOIST].units, "m3/m3");
    strcpy(out_metadata[OUT_SOIL_MOIST].description,
           "soil total moisture content");

    /* snow water equivalent in snow pack [mm] */
    strcpy(out_metadata[OUT_SWE].varname, "OUT_SWE");
    strcpy(out_metadata[OUT_SWE].long_name, "swe");
    strcpy(out_metadata[OUT_SWE].standard_name, "lwe_thickness_of_snow");
    strcpy(out_metadata[OUT_SWE].units, "mm");
    strcpy(out_metadata[OUT_SWE].description,
           "snow water equivalent in snow pack");

    /* water storage in aquifer [m] */
    strcpy(out_metadata[OUT_STORAGE_AQF].varname, "OUT_STORAGE_AQF");
    strcpy(out_metadata[OUT_STORAGE_AQF].long_name, "storage_aqf");
    strcpy(out_metadata[OUT_STORAGE_AQF].standard_name, "aquifer_water_storage");
    strcpy(out_metadata[OUT_STORAGE_AQF].units, "m");
    strcpy(out_metadata[OUT_STORAGE_AQF].description, "water storage in aquifer");

    /* snow intercepted on canopy (mm) */
    strcpy(out_metadata[OUT_INT_SNOW].varname, "OUT_INT_SNOW");
    strcpy(out_metadata[OUT_INT_SNOW].long_name, "int_snow");
    strcpy(out_metadata[OUT_INT_SNOW].standard_name, "canopy_intercepted_snow");
    strcpy(out_metadata[OUT_INT_SNOW].units, "mm");
    strcpy(out_metadata[OUT_INT_SNOW].description, "snow intercepted on canopy");

    /* rain intercepted on canopy (mm) */
    strcpy(out_metadata[OUT_INT_RAIN].varname, "OUT_INT_RAIN");
    strcpy(out_metadata[OUT_INT_RAIN].long_name, "int_rain");
    strcpy(out_metadata[OUT_INT_RAIN].standard_name, "canopy_intercepted_rain");
    strcpy(out_metadata[OUT_INT_RAIN].units, "mm");
    strcpy(out_metadata[OUT_INT_RAIN].description, "rain intercepted on canopy");

    /* vegetation water matric potential (m) [sun, shade, xylem, root] */
    strcpy(out_metadata[OUT_VEG_MATRIC].varname, "OUT_VEG_MATRIC");
    strcpy(out_metadata[OUT_VEG_MATRIC].long_name, "veg_matric");
    strcpy(out_metadata[OUT_VEG_MATRIC].standard_name, "vegetation_water_matric_potential");
    strcpy(out_metadata[OUT_VEG_MATRIC].units, "m");
    strcpy(out_metadata[OUT_VEG_MATRIC].description, "vegetation water matric potential (sun, shade, xylem, root)");

    /* total moisture interception storage in canopy [mm] */
    strcpy(out_metadata[OUT_WDEW].varname, "OUT_WDEW");
    strcpy(out_metadata[OUT_WDEW].long_name, "wdew");
    strcpy(out_metadata[OUT_WDEW].standard_name,
           "soil_moisture_storage_depth_in_canopy");
    strcpy(out_metadata[OUT_WDEW].units, "mm");
    strcpy(out_metadata[OUT_WDEW].description,
           "total moisture interception storage in canopy");

    /* snow surface frost [mm] */
    strcpy(out_metadata[OUT_SNOW_FROST].varname, "OUT_SNOW_FROST");
    strcpy(out_metadata[OUT_SNOW_FROST].long_name, "frost_snow");
    strcpy(out_metadata[OUT_SNOW_FROST].standard_name, "thickness_of_frost_snow");
    strcpy(out_metadata[OUT_SNOW_FROST].units, "mm");
    strcpy(out_metadata[OUT_SNOW_FROST].description, "depth of frost snow pack");

    /* snow surface dew [mm] */
    strcpy(out_metadata[OUT_SNOW_DEW].varname, "OUT_SNOW_DEW");
    strcpy(out_metadata[OUT_SNOW_DEW].long_name, "dew_snow");
    strcpy(out_metadata[OUT_SNOW_DEW].standard_name, "thickness_of_dew_snow");
    strcpy(out_metadata[OUT_SNOW_DEW].units, "mm");
    strcpy(out_metadata[OUT_SNOW_DEW].description, "depth of dew on snow surface");

    /* snow surface evaporation [mm] */
    strcpy(out_metadata[OUT_SNOW_EVAP].varname, "OUT_SNOW_EVAP");
    strcpy(out_metadata[OUT_SNOW_EVAP].long_name, "evap_snow");
    strcpy(out_metadata[OUT_SNOW_EVAP].standard_name, "thickness_of_evap_snow");
    strcpy(out_metadata[OUT_SNOW_EVAP].units, "mm");
    strcpy(out_metadata[OUT_SNOW_EVAP].description, "snow surface evaporation amount");

    /* snow surface sublimation [mm] */
    strcpy(out_metadata[OUT_SNOW_SUBLIM].varname, "OUT_SNOW_SUBLIM");
    strcpy(out_metadata[OUT_SNOW_SUBLIM].long_name, "sublim_snow");
    strcpy(out_metadata[OUT_SNOW_SUBLIM].standard_name, "thickness_of_sublim_snow");
    strcpy(out_metadata[OUT_SNOW_SUBLIM].units, "mm");
    strcpy(out_metadata[OUT_SNOW_SUBLIM].description, "snow surface sublimation amount");

    /* snow surface freezing [mm] */
    strcpy(out_metadata[OUT_SNOW_FRZE].varname, "OUT_SNOW_FRZE");
    strcpy(out_metadata[OUT_SNOW_FRZE].long_name, "freeze_snow");
    strcpy(out_metadata[OUT_SNOW_FRZE].standard_name, "thickness_of_freeze_snow");
    strcpy(out_metadata[OUT_SNOW_FRZE].units, "mm");
    strcpy(out_metadata[OUT_SNOW_FRZE].description, "snow surface freezing amount");

    /* water table position [m] (zwt within lowest unsaturated layer) */
    strcpy(out_metadata[OUT_ZWT].varname, "OUT_ZWT");
    strcpy(out_metadata[OUT_ZWT].long_name, "zwt");
    strcpy(out_metadata[OUT_ZWT].standard_name,
           "water_table_position_lowest_layer");
    strcpy(out_metadata[OUT_ZWT].units, "m");
    strcpy(out_metadata[OUT_ZWT].description,
           "water table position (zwt within lowest unsaturated layer)");

    // Water Balance Terms - fluxes
    /* baseflow out of the bottom layer [mm] */
    strcpy(out_metadata[OUT_BASEFLOW].varname, "OUT_BASEFLOW");
    strcpy(out_metadata[OUT_BASEFLOW].long_name, "baseflow");
    strcpy(out_metadata[OUT_BASEFLOW].standard_name, "baseflow_amount");
    strcpy(out_metadata[OUT_BASEFLOW].units, "mm");
    strcpy(out_metadata[OUT_BASEFLOW].description,
           "baseflow out of the bottom layer");

    /* net evaporation from canopy interception [mm] */
    strcpy(out_metadata[OUT_DEW_CANOP].varname, "OUT_DEW_CANOP");
    strcpy(out_metadata[OUT_DEW_CANOP].long_name, "dew_canop");
    strcpy(out_metadata[OUT_DEW_CANOP].standard_name, "canopy_dew");
    strcpy(out_metadata[OUT_DEW_CANOP].units, "mm");
    strcpy(out_metadata[OUT_DEW_CANOP].description, "dew on canopy");

    /* snow depth increasing rate */
    strcpy(out_metadata[OUT_DELDEPTH].varname, "OUT_DELDEPTH");
    strcpy(out_metadata[OUT_DELDEPTH].long_name, "deldepth");
    strcpy(out_metadata[OUT_DELDEPTH].standard_name, "snow_depth_increase_rate");
    strcpy(out_metadata[OUT_DELDEPTH].units, "m s-1");
    strcpy(out_metadata[OUT_DELDEPTH].description, "snow depth increasing rate");

    /* recharge to groundwater storage [m/s] */
    strcpy(out_metadata[OUT_RECHARGE].varname, "OUT_RECHARGE");
    strcpy(out_metadata[OUT_RECHARGE].long_name, "recharge");
    strcpy(out_metadata[OUT_RECHARGE].standard_name, "recharge_groundwater");
    strcpy(out_metadata[OUT_RECHARGE].units, "m/s");
    strcpy(out_metadata[OUT_RECHARGE].description, "recharge to groundwater storage");

    /* total net evaporation [mm] */
    strcpy(out_metadata[OUT_EVAP].varname, "OUT_EVAP");
    strcpy(out_metadata[OUT_EVAP].long_name, "evap");
    strcpy(out_metadata[OUT_EVAP].standard_name, "water_evaporation_flux_net");
    strcpy(out_metadata[OUT_EVAP].units, "mm");
    strcpy(out_metadata[OUT_EVAP].description, "total net evaporation");

    /* net evaporation from bare soil [mm] */
    strcpy(out_metadata[OUT_EVAP_BARE].varname, "OUT_EVAP_BARE");
    strcpy(out_metadata[OUT_EVAP_BARE].long_name, "evap_bare");
    strcpy(out_metadata[OUT_EVAP_BARE].standard_name,
           "water_evaporation_from_soil");
    strcpy(out_metadata[OUT_EVAP_BARE].units, "mm");
    strcpy(out_metadata[OUT_EVAP_BARE].description,
           "net evaporation from bare soil");

    /* net evaporation from canopy interception [mm] */
    strcpy(out_metadata[OUT_EVAP_CANOP].varname, "OUT_EVAP_CANOP");
    strcpy(out_metadata[OUT_EVAP_CANOP].long_name, "evap_canop");
    strcpy(out_metadata[OUT_EVAP_CANOP].standard_name,
           "water_evaporation_from_canopy");
    strcpy(out_metadata[OUT_EVAP_CANOP].units, "mm");
    strcpy(out_metadata[OUT_EVAP_CANOP].description,
           "net evaporation from canopy interception");

    /* glacier snow excess flow */
    strcpy(out_metadata[OUT_GLAC_EXCESS].varname, "OUT_GLAC_EXCESS");
    strcpy(out_metadata[OUT_GLAC_EXCESS].long_name, "glac_excess");
    strcpy(out_metadata[OUT_GLAC_EXCESS].standard_name, "glacier_excess_flow");
    strcpy(out_metadata[OUT_GLAC_EXCESS].units, "mm");
    strcpy(out_metadata[OUT_GLAC_EXCESS].description, "glacier snow excess flow");

    /* glacier mass balance [mm] */
    strcpy(out_metadata[OUT_GLAC_MBAL].varname, "OUT_GLAC_MBAL");
    strcpy(out_metadata[OUT_GLAC_MBAL].long_name, "glac_mbal");
    strcpy(out_metadata[OUT_GLAC_MBAL].standard_name, "glacier_mass_balance");
    strcpy(out_metadata[OUT_GLAC_MBAL].units, "mm");
    strcpy(out_metadata[OUT_GLAC_MBAL].description, "glacier mass balance");

    /* glacier ice accumulation from conversion of firn to ice [mm] */
    strcpy(out_metadata[OUT_GLAC_ACCUM].varname, "OUT_GLAC_ACCUM");
    strcpy(out_metadata[OUT_GLAC_ACCUM].long_name, "glac_accum");
    strcpy(out_metadata[OUT_GLAC_ACCUM].standard_name, "glacier_accumulation");
    strcpy(out_metadata[OUT_GLAC_ACCUM].units, "mm");
    strcpy(out_metadata[OUT_GLAC_ACCUM].description, "glacier ice accumulation from conversion of firn to ice");

    /* glacier ice melt [mm] */
    strcpy(out_metadata[OUT_GLAC_MELT].varname, "OUT_GLAC_MELT");
    strcpy(out_metadata[OUT_GLAC_MELT].long_name, "glac_melt");
    strcpy(out_metadata[OUT_GLAC_MELT].standard_name, "glacier_ice_melt");
    strcpy(out_metadata[OUT_GLAC_MELT].units, "mm");
    strcpy(out_metadata[OUT_GLAC_MELT].description, "glacier ice melt");

    /* glacier water inflow from snow melt, ice melt and rainfall [mm] */
    strcpy(out_metadata[OUT_GLAC_INFLOW].varname, "OUT_GLAC_INFLOW");
    strcpy(out_metadata[OUT_GLAC_INFLOW].long_name, "glac_inflow");
    strcpy(out_metadata[OUT_GLAC_INFLOW].standard_name, "glacier_water_inflow");
    strcpy(out_metadata[OUT_GLAC_INFLOW].units, "mm");
    strcpy(out_metadata[OUT_GLAC_INFLOW].description, "glacier water inflow from snow melt, ice melt and rainfall");

    /* glacier water outflow [mm] */
    strcpy(out_metadata[OUT_GLAC_OUTFLOW].varname, "OUT_GLAC_OUTFLOW");
    strcpy(out_metadata[OUT_GLAC_OUTFLOW].long_name, "glac_outflow");
    strcpy(out_metadata[OUT_GLAC_OUTFLOW].standard_name, "glacier_water_outflow");
    strcpy(out_metadata[OUT_GLAC_OUTFLOW].units, "mm");
    strcpy(out_metadata[OUT_GLAC_OUTFLOW].description, "glacier water outflow");

    /* frost on canopy */
    strcpy(out_metadata[OUT_FROST_CANOP].varname, "OUT_FROST_CANOP");
    strcpy(out_metadata[OUT_FROST_CANOP].long_name, "frost_canop");
    strcpy(out_metadata[OUT_FROST_CANOP].standard_name, "canopy_frost");
    strcpy(out_metadata[OUT_FROST_CANOP].units, "mm");
    strcpy(out_metadata[OUT_FROST_CANOP].description, "frost on canopy");

    /* soil surface dew [mm] */
    strcpy(out_metadata[OUT_SOIL_DEW].varname, "OUT_SOIL_DEW");
    strcpy(out_metadata[OUT_SOIL_DEW].long_name, "dew_soil");
    strcpy(out_metadata[OUT_SOIL_DEW].standard_name, "thickness_of_dew_soil");
    strcpy(out_metadata[OUT_SOIL_DEW].units, "mm");
    strcpy(out_metadata[OUT_SOIL_DEW].description, "depth of dew on soil surface");

    /* soil surface frost [mm] */
    strcpy(out_metadata[OUT_SOIL_FROST].varname, "OUT_SOIL_FROST");
    strcpy(out_metadata[OUT_SOIL_FROST].long_name, "frost_soil");
    strcpy(out_metadata[OUT_SOIL_FROST].standard_name, "thickness_of_frost_soil");
    strcpy(out_metadata[OUT_SOIL_FROST].units, "mm");
    strcpy(out_metadata[OUT_SOIL_FROST].description, "depth of frost on soil surface");

    /* soil evaporation from soil layer [mm] */
    strcpy(out_metadata[OUT_SOIL_EVAP].varname, "OUT_SOIL_EVAP");
    strcpy(out_metadata[OUT_SOIL_EVAP].long_name, "evap_soil");
    strcpy(out_metadata[OUT_SOIL_EVAP].standard_name, "thickness_of_evap_soil");
    strcpy(out_metadata[OUT_SOIL_EVAP].units, "mm");
    strcpy(out_metadata[OUT_SOIL_EVAP].description, "soil evaporation amount from soil layer");

    /* soil surface sublimation [mm] */
    strcpy(out_metadata[OUT_SOIL_SUBLIM].varname, "OUT_SOIL_SUBLIM");
    strcpy(out_metadata[OUT_SOIL_SUBLIM].long_name, "sublim_soil");
    strcpy(out_metadata[OUT_SOIL_SUBLIM].standard_name, "thickness_of_sublim_soil");
    strcpy(out_metadata[OUT_SOIL_SUBLIM].units, "mm");
    strcpy(out_metadata[OUT_SOIL_SUBLIM].description, "soil surface sublimation amount");

    /* net evaporation from canopy interception [mm] */
    strcpy(out_metadata[OUT_NET_EVAP].varname, "OUT_NET_EVAP");
    strcpy(out_metadata[OUT_NET_EVAP].long_name, "evap");
    strcpy(out_metadata[OUT_NET_EVAP].standard_name,
           "water_evaporation_from_canopy");
    strcpy(out_metadata[OUT_NET_EVAP].units, "mm");
    strcpy(out_metadata[OUT_NET_EVAP].description,
           "net evaporation from canopy interception");

    /* outflow of liq water from each snow pack (m/s) */
    strcpy(out_metadata[OUT_PACK_OUTFLOW].varname, "OUT_PACK_OUTFLOW");
    strcpy(out_metadata[OUT_PACK_OUTFLOW].long_name, "pack_outflow");
    strcpy(out_metadata[OUT_PACK_OUTFLOW].standard_name, "snow_pack_outflow");
    strcpy(out_metadata[OUT_PACK_OUTFLOW].units, "m s-1");
    strcpy(out_metadata[OUT_PACK_OUTFLOW].description, "outflow of liquid water from each snow pack");

    /* snow that reaches the ground through the canopy [mm] */
    strcpy(out_metadata[OUT_SNOWTHROUGHFALL].varname, "OUT_SNOWTHROUGHFALL");
    strcpy(out_metadata[OUT_SNOWTHROUGHFALL].long_name, "snow_throughfall");
    strcpy(out_metadata[OUT_SNOWTHROUGHFALL].standard_name, "thickness_of_snow_throughfall");
    strcpy(out_metadata[OUT_SNOWTHROUGHFALL].units, "mm");
    strcpy(out_metadata[OUT_SNOWTHROUGHFALL].description, "snow that reaches the ground through the canopy");

    /* rain that reaches the ground through the canopy [mm] */
    strcpy(out_metadata[OUT_RAINTHROUGHFALL].varname, "OUT_RAINTHROUGHFALL");
    strcpy(out_metadata[OUT_RAINTHROUGHFALL].long_name, "rain_throughfall");
    strcpy(out_metadata[OUT_RAINTHROUGHFALL].standard_name, "thickness_of_rain_throughfall");
    strcpy(out_metadata[OUT_RAINTHROUGHFALL].units, "mm");
    strcpy(out_metadata[OUT_RAINTHROUGHFALL].description, "rain that reaches the ground through the canopy");

    /* rain drip from canopy [mm] */
    strcpy(out_metadata[OUT_RAIN_DRIP].varname, "OUT_RAIN_DRIP");
    strcpy(out_metadata[OUT_RAIN_DRIP].long_name, "rain_drip");
    strcpy(out_metadata[OUT_RAIN_DRIP].standard_name, "thickness_of_rain_drip");
    strcpy(out_metadata[OUT_RAIN_DRIP].units, "mm");
    strcpy(out_metadata[OUT_RAIN_DRIP].description, "rain dripping from canopy");

    /* snow drip from canopy [mm] */
    strcpy(out_metadata[OUT_SNOW_DRIP].varname, "OUT_SNOW_DRIP");
    strcpy(out_metadata[OUT_SNOW_DRIP].long_name, "snow_drip");
    strcpy(out_metadata[OUT_SNOW_DRIP].standard_name, "thickness_of_snow_drip");
    strcpy(out_metadata[OUT_SNOW_DRIP].units, "mm");
    strcpy(out_metadata[OUT_SNOW_DRIP].description, "snow dripping from canopy");

    /* outflow of liquid water from the snowpack bottom [mm] */
    strcpy(out_metadata[OUT_SNOW_OUTFLOW].varname, "OUT_SNOW_OUTFLOW");
    strcpy(out_metadata[OUT_SNOW_OUTFLOW].long_name, "snow_outflow");
    strcpy(out_metadata[OUT_SNOW_OUTFLOW].standard_name, "thickness_of_snow_outflow");
    strcpy(out_metadata[OUT_SNOW_OUTFLOW].units, "mm");
    strcpy(out_metadata[OUT_SNOW_OUTFLOW].description, "outflow of liquid water from the snowpack bottom");

    /* moisture that reaches top of soil column [mm] */
    strcpy(out_metadata[OUT_INFLOW].varname, "OUT_INFLOW");
    strcpy(out_metadata[OUT_INFLOW].long_name, "inflow");
    strcpy(out_metadata[OUT_INFLOW].standard_name, "soil_column_inflow");
    strcpy(out_metadata[OUT_INFLOW].units, "mm");
    strcpy(out_metadata[OUT_INFLOW].description,
           "moisture that reaches top of soil column");

    /* incoming precipitation [mm] */
    strcpy(out_metadata[OUT_PREC].varname, "OUT_PREC");
    strcpy(out_metadata[OUT_PREC].long_name, "prec");
    strcpy(out_metadata[OUT_PREC].standard_name, "precipitation_amount");
    strcpy(out_metadata[OUT_PREC].units, "mm");
    strcpy(out_metadata[OUT_PREC].description, "incoming precipitation");

    /* rainfall [mm] */
    strcpy(out_metadata[OUT_RAINF].varname, "OUT_RAINF");
    strcpy(out_metadata[OUT_RAINF].long_name, "rainf");
    strcpy(out_metadata[OUT_RAINF].standard_name, "rainfall_amount");
    strcpy(out_metadata[OUT_RAINF].units, "mm");
    strcpy(out_metadata[OUT_RAINF].description, "liquid rainfall amount");

    /* surface runoff [mm] */
    strcpy(out_metadata[OUT_RUNOFF].varname, "OUT_RUNOFF");
    strcpy(out_metadata[OUT_RUNOFF].long_name, "runoff");
    strcpy(out_metadata[OUT_RUNOFF].standard_name, "runoff_amount");
    strcpy(out_metadata[OUT_RUNOFF].units, "mm");
    strcpy(out_metadata[OUT_RUNOFF].description, "surface runoff");

    /* snow melt [mm] */
    strcpy(out_metadata[OUT_SNOW_MELT].varname, "OUT_SNOW_MELT");
    strcpy(out_metadata[OUT_SNOW_MELT].long_name, "snow_melt");
    strcpy(out_metadata[OUT_SNOW_MELT].standard_name, "snow_melt_amount");
    strcpy(out_metadata[OUT_SNOW_MELT].units, "mm");
    strcpy(out_metadata[OUT_SNOW_MELT].description, "snow melt");

    /* snowfall [mm] */
    strcpy(out_metadata[OUT_SNOWF].varname, "OUT_SNOWF");
    strcpy(out_metadata[OUT_SNOWF].long_name, "snowf");
    strcpy(out_metadata[OUT_SNOWF].standard_name, "snowfall_lwe_amount");
    strcpy(out_metadata[OUT_SNOWF].units, "mm");
    strcpy(out_metadata[OUT_SNOWF].description, "snowfall");

    /* net sublimation from snow stored in canopy [mm] */
    strcpy(out_metadata[OUT_SUB_CANOP].varname, "OUT_SUB_CANOP");
    strcpy(out_metadata[OUT_SUB_CANOP].long_name, "sub_canop");
    strcpy(out_metadata[OUT_SUB_CANOP].standard_name,
           "sublimation_amount_from_canopy_snow");
    strcpy(out_metadata[OUT_SUB_CANOP].units, "mm");
    strcpy(out_metadata[OUT_SUB_CANOP].description,
           "net sublimation from snow stored in canopy");

    /* net transpiration from vegetation [mm] */
    strcpy(out_metadata[OUT_TRANSP_VEG].varname, "OUT_TRANSP_VEG");
    strcpy(out_metadata[OUT_TRANSP_VEG].long_name, "transp_veg");
    strcpy(out_metadata[OUT_TRANSP_VEG].standard_name, "transpiration_amount");
    strcpy(out_metadata[OUT_TRANSP_VEG].units, "mm");
    strcpy(out_metadata[OUT_TRANSP_VEG].description,
           "net transpiration from vegetation");

    /* transpiration sink term [m/s] */
    strcpy(out_metadata[OUT_TRANSP_SINK].varname, "OUT_TRANSP_SINK");
    strcpy(out_metadata[OUT_TRANSP_SINK].long_name, "transp_sink");
    strcpy(out_metadata[OUT_TRANSP_SINK].standard_name, "transpiration_sink");
    strcpy(out_metadata[OUT_TRANSP_SINK].units, "m s-1");
    strcpy(out_metadata[OUT_TRANSP_SINK].description, "transpiration sink term");

    /* vapor from canopy */
    strcpy(out_metadata[OUT_VAPOR_CANOP].varname, "OUT_VAPOR_CANOP");
    strcpy(out_metadata[OUT_VAPOR_CANOP].long_name, "vapor_canop");
    strcpy(out_metadata[OUT_VAPOR_CANOP].standard_name, "canopy_vapor");
    strcpy(out_metadata[OUT_VAPOR_CANOP].units, "mm");
    strcpy(out_metadata[OUT_VAPOR_CANOP].description, "vapor from canopy");

    /* water budget error [mm] */
    strcpy(out_metadata[OUT_WATER_ERROR].varname, "OUT_WATER_ERROR");
    strcpy(out_metadata[OUT_WATER_ERROR].long_name, "water_error");
    strcpy(out_metadata[OUT_WATER_ERROR].standard_name, "water_budget_error");
    strcpy(out_metadata[OUT_WATER_ERROR].units, "mm");
    strcpy(out_metadata[OUT_WATER_ERROR].description, "water budget error");

    // Energy Balance Terms - state variables
    /* albedo [fraction] */
    strcpy(out_metadata[OUT_ALBEDO].varname, "OUT_ALBEDO");
    strcpy(out_metadata[OUT_ALBEDO].long_name, "albedo");
    strcpy(out_metadata[OUT_ALBEDO].standard_name, "surface_albedo");
    strcpy(out_metadata[OUT_ALBEDO].units, "1");
    strcpy(out_metadata[OUT_ALBEDO].description, "albedo");

    /* bare soil surface temperature [K] */
    strcpy(out_metadata[OUT_BARESOILT].varname, "OUT_BARESOILT");
    strcpy(out_metadata[OUT_BARESOILT].long_name, "baresoilt");
    strcpy(out_metadata[OUT_BARESOILT].standard_name,
           "surface_temperature_bare_soil");
    strcpy(out_metadata[OUT_BARESOILT].units, "K");
    strcpy(out_metadata[OUT_BARESOILT].description,
           "bare soil surface temperature");

    /* depth of freezing fronts [m] */
    strcpy(out_metadata[OUT_FDEPTH].varname, "OUT_FDEPTH");
    strcpy(out_metadata[OUT_FDEPTH].long_name, "fdepth");
    strcpy(out_metadata[OUT_FDEPTH].standard_name, "freezing_front_depth");
    strcpy(out_metadata[OUT_FDEPTH].units, "m");
    strcpy(out_metadata[OUT_FDEPTH].description, "depth of freezing fronts");

    /* surface water temperature or glacier [K] */
    strcpy(out_metadata[OUT_H2OSFC_T].varname, "OUT_H2OSFC_T");
    strcpy(out_metadata[OUT_H2OSFC_T].long_name, "h2osfc_t");
    strcpy(out_metadata[OUT_H2OSFC_T].standard_name, "surface_water_temperature");
    strcpy(out_metadata[OUT_H2OSFC_T].units, "K");
    strcpy(out_metadata[OUT_H2OSFC_T].description, "surface water temperature or glacier");

    /* average surface temperature [K] */
    strcpy(out_metadata[OUT_SURF_TEMP].varname, "OUT_SURF_TEMP");
    strcpy(out_metadata[OUT_SURF_TEMP].long_name, "surf_temp");
    strcpy(out_metadata[OUT_SURF_TEMP].standard_name,
           "surface_average_temperature");
    strcpy(out_metadata[OUT_SURF_TEMP].units, "K");
    strcpy(out_metadata[OUT_SURF_TEMP].description,
           "average surface temperature");

    /* snow albedo [fraction] */
    strcpy(out_metadata[OUT_SALBEDO].varname, "OUT_SALBEDO");
    strcpy(out_metadata[OUT_SALBEDO].long_name, "salbedo");
    strcpy(out_metadata[OUT_SALBEDO].standard_name, "snow_albedo");
    strcpy(out_metadata[OUT_SALBEDO].units, "1");
    strcpy(out_metadata[OUT_SALBEDO].description, "snow albedo");

    /* snow pack temperature [K] */
    strcpy(out_metadata[OUT_SNOW_PACK_TEMP].varname, "OUT_SNOW_PACK_TEMP");
    strcpy(out_metadata[OUT_SNOW_PACK_TEMP].long_name, "snow_pack_temp");
    strcpy(out_metadata[OUT_SNOW_PACK_TEMP].standard_name,
           "snow_pack_temperature");
    strcpy(out_metadata[OUT_SNOW_PACK_TEMP].units, "K");
    strcpy(out_metadata[OUT_SNOW_PACK_TEMP].description,
           "snow pack temperature");

    /* soil temperature [K] for each soil layer */
    strcpy(out_metadata[OUT_SOIL_TEMP].varname, "OUT_SOIL_TEMP");
    strcpy(out_metadata[OUT_SOIL_TEMP].long_name, "soil_temp");
    strcpy(out_metadata[OUT_SOIL_TEMP].standard_name, "soil_temperature");
    strcpy(out_metadata[OUT_SOIL_TEMP].units, "k");
    strcpy(out_metadata[OUT_SOIL_TEMP].description,
           "soil temperature for each soil layer");

    /* frozen soil present flag */
    strcpy(out_metadata[OUT_TRND_FBFLAG].varname, "OUT_TRND_FBFLAG");
    strcpy(out_metadata[OUT_TRND_FBFLAG].long_name, "trnd_fbflag");
    strcpy(out_metadata[OUT_TRND_FBFLAG].standard_name, "frozen_soil_flag");
    strcpy(out_metadata[OUT_TRND_FBFLAG].units, "1");
    strcpy(out_metadata[OUT_TRND_FBFLAG].description, "frozen soil present flag");

    /* Tcanopy flag */
    strcpy(out_metadata[OUT_TCAN_FBFLAG].varname, "OUT_TCAN_FBFLAG");
    strcpy(out_metadata[OUT_TCAN_FBFLAG].long_name, "tcan_fbflag");
    strcpy(out_metadata[OUT_TCAN_FBFLAG].standard_name,
           "canopy_temperature_flag");
    strcpy(out_metadata[OUT_TCAN_FBFLAG].units, "1");
    strcpy(out_metadata[OUT_TCAN_FBFLAG].description,
           "Canopy temperature fallback flag");

    /* depth of thawing fronts [m] */
    strcpy(out_metadata[OUT_TDEPTH].varname, "OUT_TDEPTH");
    strcpy(out_metadata[OUT_TDEPTH].long_name, "tdepth");
    strcpy(out_metadata[OUT_TDEPTH].standard_name, "thawing_front_depth");
    strcpy(out_metadata[OUT_TDEPTH].units, "m");
    strcpy(out_metadata[OUT_TDEPTH].description, "depth of thawing fronts");

    /* average vegetation canopy temperature [C] */
    strcpy(out_metadata[OUT_VEGT].varname, "OUT_VEGT");
    strcpy(out_metadata[OUT_VEGT].long_name, "vegt");
    strcpy(out_metadata[OUT_VEGT].standard_name, "canopy_temperature");
    strcpy(out_metadata[OUT_VEGT].units, "k");
    strcpy(out_metadata[OUT_VEGT].description,
           "average vegetation canopy temperature");

    /* temperature of the canopy [K] */
    strcpy(out_metadata[OUT_VEGTAIR].varname, "OUT_VEGTAIR");
    strcpy(out_metadata[OUT_VEGTAIR].long_name, "vegtaor");
    strcpy(out_metadata[OUT_VEGTAIR].standard_name, "canopy_air_temperature");
    strcpy(out_metadata[OUT_VEGTAIR].units, "K");
    strcpy(out_metadata[OUT_VEGTAIR].description, "temperature of the canopy");

    /* temperature of the stem [K] */
    strcpy(out_metadata[OUT_VEGTSTEM].varname, "OUT_VEGTSTEM");
    strcpy(out_metadata[OUT_VEGTSTEM].long_name, "vegtstem");
    strcpy(out_metadata[OUT_VEGTSTEM].standard_name, "stem_temperature");
    strcpy(out_metadata[OUT_VEGTSTEM].units, "K");
    strcpy(out_metadata[OUT_VEGTSTEM].description, "temperature of the stem");

    // Energy Balance Terms - fluxes
    /* advected energy [W m-2] */
    strcpy(out_metadata[OUT_ADVECTION].varname, "OUT_ADVECTION");
    strcpy(out_metadata[OUT_ADVECTION].long_name, "advection");
    strcpy(out_metadata[OUT_ADVECTION].standard_name, "advected_energy");
    strcpy(out_metadata[OUT_ADVECTION].units, "W m-2");
    strcpy(out_metadata[OUT_ADVECTION].description, "advected energy ");

    /* advective flux from understory vegetation (Wm-2) */
    strcpy(out_metadata[OUT_ADVECTSUB].varname, "OUT_ADVECTSUB");
    strcpy(out_metadata[OUT_ADVECTSUB].long_name, "advectsub");
    strcpy(out_metadata[OUT_ADVECTSUB].standard_name, "understory_advective_flux");
    strcpy(out_metadata[OUT_ADVECTSUB].units, "W m-2");
    strcpy(out_metadata[OUT_ADVECTSUB].description, "advective flux from understory vegetation");

    /* advective flux from bare ground (Wm-2) */
    strcpy(out_metadata[OUT_ADVECTGRND].varname, "OUT_ADVECTGRND");
    strcpy(out_metadata[OUT_ADVECTGRND].long_name, "advectgrnd");
    strcpy(out_metadata[OUT_ADVECTGRND].standard_name, "ground_advective_flux");
    strcpy(out_metadata[OUT_ADVECTGRND].units, "W m-2");
    strcpy(out_metadata[OUT_ADVECTGRND].description, "advective flux from bare ground");

    /* advective flux from overstory vegetation (Wm-2) */
    strcpy(out_metadata[OUT_ADVECTOVER].varname, "OUT_ADVECTOVER");
    strcpy(out_metadata[OUT_ADVECTOVER].long_name, "advectover");
    strcpy(out_metadata[OUT_ADVECTOVER].standard_name, "overstory_advective_flux");
    strcpy(out_metadata[OUT_ADVECTOVER].units, "W m-2");
    strcpy(out_metadata[OUT_ADVECTOVER].description, "advective flux from overstory vegetation");

    /* energy budget error [W m-2] */
    strcpy(out_metadata[OUT_ENERGY_ERROR].varname, "OUT_ENERGY_ERROR");
    strcpy(out_metadata[OUT_ENERGY_ERROR].long_name, "energy_error");
    strcpy(out_metadata[OUT_ENERGY_ERROR].standard_name, "energy_budget_error");
    strcpy(out_metadata[OUT_ENERGY_ERROR].units, "W m-2");
    strcpy(out_metadata[OUT_ENERGY_ERROR].description, "energy budget error");

    /* net heat flux into ground [W m-2] */
    strcpy(out_metadata[OUT_GRND_FLUX].varname, "OUT_GRND_FLUX");
    strcpy(out_metadata[OUT_GRND_FLUX].long_name, "grnd_flux");
    strcpy(out_metadata[OUT_GRND_FLUX].standard_name,
           "downward_heat_flux_in_soil");
    strcpy(out_metadata[OUT_GRND_FLUX].units, "W m-2");
    strcpy(out_metadata[OUT_GRND_FLUX].description,
           "net heat flux into ground");

    /* net upward latent heat flux [W m-2] */
    strcpy(out_metadata[OUT_LATENT].varname, "OUT_LATENT");
    strcpy(out_metadata[OUT_LATENT].long_name, "latent");
    strcpy(out_metadata[OUT_LATENT].standard_name,
           "surface_upward_latent_heat_flux");
    strcpy(out_metadata[OUT_LATENT].units, "W m-2");
    strcpy(out_metadata[OUT_LATENT].description, "net upward latent heat flux");

    /* net downward longwave flux [W m-2] */
    strcpy(out_metadata[OUT_LWNET].varname, "OUT_LWNET");
    strcpy(out_metadata[OUT_LWNET].long_name, "lwnet");
    strcpy(out_metadata[OUT_LWNET].standard_name,
           "net_downward_longwave_flux_at_surface");
    strcpy(out_metadata[OUT_LWNET].units, "W m-2");
    strcpy(out_metadata[OUT_LWNET].description, "net downward longwave flux");

    /* net downward shortwave flux [W m-2] */
    strcpy(out_metadata[OUT_SWNET].varname, "OUT_SWNET");
    strcpy(out_metadata[OUT_SWNET].long_name, "swnet");
    strcpy(out_metadata[OUT_SWNET].standard_name,
           "net_downward_shortwave_flux_at_surface");
    strcpy(out_metadata[OUT_SWNET].units, "W m-2");
    strcpy(out_metadata[OUT_SWNET].description, "net downward shortwave flux");

    /* net upward sensible heat flux [W m-2] */
    strcpy(out_metadata[OUT_SENSIBLE].varname, "OUT_SENSIBLE");
    strcpy(out_metadata[OUT_SENSIBLE].long_name, "sensible");
    strcpy(out_metadata[OUT_SENSIBLE].standard_name,
           "surface_upward_net_sensible_heat_flux");
    strcpy(out_metadata[OUT_SENSIBLE].units, "W m-2");
    strcpy(out_metadata[OUT_SENSIBLE].description,
           "net upward sensible heat flux");

    /* total absorbed solar radiation by snow for each layer [W/m2] */
    strcpy(out_metadata[OUT_SW_SNOW].varname, "OUT_SW_SNOW");
    strcpy(out_metadata[OUT_SW_SNOW].long_name, "sw_snow");
    strcpy(out_metadata[OUT_SW_SNOW].standard_name, 
           "snow_absorbed_solar_radiation");
    strcpy(out_metadata[OUT_SW_SNOW].units, "W m-2");
    strcpy(out_metadata[OUT_SW_SNOW].description, 
           "total absorbed solar radiation by snow for each layer");

    /* emitted longwave flux from ground (Wm-2) */
    strcpy(out_metadata[OUT_EMISS_LWGRND].varname, "OUT_EMISS_LWGRND");
    strcpy(out_metadata[OUT_EMISS_LWGRND].long_name, "emiss_lwgrnd");
    strcpy(out_metadata[OUT_EMISS_LWGRND].standard_name, 
           "ground_emitted_longwave_flux");
    strcpy(out_metadata[OUT_EMISS_LWGRND].units, "W m-2");
    strcpy(out_metadata[OUT_EMISS_LWGRND].description, 
           "emitted longwave flux from ground");

    /* emitted longwave flux from understory (Wm-2) */
    strcpy(out_metadata[OUT_EMISS_LWSUB].varname, "OUT_EMISS_LWSUB");
    strcpy(out_metadata[OUT_EMISS_LWSUB].long_name, "emiss_lwsub");
    strcpy(out_metadata[OUT_EMISS_LWSUB].standard_name, 
           "understory_emitted_longwave_flux");
    strcpy(out_metadata[OUT_EMISS_LWSUB].units, "W m-2");
    strcpy(out_metadata[OUT_EMISS_LWSUB].description, 
           "emitted longwave flux from understory");

    /* emitted longwave flux from surface (Wm-2) */
    strcpy(out_metadata[OUT_EMISS_LWSURF].varname, "OUT_EMISS_LWSURF");
    strcpy(out_metadata[OUT_EMISS_LWSURF].long_name, "emiss_lwsurf");
    strcpy(out_metadata[OUT_EMISS_LWSURF].standard_name, 
          "surface_emitted_longwave_flux");
    strcpy(out_metadata[OUT_EMISS_LWSURF].units, "W m-2");
    strcpy(out_metadata[OUT_EMISS_LWSURF].description, 
          "emitted longwave flux from surface");

    /* reflected shortwave flux from ground (Wm-2) */
    strcpy(out_metadata[OUT_REFL_SWGRND].varname, "OUT_REFL_SWGRND");
    strcpy(out_metadata[OUT_REFL_SWGRND].long_name, "refl_swgrnd");
    strcpy(out_metadata[OUT_REFL_SWGRND].standard_name, 
           "ground_reflected_shortwave_flux");
    strcpy(out_metadata[OUT_REFL_SWGRND].units, "W m-2");
    strcpy(out_metadata[OUT_REFL_SWGRND].description, 
           "reflected shortwave flux from ground");

    /* reflected shortwave flux from understory (Wm-2) */
    strcpy(out_metadata[OUT_REFL_SWSUB].varname, "OUT_REFL_SWSUB");
    strcpy(out_metadata[OUT_REFL_SWSUB].long_name, "refl_swsub");
    strcpy(out_metadata[OUT_REFL_SWSUB].standard_name, 
          "understory_reflected_shortwave_flux");
    strcpy(out_metadata[OUT_REFL_SWSUB].units, "W m-2");
    strcpy(out_metadata[OUT_REFL_SWSUB].description, 
          "reflected shortwave flux from understory");

    /* reflected shortwave flux from surface (Wm-2) */
    strcpy(out_metadata[OUT_REFL_SWSURF].varname, "OUT_REFL_SWSURF");
    strcpy(out_metadata[OUT_REFL_SWSURF].long_name, "refl_swsurf");
    strcpy(out_metadata[OUT_REFL_SWSURF].standard_name, 
          "surface_reflected_shortwave_flux");
    strcpy(out_metadata[OUT_REFL_SWSURF].units, "W m-2");
    strcpy(out_metadata[OUT_REFL_SWSURF].description, 
          "reflected shortwave flux from surface");

    // Miscellaneous Terms
    /* ground surface resistance to evaporation [m/s] */
    strcpy(out_metadata[OUT_RA_EVAP].varname, "OUT_RA_EVAP");
    strcpy(out_metadata[OUT_RA_EVAP].long_name, "ra_evap");
    strcpy(out_metadata[OUT_RA_EVAP].standard_name, 
          "ground_surface_resistance_to_evaporation");
    strcpy(out_metadata[OUT_RA_EVAP].units, "m s-1");
    strcpy(out_metadata[OUT_RA_EVAP].description, 
          "ground surface resistance to evaporation");

    /* bare ground surface resistance [m/s] */
    strcpy(out_metadata[OUT_RA_GRND].varname, "OUT_RA_GRND");
    strcpy(out_metadata[OUT_RA_GRND].long_name, "ra_grnd");
    strcpy(out_metadata[OUT_RA_GRND].standard_name, "bare_ground_surface_resistance");
    strcpy(out_metadata[OUT_RA_GRND].units, "m s-1");
    strcpy(out_metadata[OUT_RA_GRND].description, "bare ground surface resistance");

    /* canopy leaf resistance to transpiration [s/m] */
    strcpy(out_metadata[OUT_RA_LEAF].varname, "OUT_RA_LEAF");
    strcpy(out_metadata[OUT_RA_LEAF].long_name, "ra_leaf");
    strcpy(out_metadata[OUT_RA_LEAF].standard_name, "canopy_leaf_resistance_to_transpiration");
    strcpy(out_metadata[OUT_RA_LEAF].units, "s m-1");
    strcpy(out_metadata[OUT_RA_LEAF].description, "canopy leaf resistance to transpiration");

    /* canopy aerodynamic resistance [s/m] */
    strcpy(out_metadata[OUT_RA_OVER].varname, "OUT_RA_OVER");
    strcpy(out_metadata[OUT_RA_OVER].long_name, "ra_over");
    strcpy(out_metadata[OUT_RA_OVER].standard_name, "canopy_aerodynamic_resistance");
    strcpy(out_metadata[OUT_RA_OVER].units, "s m-1");
    strcpy(out_metadata[OUT_RA_OVER].description, "canopy aerodynamic resistance");

    /* canopy ground surface resistance [s/m] */
    strcpy(out_metadata[OUT_RA_SUB].varname, "OUT_RA_SUB");
    strcpy(out_metadata[OUT_RA_SUB].long_name, "ra_sub");
    strcpy(out_metadata[OUT_RA_SUB].standard_name, "canopy_ground_surface_resistance");
    strcpy(out_metadata[OUT_RA_SUB].units, "s m-1");
    strcpy(out_metadata[OUT_RA_SUB].description, "canopy ground surface resistance");

    /* veg stem aerodynamic resistance [s/m] */
    strcpy(out_metadata[OUT_RA_STEM].varname, "OUT_RA_STEM");
    strcpy(out_metadata[OUT_RA_STEM].long_name, "ra_stem");
    strcpy(out_metadata[OUT_RA_STEM].standard_name, "stem_aerodynamic_resistance");
    strcpy(out_metadata[OUT_RA_STEM].units, "s m-1");
    strcpy(out_metadata[OUT_RA_STEM].description, "veg stem aerodynamic resistance");

    /* air temperature [k] */
    strcpy(out_metadata[OUT_AIR_TEMP].varname, "OUT_AIR_TEMP");
    strcpy(out_metadata[OUT_AIR_TEMP].long_name, "air_temp");
    strcpy(out_metadata[OUT_AIR_TEMP].standard_name, "air_temperature");
    strcpy(out_metadata[OUT_AIR_TEMP].units, "k");
    strcpy(out_metadata[OUT_AIR_TEMP].description, "air temperature");

    /* atmospheric CO2 concentration [ppm] */
    strcpy(out_metadata[OUT_CATM].varname, "OUT_CATM");
    strcpy(out_metadata[OUT_CATM].long_name, "catm");
    strcpy(out_metadata[OUT_CATM].standard_name,
           "concentration_of_carbon_dioxide_in_air");
    strcpy(out_metadata[OUT_CATM].units, "ppm");
    strcpy(out_metadata[OUT_CATM].description, "atmospheric CO2 concentration");

    /* near-surface atmospheric density [kg m-3] */
    strcpy(out_metadata[OUT_DENSITY].varname, "OUT_DENSITY");
    strcpy(out_metadata[OUT_DENSITY].long_name, "density");
    strcpy(out_metadata[OUT_DENSITY].standard_name, "air_density");
    strcpy(out_metadata[OUT_DENSITY].units, "kg m-3");
    strcpy(out_metadata[OUT_DENSITY].description,
           "near-surface atmospheric density");

    /* fractional area covered by plant canopy [fraction] */
    strcpy(out_metadata[OUT_FCANOPY].varname, "OUT_FCANOPY");
    strcpy(out_metadata[OUT_FCANOPY].long_name, "fcanopy");
    strcpy(out_metadata[OUT_FCANOPY].standard_name,
           "canopy_cover_area_fraction");
    strcpy(out_metadata[OUT_FCANOPY].units, "1");
    strcpy(out_metadata[OUT_FCANOPY].description,
           "fractional area covered by plant canopy");

    /* fraction of incoming shortwave that is direct [fraction] */
    strcpy(out_metadata[OUT_FDIR].varname, "OUT_FDIR");
    strcpy(out_metadata[OUT_FDIR].long_name, "fdir");
    strcpy(out_metadata[OUT_FDIR].standard_name,
           "fraction_of_incoming_shorwave_radiation_that_is_direct");
    strcpy(out_metadata[OUT_FDIR].units, "1");
    strcpy(out_metadata[OUT_FDIR].description,
           "fraction of incoming shortwave that is direct");

    /* Cosine of solar zenith angle [fraction] */
    strcpy(out_metadata[OUT_COSZEN].varname, "OUT_COSZEN");
    strcpy(out_metadata[OUT_COSZEN].long_name, "coszen");
    strcpy(out_metadata[OUT_COSZEN].standard_name,
           "Cosine of solar zenith angle");
    strcpy(out_metadata[OUT_COSZEN].units, "1");
    strcpy(out_metadata[OUT_COSZEN].description,
           "Cosine of solar zenith angle");

    /* leaf area index [1] */
    strcpy(out_metadata[OUT_LAI].varname, "OUT_LAI");
    strcpy(out_metadata[OUT_LAI].long_name, "lai");
    strcpy(out_metadata[OUT_LAI].standard_name, "leaf_area_index");
    strcpy(out_metadata[OUT_LAI].units, "1");
    strcpy(out_metadata[OUT_LAI].description, "leaf area index");

    /* incoming longwave [W m-2] */
    strcpy(out_metadata[OUT_LWDOWN].varname, "OUT_LWDOWN");
    strcpy(out_metadata[OUT_LWDOWN].long_name, "lwdown");
    strcpy(out_metadata[OUT_LWDOWN].standard_name,
           "downwelling_longwave_flux_in_air");
    strcpy(out_metadata[OUT_LWDOWN].units, "W m-2");
    strcpy(out_metadata[OUT_LWDOWN].description, "incoming longwave");

    /* incoming photosynthetically active radiation [W m-2] */
    strcpy(out_metadata[OUT_PAR].varname, "OUT_PAR");
    strcpy(out_metadata[OUT_PAR].long_name, "par");
    strcpy(out_metadata[OUT_PAR].standard_name,
           "surface_downwelling_photosynthetic_radiative_flux_in_air");
    strcpy(out_metadata[OUT_PAR].units, "W m-2");
    strcpy(out_metadata[OUT_PAR].description,
           "incoming photosynthetically active radiation");

    /* near surface atmospheric pressure [kPa] */
    strcpy(out_metadata[OUT_PRESSURE].varname, "OUT_PRESSURE");
    strcpy(out_metadata[OUT_PRESSURE].long_name, "pressure");
    strcpy(out_metadata[OUT_PRESSURE].standard_name, "surface_air_pressure");
    strcpy(out_metadata[OUT_PRESSURE].units, "kPa");
    strcpy(out_metadata[OUT_PRESSURE].description,
           "near surface atmospheric pressure");

    /* specific humidity [kg/kg] */
    strcpy(out_metadata[OUT_QAIR].varname, "OUT_QAIR");
    strcpy(out_metadata[OUT_QAIR].long_name, "qair");
    strcpy(out_metadata[OUT_QAIR].standard_name, "specific_humidity");
    strcpy(out_metadata[OUT_QAIR].units, "1");
    strcpy(out_metadata[OUT_QAIR].description, "specific humidity");

    /* relative humidity [fraction]*/
    strcpy(out_metadata[OUT_REL_HUMID].varname, "OUT_REL_HUMID");
    strcpy(out_metadata[OUT_REL_HUMID].long_name, "rel_humid");
    strcpy(out_metadata[OUT_REL_HUMID].standard_name, "relative_humidity");
    strcpy(out_metadata[OUT_REL_HUMID].units, "1");
    strcpy(out_metadata[OUT_REL_HUMID].description, "relative humidity");

    /* incoming shortwave [W m-2] */
    strcpy(out_metadata[OUT_SWDOWN].varname, "OUT_SWDOWN");
    strcpy(out_metadata[OUT_SWDOWN].long_name, "swdown");
    strcpy(out_metadata[OUT_SWDOWN].standard_name, "incoming shortwave");
    strcpy(out_metadata[OUT_SWDOWN].units, "W m-2");
    strcpy(out_metadata[OUT_SWDOWN].description, "incoming shortwave");

    /* near surface vapor pressure [kPa] */
    strcpy(out_metadata[OUT_VP].varname, "OUT_VP");
    strcpy(out_metadata[OUT_VP].long_name, "vp");
    strcpy(out_metadata[OUT_VP].standard_name, "water_vapor_pressure");
    strcpy(out_metadata[OUT_VP].units, "kPa");
    strcpy(out_metadata[OUT_VP].description, "near surface vapor pressure");

    /* near surface vapor pressure deficit [kPa] */
    strcpy(out_metadata[OUT_VPD].varname, "OUT_VPD");
    strcpy(out_metadata[OUT_VPD].long_name, "vpd");
    strcpy(out_metadata[OUT_VPD].standard_name,
           "water_vapor_saturation_deficit");
    strcpy(out_metadata[OUT_VPD].units, "kPa");
    strcpy(out_metadata[OUT_VPD].description,
           "near surface vapor pressure deficit");

    /* near surface wind speed [m/s] */
    strcpy(out_metadata[OUT_WIND].varname, "OUT_WIND");
    strcpy(out_metadata[OUT_WIND].long_name, "wind");
    strcpy(out_metadata[OUT_WIND].standard_name, "wind_speed");
    strcpy(out_metadata[OUT_WIND].units, "m s-1");
    strcpy(out_metadata[OUT_WIND].description, "near surface wind speed");                                                                                             

    /* Wall time spent inside vic_run [seconds] */
    strcpy(out_metadata[OUT_TIME_VICRUN_WALL].varname, "OUT_TIME_VICRUN_WALL");
    strcpy(out_metadata[OUT_TIME_VICRUN_WALL].long_name, "time_vicrun_wall");
    strcpy(out_metadata[OUT_TIME_VICRUN_WALL].standard_name,
           "vic_run_wall_time");
    strcpy(out_metadata[OUT_TIME_VICRUN_WALL].units, "seconds");
    strcpy(out_metadata[OUT_TIME_VICRUN_WALL].description,
           "Wall time spent inside vic_run");

    /* CPU time spent inside vic_run [seconds] */
    strcpy(out_metadata[OUT_TIME_VICRUN_CPU].varname, "OUT_TIME_VICRUN_CPU");
    strcpy(out_metadata[OUT_TIME_VICRUN_CPU].long_name, "time_vicrun_cpu");
    strcpy(out_metadata[OUT_TIME_VICRUN_CPU].standard_name, "vic_run_cpu_time");
    strcpy(out_metadata[OUT_TIME_VICRUN_CPU].units, "seconds");
    strcpy(out_metadata[OUT_TIME_VICRUN_CPU].description,
           "CPU time spent inside vic_run");

    out_metadata[OUT_SOIL_ICE].nelem = MAX_SOILS;
    out_metadata[OUT_SOIL_ICE].elem_type = OUT_ELEM_SOIL;
    out_metadata[OUT_SOIL_LIQ].nelem = MAX_SOILS;
    out_metadata[OUT_SOIL_LIQ].elem_type = OUT_ELEM_SOIL;
    out_metadata[OUT_SOIL_MOIST].nelem = MAX_SOILS;
    out_metadata[OUT_SOIL_MOIST].elem_type = OUT_ELEM_SOIL;
    out_metadata[OUT_SOIL_TEMP].nelem = MAX_SOILS;
    out_metadata[OUT_SOIL_TEMP].elem_type = OUT_ELEM_SOIL;
    out_metadata[OUT_SOIL_POROSITY].nelem = MAX_SOILS;
    out_metadata[OUT_SOIL_POROSITY].elem_type = OUT_ELEM_SOIL;
    out_metadata[OUT_TRANSP_SINK].nelem = MAX_SOILS;
    out_metadata[OUT_TRANSP_SINK].elem_type = OUT_ELEM_SOIL;
    out_metadata[OUT_MATRIC].nelem = MAX_SOILS;
    out_metadata[OUT_MATRIC].elem_type = OUT_ELEM_SOIL;
    out_metadata[OUT_SNOW_PACK_ICE].nelem = MAX_SNOWS;
    out_metadata[OUT_SNOW_PACK_ICE].elem_type = OUT_ELEM_SNOW;
    out_metadata[OUT_SNOW_PACK_LIQ].nelem = MAX_SNOWS;
    out_metadata[OUT_SNOW_PACK_LIQ].elem_type = OUT_ELEM_SNOW;
    out_metadata[OUT_SNOW_PACK_TEMP].nelem = MAX_SNOWS;
    out_metadata[OUT_SNOW_PACK_TEMP].elem_type = OUT_ELEM_SNOW;
    out_metadata[OUT_SNOW_ICEFRAC].nelem = MAX_SNOWS;
    out_metadata[OUT_SNOW_ICEFRAC].elem_type = OUT_ELEM_SNOW;
    out_metadata[OUT_SNOW_LIQFRAC].nelem = MAX_SNOWS;
    out_metadata[OUT_SNOW_LIQFRAC].elem_type = OUT_ELEM_SNOW;
    out_metadata[OUT_SNOW_POROSITY].nelem = MAX_SNOWS;
    out_metadata[OUT_SNOW_POROSITY].elem_type = OUT_ELEM_SNOW; 
    out_metadata[OUT_SNOW_RADIUS].nelem = MAX_SNOWS;
    out_metadata[OUT_SNOW_RADIUS].elem_type = OUT_ELEM_SNOW;
    out_metadata[OUT_PACK_OUTFLOW].nelem = MAX_SNOWS;
    out_metadata[OUT_PACK_OUTFLOW].elem_type = OUT_ELEM_SNOW;
    out_metadata[OUT_SNOW_MELT].nelem = MAX_SNOWS;
    out_metadata[OUT_SNOW_MELT].elem_type = OUT_ELEM_SNOW;
    out_metadata[OUT_SNOW_FRZE].nelem = MAX_SNOWS;
    out_metadata[OUT_SNOW_FRZE].elem_type = OUT_ELEM_SNOW;
    out_metadata[OUT_VEG_MATRIC].nelem = 4;
    out_metadata[OUT_VEG_MATRIC].elem_type = OUT_ELEM_VEGMAT;
    out_metadata[OUT_RA_OVER].nelem = 3;
    out_metadata[OUT_RA_OVER].elem_type = OUT_ELEM_TURBUL;
    out_metadata[OUT_RA_SUB].nelem = 3;
    out_metadata[OUT_RA_SUB].elem_type = OUT_ELEM_TURBUL;
    out_metadata[OUT_RA_GRND].nelem = 3;
    out_metadata[OUT_RA_GRND].elem_type = OUT_ELEM_TURBUL;
}
