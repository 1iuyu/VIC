/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine initalizes all options before they are called by
 * the model.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    Initialize all options before they are called by the
 *           model.
 *****************************************************************************/
void
initialize_options()
{
    extern option_struct options;

    /** Initialize model option flags **/

    // simulation modes
    options.AboveTreelineVeg = -1;
    options.AERO_RESIST_CANSNOW = AR_406_FULL;
    options.BLOWING = false;
    options.BLOWING_VAR_THRESHOLD = true;
    options.BLOWING_CALC_PROB = true;
    options.BLOWING_SIMPLE = false;
    options.BLOWING_FETCH = true;
    options.BLOWING_SPATIAL_WIND = true;
    options.CARBON = false;
    options.CLOSE_ENERGY = false;
    options.COMPUTE_TREELINE = false;
    options.CONTINUEONERROR = true;
    options.CORRPREC = false;
    options.EQUAL_AREA = false;
    options.EXP_TRANS = true;
    options.FROZEN_SOIL = false;
    options.FULL_ENERGY = false;
    options.GRND_FLUX_TYPE = GF_410;
    options.IMPLICIT = true;
    options.LAKES = false;
    options.LAKE_PROFILE = false;
    options.NOFLUX = false;
    options.QUICK_FLUX = true;
    options.QUICK_SOLVE = false;
    options.RC_MODE = RC_JARVIS;
    options.SHARE_LAYER_MOIST = true;
    options.SNOW_DENSITY = DENS_BRAS;
    options.SPATIAL_FROST = false;
    options.SPATIAL_SNOW = false;
    options.TFALLBACK = true;
    options.GLACIER_DYNAMICS = true;
    // Model dimensions
    options.Ncanopy = 3;
    options.Nfrost = 1;
    options.Nlakebasnode = MAX_LAKE_BASIN_NODES;
    options.Nlakenode = MAX_LAKE_NODES;
    options.Nlayer = 3;
    options.Nnode = 3;
    options.ROOT_ZONES = 0;
    options.SNOW_BAND = 1;
    // input options
    options.ALB_SRC = FROM_VEGLIB;
    options.BASEFLOW = ARNO;
    options.FCAN_SRC = FROM_DEFAULT;
    options.GRID_DECIMAL = 2;
    options.JULY_TAVG_SUPPLIED = false;
    options.LAI_SRC = FROM_VEGLIB;
    options.ORGANIC_FRACT = false;
    options.BULK_DENSITY_COMB = false;
    options.MAX_SNOW_ALBEDO = false;
    options.VEGLIB_FCAN = false;
    options.VEGLIB_PHOTO = false;
    options.VEGPARAM_ALB = false;
    options.VEGPARAM_FCAN = false;
    options.VEGPARAM_LAI = false;
    // state options
    options.STATE_FORMAT = UNSET_FILE_FORMAT;
    options.STATENAME_CESM = false;
    options.INIT_STATE = false;
    options.SAVE_STATE = false;
    // output options
    options.Noutstreams = 2;
    options.PRT_HEADER = false;
}
