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
    options.BLOWING = false;
    options.BLOWING_VAR_THRESHOLD = true;
    options.BLOWING_CALC_PROB = true;
    options.BLOWING_SIMPLE = false;
    options.BLOWING_FETCH = true;
    options.BLOWING_SPATIAL_WIND = true;
    options.CARBON = false;
    options.CONTINUEONERROR = true;
    options.CORRPREC = false;
    options.EQUAL_AREA = false;
    options.FROZEN_SOIL = false;
    options.NOFLUX = false;
    options.RC_MODE = RC_JARVIS;
    options.SNOW_DENSITY = DENS_BRAS;
    options.SNOW_AGING = BATS;
    options.TFALLBACK = true;
    
    // Model dimensions
    options.Ncanopy = 3;
    options.Nlayer = 3;
    options.Nnode = 3;
    options.Nswband = 2;
    options.SNOW_BAND = 1;
    options.GLACIER_ID = 0;
    // input options
    options.ALB_SRC = FROM_VEGLIB;
    options.BASEFLOW = ARNO;
    options.SOIL_TRANSP = CLM;
    options.FCAN_SRC = FROM_DEFAULT;
    options.GRID_DECIMAL = 2;
    options.LAI_SRC = FROM_VEGLIB;
    options.ORGANIC_FRACT = false;
    options.BULK_DENSITY_COMB = false;
    options.ROUT_PARAM = false;
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
    options.ROUT_PARAM = false;
}
