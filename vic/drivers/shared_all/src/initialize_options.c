/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine initalizes all options before they are called by
 * the model.
 *****************************************************************************/

#include "vic_driver_shared_all.h"

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
    options.ACTIVE_LAYER = true;
    options.CARBON = false;
    options.CONTINUEONERROR = true;
    options.CORRPREC = false;
    options.FROZEN_SOIL = true;
    options.NOFLUX = false;
    options.BIOMASST = true;
    options.SNOW_DENSITY = DENS_BRAS;
    options.SNOW_AGING = BATS;
    options.CANOPY_INTERCEP = NOAH;
    options.TFALLBACK = true;
    options.SWRC = SWRC_VAN_GENUCHTEN;
    // Model dimensions
    options.Nlayer = 3;
    options.Nswband = 2;
    options.SNOW_BAND = 1;
    options.GLACIER_ID = 0;
    // input options
    options.AERO_RESIST = AR_MEIER;
    options.FCAN_SRC = FROM_DEFAULT;
    options.GRID_DECIMAL = 2;
    options.LAI_SRC = FROM_VEGLIB;
    options.PARAM_FROM_SOIL = true;
    options.ROUT = false;
    options.VEGLIB_FCAN = false;
    options.VEGPARAM_FCAN = false;
    options.VEGPARAM_LAI = false;
    // state options
    options.STATE_FORMAT = UNSET_FILE_FORMAT;
    options.INIT_STATE = false;
    options.SAVE_STATE = false;
    // output options
    options.Noutstreams = 2;
    options.PRT_HEADER = false;
}
