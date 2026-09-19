/******************************************************************************
 * @section DESCRIPTION
 *
 * Allocate and free memory for the force data struct
 *****************************************************************************/

#include "vic_driver_shared_image.h"

/******************************************************************************
 * @brief    Allocate memory for the atmos data structure.
 *****************************************************************************/
void
alloc_atmos(atmos_data_struct *atmos)
{
    extern option_struct options;

    atmos->air_temp = calloc(NR + 1, sizeof(*(atmos->air_temp)));
    check_alloc_status(atmos->air_temp, "Memory allocation error.");

    atmos->longwave = calloc(NR + 1, sizeof(*(atmos->longwave)));
    check_alloc_status(atmos->longwave, "Memory allocation error.");

    atmos->prec = calloc(NR + 1, sizeof(*(atmos->prec)));
    check_alloc_status(atmos->prec, "Memory allocation error.");

    atmos->pressure = calloc(NR + 1, sizeof(*(atmos->pressure)));
    check_alloc_status(atmos->pressure, "Memory allocation error.");

    atmos->shortwave = calloc(NR + 1, sizeof(*(atmos->shortwave)));
    check_alloc_status(atmos->shortwave, "Memory allocation error.");

    atmos->snowf = calloc(NR + 1, sizeof(*(atmos->snowf)));
    check_alloc_status(atmos->snowf, "Memory allocation error.");

    atmos->rainf = calloc(NR + 1, sizeof(*(atmos->rainf)));
    check_alloc_status(atmos->prec, "Memory allocation error.");

    atmos->Qair = calloc(NR + 1, sizeof(*(atmos->Qair)));
    check_alloc_status(atmos->Qair, "Memory allocation error.");

    atmos->rel_humid = calloc(NR + 1, sizeof(*(atmos->rel_humid)));
    check_alloc_status(atmos->rel_humid, "Memory allocation error.");

    atmos->wind = calloc(NR + 1, sizeof(*(atmos->wind)));
    check_alloc_status(atmos->wind, "Memory allocation error.");

    if (options.CARBON) {
        atmos->Catm = calloc(NR + 1, sizeof(*(atmos->Catm)));
        check_alloc_status(atmos->Catm, "Memory allocation error.");

        atmos->fdir = calloc(NR + 1, sizeof(*(atmos->fdir)));
        check_alloc_status(atmos->fdir, "Memory allocation error.");

        atmos->par = calloc(NR + 1, sizeof(*(atmos->par)));
        check_alloc_status(atmos->par, "Memory allocation error.");
    }
    if (options.ROUT) {
        atmos->channel_in = calloc(NR + 1, sizeof(*(atmos->channel_in)));
        check_alloc_status(atmos->channel_in, "Memory allocation error.");
    }
}

/******************************************************************************
 * @brief    Free memory for the atmos data structure.
 *****************************************************************************/
void
free_atmos(atmos_data_struct *atmos)

{
    extern option_struct options;

    if (atmos == NULL) {
        return;
    }

    free(atmos->air_temp);
    free(atmos->longwave);
    free(atmos->prec);
    free(atmos->pressure);
    free(atmos->shortwave);
    free(atmos->snowf);
    free(atmos->rainf);
    free(atmos->Qair);
    free(atmos->rel_humid);
    free(atmos->wind);
    if (options.CARBON) {
        free(atmos->Catm);
        free(atmos->fdir);
        free(atmos->par);
    }
    if (options.ROUT) {
        free(atmos->channel_in);
    }
}

