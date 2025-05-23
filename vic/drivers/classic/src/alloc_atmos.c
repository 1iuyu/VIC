/******************************************************************************
 * @section DESCRIPTION
 *
 * Allocate and free memory for the atmos data struct
 *****************************************************************************/

#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    Allocate memory for the atmos data structure.
 *****************************************************************************/
void
alloc_atmos(int                 nrecs,
            force_data_struct **force)
{
    extern option_struct options;
    int                  i;

    *force = calloc(nrecs, sizeof(force_data_struct));
    check_alloc_status(*force, "Memory allocation error.");

    for (i = 0; i < nrecs; i++) {
        (*force)[i].air_temp = calloc(NR + 1, sizeof(*(*force)[i].air_temp));
        check_alloc_status((*force)[i].air_temp, "Memory allocation error.");
        (*force)[i].density = calloc(NR + 1, sizeof(*(*force)[i].density));
        check_alloc_status((*force)[i].density, "Memory allocation error.");
        (*force)[i].longwave = calloc(NR + 1, sizeof(*(*force)[i].longwave));
        check_alloc_status((*force)[i].longwave, "Memory allocation error.");
        (*force)[i].prec = calloc(NR + 1, sizeof(*(*force)[i].prec));
        check_alloc_status((*force)[i].prec, "Memory allocation error.");
        (*force)[i].rainf = calloc(NR + 1, sizeof(*(*force)[i].rainf));
        check_alloc_status((*force)[i].rainf, "Memory allocation error.");
        (*force)[i].snowf = calloc(NR + 1, sizeof(*(*force)[i].snowf));
        check_alloc_status((*force)[i].snowf, "Memory allocation error.");
        (*force)[i].pressure = calloc(NR + 1, sizeof(*(*force)[i].pressure));
        check_alloc_status((*force)[i].pressure, "Memory allocation error.");
        (*force)[i].shortwave = calloc(NR + 1, sizeof(*(*force)[i].shortwave));
        check_alloc_status((*force)[i].shortwave, "Memory allocation error.");
        (*force)[i].snowflag = calloc(NR + 1, sizeof(*(*force)[i].snowflag));
        check_alloc_status((*force)[i].snowflag, "Memory allocation error.");
        (*force)[i].vp = calloc(NR + 1, sizeof(*(*force)[i].vp));
        check_alloc_status((*force)[i].vp, "Memory allocation error.");
        (*force)[i].vpd = calloc(NR + 1, sizeof(*(*force)[i].vpd));
        check_alloc_status((*force)[i].vpd, "Memory allocation error.");
        (*force)[i].wind = calloc(NR + 1, sizeof(*(*force)[i].wind));
        check_alloc_status((*force)[i].wind, "Memory allocation error.");
        if (options.LAKES) {
            (*force)[i].channel_in =
                calloc(NR + 1, sizeof(*(*force)[i].channel_in));
            check_alloc_status((*force)[i].channel_in,
                               "Memory allocation error.");
        }
        if (options.CARBON) {
            (*force)[i].Catm = calloc(NR + 1, sizeof(*(*force)[i].Catm));
            check_alloc_status((*force)[i].Catm, "Memory allocation error.");
            (*force)[i].coszen = calloc(NR + 1, sizeof(*(*force)[i].coszen));
            check_alloc_status((*force)[i].coszen, "Memory allocation error.");
            (*force)[i].fdir = calloc(NR + 1, sizeof(*(*force)[i].fdir));
            check_alloc_status((*force)[i].fdir, "Memory allocation error.");
            (*force)[i].par = calloc(NR + 1, sizeof(*(*force)[i].par));
            check_alloc_status((*force)[i].par, "Memory allocation error.");
        }
    }
}

/******************************************************************************
 * @brief    Free memory for the atmos data structure.
 *****************************************************************************/
void
free_atmos(int                 nrecs,
           force_data_struct **force)
{
    extern option_struct options;
    int                  i;

    if (*force == NULL) {
        return;
    }

    for (i = 0; i < nrecs; i++) {
        free((*force)[i].air_temp);
        free((*force)[i].density);
        free((*force)[i].longwave);
        free((*force)[i].prec);
        free((*force)[i].rainf);
        free((*force)[i].snowf);
        free((*force)[i].pressure);
        free((*force)[i].shortwave);
        free((*force)[i].snowflag);
        free((*force)[i].vp);
        free((*force)[i].vpd);
        free((*force)[i].wind);
        if (options.LAKES) {
            free((*force)[i].channel_in);
        }
        if (options.CARBON) {
            free((*force)[i].Catm);
            free((*force)[i].coszen);
            free((*force)[i].fdir);
            free((*force)[i].par);
        }
    }

    free(*force);
}
