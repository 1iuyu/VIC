/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine reads in a library of vegetation parameters for all vegetation
 * classes used in the model.  The veg class number is used to reference the
 * information in this library.
 *****************************************************************************/

#include "vic_driver_image.h"

/******************************************************************************
 * @brief    Read in a library of vegetation parameters for all vegetation
 *           classes used in the model.
 *****************************************************************************/
veg_lib_struct *
read_veglib(FILE   *veglib,
            size_t *Ntype)
{
    extern option_struct       options;

    veg_lib_struct            *temp;
    size_t                     i, j;
    size_t                     Nveg_type;
    char                       str[MAXSTRING];
    char                       tmpstr[MAXSTRING];
    double                     tmp_double;

    rewind(veglib);
    fgets(str, MAXSTRING, veglib);
    Nveg_type = 0;
    while (!feof(veglib)) {
        if (str[0] <= 57 && str[0] >= 48) {
            Nveg_type++;
        }
        fgets(str, MAXSTRING, veglib);
    }
    rewind(veglib);

    // +1 for bare soil
    temp = calloc(Nveg_type + 1, sizeof(*temp));
    options.NVEGTYPES = Nveg_type + 1;

    fscanf(veglib, "%s", str);
    i = 0;
    while (!feof(veglib)) {
        if (str[0] <= 57 && str[0] >= 48) {
            temp[i].NVegLibTypes = Nveg_type;
            temp[i].veg_class = atoi(str);
            fscanf(veglib, "%lf", &temp[i].Canopy_Upper);
            fscanf(veglib, "%lf", &temp[i].Canopy_Lower);
            fscanf(veglib, "%lf", &temp[i].Canopy_Radius);
            fscanf(veglib, "%lf", &temp[i].COI);
            fscanf(veglib, "%lf", &temp[i].c_biomass);
            fscanf(veglib, "%lf", &temp[i].d_leaf);
            fscanf(veglib, "%lf", &temp[i].root_a);      /**< Empirical parameter a in eqa(2) */
            fscanf(veglib, "%lf", &temp[i].root_b);      /**< Empirical parameter b in eqa(2) */
            fscanf(veglib, "%lf", &temp[i].root_d);      /**< Maximum root depth (m) */
            fscanf(veglib, "%lf", &temp[i].liq_bioms);
            fscanf(veglib, "%lf", &temp[i].slatop);
            fscanf(veglib, "%lf", &temp[i].stem_num);
            fscanf(veglib, "%lf", &temp[i].Z0sub_LAImax);
            fscanf(veglib, "%lf", &temp[i].Z0sub_Cs);
            fscanf(veglib, "%lf", &temp[i].Z0sub_Cr);
            fscanf(veglib, "%lf", &temp[i].Z0sub_c);
            fscanf(veglib, "%lf", &temp[i].Z0sub_cw);
            fscanf(veglib, "%lf", &temp[i].smpsc);
            fscanf(veglib, "%lf", &temp[i].smpso);

            for (j = 0; j < MONTHS_PER_YEAR; j++) {
                fscanf(veglib, "%lf", &temp[i].LAI[j]);
                if (temp[i].LAI[j] < 0) {
                    log_err("leaf area index must be > 0 " "(%f)",
                             temp[i].LAI[j]);
                }
            }
            for (j = 0; j < MONTHS_PER_YEAR; j++) {
                fscanf(veglib, "%lf", &temp[i].SAI[j]);
                if (temp[i].SAI[j] < 0) {
                    log_err("stem area index must be > 0 " "(%f)",
                             temp[i].SAI[j]);
                }
            }
            /* Default values of fcanopy */
            for (j = 0; j < MONTHS_PER_YEAR; j++) {
                temp[i].fcanopy[j] = 1.00;
            }
            if (options.VEGLIB_FCAN) {
                for (j = 0; j < MONTHS_PER_YEAR; j++) {
                    fscanf(veglib, "%lf", &tmp_double);
                    if (options.FCAN_SRC != FROM_DEFAULT) {
                        temp[i].fcanopy[j] = tmp_double;
                        if (temp[i].fcanopy[j] < 0 ||
                            temp[i].fcanopy[j] > 1) {
                            log_err(
                                "Veg cover fraction must be between 0 and 1 " "(%f)",
                                temp[i].fcanopy[j]);
                        }
                    }
                }
            }
            for (j = 0; j < MAX_SWBANDS; j++) {
                fscanf(veglib, "%lf", &temp[i].reflleaf[j]);
            }
            for (j = 0; j < MAX_SWBANDS; j++) {
                fscanf(veglib, "%lf", &temp[i].reflstem[j]);
            }
            for (j = 0; j < MAX_SWBANDS; j++) {
                fscanf(veglib, "%lf", &temp[i].transleaf[j]);
            }
            for (j = 0; j < MAX_SWBANDS; j++) {
                fscanf(veglib, "%lf", &temp[i].transstem[j]);
            }
            fscanf(veglib, "%lf", &temp[i].trunk_dia);
            /* Carbon-cycling parameters */
            if (options.VEGLIB_PHOTO) {
                fscanf(veglib, "%s", tmpstr); /* photosynthetic pathway */
                if (!strcmp(tmpstr, "0") || !strcmp(tmpstr, "C3")) {
                    temp[i].Ctype = PHOTO_C3;
                }
                else if (!strcmp(tmpstr, "1") || !strcmp(tmpstr, "C4")) {
                    temp[i].Ctype = PHOTO_C4;
                }
                if (!strcmp(tmpstr, "C3") || !strcmp(tmpstr, "C4")) {
                    log_warn("Use of strings \"C3\" and \"C4\" as values of "
                             "Ctype is deprecated.  Please replace these with "
                             "\"0\" and \"1\", respectively");
                }
                fscanf(veglib, "%lf", &temp[i].m_bb);
                fscanf(veglib, "%lf", &temp[i].matric50);
                fscanf(veglib, "%lf", &temp[i].kseg_max);
                fscanf(veglib, "%lf", &temp[i].kcano_max);
                fscanf(veglib, "%lf", &temp[i].kroot_max);
                fscanf(veglib, "%lf", &temp[i].theta_cj);
                fscanf(veglib, "%lf", &temp[i].leaf_CN);
                fscanf(veglib, "%lf", &temp[i].SLA_top);
                fscanf(veglib, "%lf", &temp[i].fN_rub);
                fscanf(veglib, "%lf", &temp[i].medlynslope);
                fscanf(veglib, "%lf", &temp[i].medlynint);
            }

            fgets(str, MAXSTRING, veglib); /* skip over end of line comments */
            i++;
        }
        else {
            fgets(str, MAXSTRING, veglib);
        }
        fscanf(veglib, "%s", str);
    }
    if (i != Nveg_type) {
        log_err("Problem reading vegetation library file - make sure "
                "the file has the right number of columns.");
    }
    *Ntype = Nveg_type;

    // Assign properties of bare soil to default bare soil tile
    temp[i].NVegLibTypes = Nveg_type;
    temp[i].veg_class = Nveg_type + 1;
    for (j = 0; j < MONTHS_PER_YEAR; j++) {
        temp[i].LAI[j] = 0.0;
        temp[i].SAI[j] = 0.0;
        temp[i].fcanopy[j] = MIN_FCANOPY;
        // These will be assigned in read_soilparam.c

    }
    temp[i].trunk_dia = 0.0;
    if (options.VEGLIB_PHOTO) {
        temp[i].Ctype = PHOTO_C3;
    }

    return temp;
}

/******************************************************************************
 * @brief    This routine frees the veglib structure.
 *****************************************************************************/
void
free_veglib(veg_lib_struct **veg_lib)
{
    free((char*)(*veg_lib));
}
