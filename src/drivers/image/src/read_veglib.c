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
            fscanf(veglib, "%lf", &temp[i].trunk_dia);

            fscanf(veglib, "%s", tmpstr); /* landunit type */
            if (!strcmp(tmpstr, "0") || !strcmp(tmpstr, "SOIL")) {
                temp[i].Landtype = LAND_SOIL;
            }
            else if (!strcmp(tmpstr, "1") || !strcmp(tmpstr, "GLAC")) {
                temp[i].Landtype = LAND_GLAC;
            }
            else if (!strcmp(tmpstr, "2") || !strcmp(tmpstr, "WET")) {
                temp[i].Landtype = LAND_WET;
            }
            else if (!strcmp(tmpstr, "3") || !strcmp(tmpstr, "URBAN")) {
                temp[i].Landtype = LAND_URBAN;
            }

            /* 兼容性警告：字符串形式已弃用，建议使用数字 */
            if (!strcmp(tmpstr, "VEG") || !strcmp(tmpstr, "GLAC") ||
                !strcmp(tmpstr, "WET") || !strcmp(tmpstr, "URBAN")) {
                log_warn("Use of strings (e.g., \"VEG\", \"GLAC\") as values of "
                        "landtype is deprecated.  Please replace these with "
                        "numeric codes: 0 (SOIL), 1 (GLAC), 2 (WET), or 3 (URBAN)");
            }

            /* Carbon-cycling parameters */
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
            fscanf(veglib, "%lf", &temp[i].froot_leaf);
            fscanf(veglib, "%lf", &temp[i].theta_cj);
            fscanf(veglib, "%lf", &temp[i].kcano_max);
            fscanf(veglib, "%lf", &temp[i].kroot_max);
            fscanf(veglib, "%lf", &temp[i].matric50);
            fscanf(veglib, "%lf", &temp[i].leaf_CN);
            fscanf(veglib, "%lf", &temp[i].SLA_top);
            fscanf(veglib, "%lf", &temp[i].fN_rub);
            fscanf(veglib, "%lf", &temp[i].medlynslope);
            fscanf(veglib, "%lf", &temp[i].medlynint);

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
    }
    temp[i].Canopy_Upper = 0.0;
    temp[i].Canopy_Lower = 0.0;
    temp[i].Canopy_Radius = 0.0;
    temp[i].COI = 0.0;
    temp[i].c_biomass = 0.0;
    temp[i].d_leaf = 0.0;
    temp[i].root_a = 0.0;
    temp[i].root_b = 0.0;
    temp[i].root_d = 0.0;
    temp[i].liq_bioms = 0.0;
    temp[i].slatop = 0.0;
    temp[i].stem_num = 0.0;
    temp[i].Z0sub_c = 0.0;
    temp[i].Z0sub_Cr = 0.0;
    temp[i].Z0sub_Cs = 0.0;
    temp[i].Z0sub_cw = 0.0;
    temp[i].Z0sub_LAImax = 0.0;
    temp[i].smpsc = 0.0;
    temp[i].smpso = 0.0;
    temp[i].trunk_dia = 0.0;
    temp[i].Landtype = LAND_SOIL;
    temp[i].Ctype = PHOTO_C3;
    temp[i].froot_leaf = 0.0;
    temp[i].theta_cj = 0.0;
    temp[i].kcano_max = 0.0;
    temp[i].kroot_max = 0.0;
    temp[i].matric50 = 0.0;
    temp[i].leaf_CN = 0.0;
    temp[i].SLA_top = 0.0;
    temp[i].fN_rub = 0.0;
    temp[i].medlynint = 0.0;
    temp[i].medlynslope = 0.0;

    return (temp);
}

/******************************************************************************
 * @brief    This routine frees the veglib structure.
 *****************************************************************************/
void
free_veglib(veg_lib_struct **veg_lib)
{
    free(*veg_lib);
}
