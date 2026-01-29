/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine reads in vegetation parameters for the current grid cell.
 *
 * It also relates each vegetation class in the cell to the appropriate
 * parameters in the vegetation library.
 *****************************************************************************/

#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    Read vegetation parameters.
 *****************************************************************************/
veg_con_struct *
read_vegparam(FILE  *vegparam,
              int    gridcel,
              size_t Nveg_type)
{
    void ttrim(char *string);
    extern veg_lib_struct   *veg_lib;
    extern option_struct     options;

    veg_con_struct          *temp;
    size_t                   j;
    int                      vegetat_type_num;
    int                      vegcel, i, k, skip, veg_class;
    int                      MaxVeg;
    int                      Nfields, NfieldsMax;
    double                   Cv_sum;
    char                     str[MAX_VEGPARAM_LINE_LENGTH];
    char                     line[MAXSTRING];
    char                     tmpline[MAXSTRING];
    const char               delimiters[] = " \t";
    char                    *token;
    char                    *vegarr[MAX_VEGPARAM_LINE_LENGTH];
    size_t                   length;
    size_t                   cidx;
    double                   tmp;

    skip = 1;
    if (options.VEGPARAM_LAI) {
        skip++;
    }
    if (options.VEGPARAM_FCAN) {
        skip++;
    }
    if (options.VEGPARAM_ALB) {
        skip++;
    }

    while ((fscanf(vegparam, "%d %d", &vegcel,
                   &vegetat_type_num) == 2) && vegcel != gridcel) {
        if (vegetat_type_num < 0) {
            log_err("number of vegetation tiles (%i) given for cell %i "
                    "is < 0.", vegetat_type_num, vegcel);
        }
        for (i = 0; i <= vegetat_type_num * skip; i++) {
            if (fgets(str, MAX_VEGPARAM_LINE_LENGTH, vegparam) == NULL) {
                log_err("unexpected EOF for cell %i while reading root zones "
                        "and LAI", vegcel);
            }
        }
    }
    fgets(str, MAX_VEGPARAM_LINE_LENGTH, vegparam); // read newline at end of veg class line to advance to next line
    if (vegcel != gridcel) {
        log_err("Grid cell %d not found", gridcel);
    }

    // Make sure to allocate extra memory for bare soil tile
    // and optionally an above-treeline veg tile
    MaxVeg = vegetat_type_num + 1;

    /** Allocate memory for vegetation grid cell parameters **/
    temp = calloc(MaxVeg, sizeof(*temp));
    Cv_sum = 0.0;

    for (i = 0; i < vegetat_type_num; i++) {

        temp[i].vegetat_type_num = vegetat_type_num;

        /* Upper boundaries of canopy layers, expressed in terms of fraction of total LAI  */
        if (options.CARBON) {
            temp[i].CanopLayerBnd = calloc(options.Ncanopy,
                                           sizeof(*(temp[i].CanopLayerBnd)));
            for (cidx = 0; cidx < options.Ncanopy; cidx++) {
                /* apportion LAI equally among layers */
                temp[i].CanopLayerBnd[cidx] =
                    (double) ((cidx + 1)) / (double) (options.Ncanopy);
            }
        }

        // Read the root zones line
        if (fgets(line, MAXSTRING, vegparam) == NULL) {
            log_err("unexpected EOF for cell %i while reading "
                    "vegetat_type_num %d", vegcel, vegetat_type_num);
        }
        strcpy(tmpline, line);
        ttrim(tmpline);
        token = strtok(tmpline, delimiters); /*  token => veg_class, move 'line' pointer to next field */
        Nfields = 0;
        vegarr[Nfields] =
            calloc(MAX_VEGPARAM_LINE_LENGTH, sizeof(*(vegarr[Nfields])));
        strcpy(vegarr[Nfields], token);
        Nfields++;

        token = strtok(NULL, delimiters);
        while (token != NULL && (length = strlen(token)) == 0) {
            token = strtok(NULL, delimiters);
        }
        while (token != NULL) {
            vegarr[Nfields] = calloc(MAX_VEGPARAM_LINE_LENGTH,
                                     sizeof(*(vegarr[Nfields])));
            strcpy(vegarr[Nfields], token);
            Nfields++;
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
        }

        NfieldsMax = 6; /* Number of expected fields this line. (veg_class, Cv, a, b, d, bandindex) */
        if (options.BLOWING) {
            NfieldsMax += 3;
        }
        if (Nfields != NfieldsMax) {
            log_err("Cell %d - expecting %d fields but found %d in veg line %s",
                    gridcel, NfieldsMax, Nfields, line);
        }

        temp[i].veg_class = atoi(vegarr[0]);
        temp[i].Cv = atof(vegarr[1]);
        temp[i].a = atof(vegarr[2]);
        temp[i].b = atof(vegarr[3]);
        temp[i].d = atof(vegarr[4]);
        temp[i].BandIndex = atof(vegarr[5]);

        veg_class = MISSING;
        for (j = 0; j < Nveg_type; j++) {
            if (temp[i].veg_class == veg_lib[j].veg_class) {
                veg_class = j;
            }
        }
        if (veg_class == MISSING) {
            log_err("The vegetation class id %i in vegetation tile %i from "
                    "cell %i is not defined in the vegetation library file.",
                    temp[i].veg_class, i, gridcel);
        }
        else {
            temp[i].veg_class = veg_class;
        }

        if (veg_class == options.GLACIER_ID) {
            temp[i].IS_GLAC = true;
        }
        else {
            temp[i].IS_GLAC = false;

        }
        
        Cv_sum += temp[i].Cv;

        for (k = 0; k < Nfields; k++) {
            free(vegarr[k]);
        }

        for (j = 0; j < MONTHS_PER_YEAR; j++) {
            temp[i].albedo[j] = veg_lib[temp[i].veg_class].albedo[j];
            temp[i].displacement[j] =
                veg_lib[temp[i].veg_class].displacement[j];
            temp[i].fcanopy[j] = veg_lib[temp[i].veg_class].fcanopy[j];
            temp[i].LAI[j] = veg_lib[temp[i].veg_class].LAI[j];
            temp[i].SAI[j] = veg_lib[temp[i].veg_class].SAI[j];
            temp[i].roughness[j] = veg_lib[temp[i].veg_class].roughness[j];
        }

        if (options.VEGPARAM_LAI) {
            // Read the LAI line
            if (fgets(line, MAXSTRING, vegparam) == NULL) {
                log_err("Unexpected EOF for cell %i while reading LAI for "
                        "vegetat_type_num %d", vegcel, vegetat_type_num);
            }
            Nfields = 0;
            vegarr[Nfields] = calloc(MAX_VEGPARAM_LINE_LENGTH,
                                     sizeof(*(vegarr[Nfields])));
            strcpy(tmpline, line);
            ttrim(tmpline);
            token = strtok(tmpline, delimiters);
            strcpy(vegarr[Nfields], token);
            Nfields++;

            while ((token = strtok(NULL, delimiters)) != NULL) {
                vegarr[Nfields] = calloc(MAX_VEGPARAM_LINE_LENGTH,
                                         sizeof(*(vegarr[Nfields])));
                strcpy(vegarr[Nfields], token);
                Nfields++;
            }
            NfieldsMax = MONTHS_PER_YEAR; /* For LAI */
            if (Nfields != NfieldsMax) {
                log_err("cell %d - expecting %d LAI values but found "
                        "%d in line %s", gridcel, NfieldsMax, Nfields, line);
            }

            if (options.LAI_SRC == FROM_VEGPARAM) {
                for (j = 0; j < MONTHS_PER_YEAR; j++) {
                    tmp = atof(vegarr[j]);
                    if (tmp != NODATA_VH) {
                        temp[i].LAI[j] = tmp;
                    }
                }
            }
            for (k = 0; k < Nfields; k++) {
                free(vegarr[k]);
            }
        }

        if (options.VEGPARAM_FCAN) {
            // Read the fcanopy line
            if (fgets(line, MAXSTRING, vegparam) == NULL) {
                log_err("unexpected EOF for cell %i while reading fcanopy "
                        "for vegetat_type_num %d", vegcel, vegetat_type_num);
            }
            Nfields = 0;
            vegarr[Nfields] = calloc(MAX_VEGPARAM_LINE_LENGTH,
                                     sizeof(*(vegarr[Nfields])));
            strcpy(tmpline, line);
            ttrim(tmpline);
            token = strtok(tmpline, delimiters);
            strcpy(vegarr[Nfields], token);
            Nfields++;

            while ((token = strtok(NULL, delimiters)) != NULL) {
                vegarr[Nfields] = calloc(MAX_VEGPARAM_LINE_LENGTH,
                                         sizeof(*(vegarr[Nfields])));
                strcpy(vegarr[Nfields], token);
                Nfields++;
            }
            NfieldsMax = MONTHS_PER_YEAR; /* For fcanopy */
            if (Nfields != NfieldsMax) {
                log_err("cell %d - expecting %d fcanopy values but found %d "
                        "in line %s", gridcel, NfieldsMax, Nfields, line);
            }

            if (options.FCAN_SRC == FROM_VEGPARAM) {
                for (j = 0; j < MONTHS_PER_YEAR; j++) {
                    tmp = atof(vegarr[j]);
                    if (tmp != NODATA_VH) {
                        temp[i].fcanopy[j] = tmp;
                    }
                }
            }
            for (k = 0; k < Nfields; k++) {
                free(vegarr[k]);
            }
        }

        if (options.VEGPARAM_ALB) {
            // Read the albedo line
            if (fgets(line, MAXSTRING, vegparam) == NULL) {
                log_err("unexpected EOF for cell %i while reading albedo for "
                        "vegetat_type_num %d", vegcel, vegetat_type_num);
            }
            Nfields = 0;
            vegarr[Nfields] = calloc(MAX_VEGPARAM_LINE_LENGTH,
                                     sizeof(*(vegarr[Nfields])));
            strcpy(tmpline, line);
            ttrim(tmpline);
            token = strtok(tmpline, delimiters);
            strcpy(vegarr[Nfields], token);
            Nfields++;

            while ((token = strtok(NULL, delimiters)) != NULL) {
                vegarr[Nfields] = calloc(MAX_VEGPARAM_LINE_LENGTH,
                                         sizeof(*(vegarr[Nfields])));
                strcpy(vegarr[Nfields], token);
                Nfields++;
            }
            NfieldsMax = MONTHS_PER_YEAR; /* For albedo */
            if (Nfields != NfieldsMax) {
                log_err("cell %d - expecting %d albedo values but found %d in "
                        "line %s", gridcel, NfieldsMax, Nfields, line);
            }

            if (options.ALB_SRC == FROM_VEGPARAM) {
                for (j = 0; j < MONTHS_PER_YEAR; j++) {
                    tmp = atof(vegarr[j]);
                    if (tmp != NODATA_VH) {
                        temp[i].albedo[j] = tmp;
                    }
                }
            }
            for (k = 0; k < Nfields; k++) {
                free(vegarr[k]);
            }
        }
    }

    // Determine if we have bare soil
    if (Cv_sum > 1.0) {
        log_warn("Cv_sum exceeds 1.0 (%f) at grid cell %d, fractions being "
                 "adjusted to equal 1", Cv_sum, gridcel);
        for (j = 0; j < (size_t)vegetat_type_num; j++) {
            temp[j].Cv = temp[j].Cv / Cv_sum;
        }
        Cv_sum = 1.;
    }
    else if (Cv_sum > 0.99 && Cv_sum < 1.0) {
        log_warn("Cv > 0.99 and Cv < 1.0 at grid cell %d, model "
                 "assuming that bare soil is not to be run - fractions being "
                 "adjusted to equal 1",
                 gridcel);
        for (j = 0; j < (size_t)vegetat_type_num; j++) {
            temp[j].Cv = temp[j].Cv / Cv_sum;
        }
        Cv_sum = 1.;
    }

    // Default bare soil tile - not specified in vegparam file
    i = vegetat_type_num;
    temp[i].veg_class = Nveg_type; // Create a veg_class ID for bare soil, which is not mentioned in the veg library
    temp[i].Cv = 1.0 - Cv_sum;
    if (temp[i].Cv < 0) {
        temp[i].Cv = 0;
    }
    for (j = 0; j < MONTHS_PER_YEAR; j++) {
        temp[i].albedo[j] = veg_lib[temp[i].veg_class].albedo[j];
        temp[i].displacement[j] =
            veg_lib[temp[i].veg_class].displacement[j];
        temp[i].fcanopy[j] = veg_lib[temp[i].veg_class].fcanopy[j];
        temp[i].LAI[j] = veg_lib[temp[i].veg_class].LAI[j];
        temp[i].roughness[j] = veg_lib[temp[i].veg_class].roughness[j];
    }

    return temp;
}

/* trim trailing newlines */

#define END '\0'
#define NEW '\n'

/******************************************************************************
 * @brief    trim trailing newlines
 *****************************************************************************/
void
ttrim(char *c)
{
    while ((*c++ != END)) {
        ;
    }
    for (--c; *--c == NEW; *c = END) {
        ;
    }
}
