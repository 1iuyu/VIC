/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine reads soil parameters for each grid cell.
 *****************************************************************************/

#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    This routine reads soil parameters for each grid cell.
 *****************************************************************************/
void
read_soilparam(FILE            *soilparam,
               soil_con_struct *temp,
               bool            *RUN_MODEL,
               bool            *MODEL_DONE)
{
    void ttrim(char *string);
    extern option_struct     options;
    extern veg_lib_struct   *veg_lib;
    extern parameters_struct param;

    char                     line[MAXSTRING];
    char                     tmpline[MAXSTRING];
    const char               delimiters[] = " \t";
    char                    *token;
    size_t                   layer;
    int                      tempint, j;
    double                   off_gmt;
    size_t                   length;
    int                      Nbands, band;
    int                      flag;
    size_t                   k;

    /** Read plain ASCII soil parameter file **/
    if ((fscanf(soilparam, "%d", &flag)) != EOF) {
        if (flag) {
            *RUN_MODEL = true;
        }
        else {
            *RUN_MODEL = false;
        }

        if (fgets(line, MAXSTRING, soilparam) == NULL) {
            log_err("Unexpected EOF while reading soil file");
        }
    }
    else {
        *MODEL_DONE = true;
        *RUN_MODEL = false;
    }

    if (!(*MODEL_DONE) && (*RUN_MODEL)) {
        strcpy(tmpline, line);
        ttrim(tmpline);
        if ((token = strtok(tmpline, delimiters)) == NULL) {
            log_err("Can't find values for CELL NUMBER in soil file");
        }
        sscanf(token, "%d", &(temp->gridcel));
        token = strtok(NULL, delimiters);
        while (token != NULL && (length = strlen(token)) == 0) {
            token = strtok(NULL, delimiters);
        }
        if (token == NULL) {
            log_err("Can't find values for CELL LATITUDE in soil file");
        }
        sscanf(token, "%lf", &(temp->lat));
        token = strtok(NULL, delimiters);
        while (token != NULL && (length = strlen(token)) == 0) {
            token = strtok(NULL, delimiters);
        }
        if (token == NULL) {
            log_err("Can't find values for CELL LONGITUDE in soil file");
        }
        sscanf(token, "%lf", &(temp->lng));

        /* read infiltration parameter */
        token = strtok(NULL, delimiters);
        while (token != NULL && (length = strlen(token)) == 0) {
            token = strtok(NULL, delimiters);
        }
        if (token == NULL) {
            log_err("Can't find values for INFILTRATION in soil file");
        }
        sscanf(token, "%lf", &(temp->b_infilt));
        if (temp->b_infilt <= 0) {
            log_err("b_infilt (%f) in soil file is <= 0; b_infilt must "
                    "be positive", temp->b_infilt);
        }

        /* read heterogeniety parameter for infiltration */
        token = strtok(NULL, delimiters);
        while (token != NULL && (length = strlen(token)) == 0) {
            token = strtok(NULL, delimiters);
        }
        if (token == NULL) {
            log_err("Can't find values for dynamic VIC heterogeniety parameter "
                    "in soil file");
        }
        sscanf(token, "%lf", &(temp->b_dynamic));

        /* read mean capilary drive (m) for dynamic VIC runoff */
        token = strtok(NULL, delimiters);
        while (token != NULL && (length = strlen(token)) == 0) {
            token = strtok(NULL, delimiters);
        }
        if (token == NULL) {
            log_err("Can't find values for MEAN CAPILARY DRIVE in "
                    "soil file");
        }
        sscanf(token, "%lf", &(temp->capil_drive));

        /* read bexp for each layer */
        for (layer = 0; layer < options.Nlayer; layer++) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for B_EXP for layer %zu in "
                        "soil file", layer);
            }
            sscanf(token, "%lf", &(temp->expt_node)[layer]);
        }

        /* read expt for each layer */
        for (layer = 0; layer < options.Nlayer; layer++) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for EXPT for layer %zu in "
                        "soil file", layer);
            }
            sscanf(token, "%lf", &(temp->expt_node)[layer]);
        }

        /* read layer saturated hydraulic conductivity */
        for (layer = 0; layer < options.Nlayer; layer++) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for SATURATED HYDRAULIC "
                        "CONDUCTIVITY for layer %zu in soil file", layer);
            }
            sscanf(token, "%lf", &(temp->Ksat_node)[layer]);
        }

        /* read layer saturated hydraulic diffusivity */
        for (layer = 0; layer < options.Nlayer; layer++) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for Dsat for layer %zu in "
                        "soil file", layer);
            }
            sscanf(token, "%lf", &(temp->Dsat)[layer]);
        }

        /* read cell mean elevation */
        token = strtok(NULL, delimiters);
        while (token != NULL && (length = strlen(token)) == 0) {
            token = strtok(NULL, delimiters);
        }
        if (token == NULL) {
            log_err("Can't find values for CELL MEAN ELEVATION in soil "
                    "file");
        }
        sscanf(token, "%lf", &(temp->elevation));

        /* soil layer thicknesses */
        for (layer = 0; layer < options.Nlayer; layer++) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for LAYER THICKNESS for "
                        "layer %zu in soil file", layer);
            }
            sscanf(token, "%lf", &(temp->depth)[layer]);
        }
        /* round soil layer thicknesses to nearest mm */
        for (layer = 0; layer < options.Nlayer; layer++) {
            temp->depth[layer] =
                round(temp->depth[layer] * MM_PER_M) / MM_PER_M;
        }

        /* read soil damping depth */
        token = strtok(NULL, delimiters);
        while (token != NULL && (length = strlen(token)) == 0) {
            token = strtok(NULL, delimiters);
        }
        if (token == NULL) {
            log_err("Can't find values for SOIL DAMPING DEPTH in soil "
                    "file");
        }
        sscanf(token, "%lf", &(temp->dp));

        /* read layer bulk density */
        for (layer = 0; layer < options.Nlayer; layer++) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for mineral BULK DENSITY "
                        "for layer %zu in soil file", layer);
            }
            sscanf(token, "%lf", &(temp->bulk_dens_min)[layer]);
            if (temp->bulk_dens_min[layer] <= 0) {
                log_err("layer %zu mineral bulk density (%f) must "
                        "be > 0", layer, temp->bulk_dens_min[layer]);
            }
        }

        /* read layer soil density */
        for (layer = 0; layer < options.Nlayer; layer++) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for mineral SOIL DENSITY "
                        "for layer %zu in soil file", layer);
            }
            sscanf(token, "%lf", &(temp->soil_dens_min)[layer]);
            if (temp->soil_dens_min[layer] <= 0) {
                log_err("layer %zu mineral soil density (%f) must "
                        "be > 0", layer, temp->soil_dens_min[layer]);
            }
            if (temp->bulk_dens_min[layer] >= temp->soil_dens_min[layer]) {
                log_err("layer %zu mineral bulk density (%f) must "
                        "be less than mineral soil density (%f)", layer,
                        temp->bulk_dens_min[layer],
                        temp->soil_dens_min[layer]);
            }
        }

        if (options.ORGANIC_FRACT) {
            /* read layer organic content */
            for (layer = 0; layer < options.Nlayer; layer++) {
                token = strtok(NULL, delimiters);
                while (token != NULL && (length = strlen(token)) == 0) {
                    token = strtok(NULL, delimiters);
                }
                if (token == NULL) {
                    log_err("Can't find values for ORGANIC CONTENT for "
                            "layer %zu in soil file", layer);
                }
                sscanf(token, "%lf", &(temp->organic)[layer]);
                if (temp->organic[layer] > 1. || temp->organic[layer] < 0) {
                    log_err("Need valid volumetric organic soil "
                            "fraction when options.ORGANIC_FRACT is set "
                            "to true.  %f is not acceptable.",
                            temp->organic[layer]);
                }
            }

            /* read layer bulk density */
            for (layer = 0; layer < options.Nlayer; layer++) {
                token = strtok(NULL, delimiters);
                while (token != NULL && (length = strlen(token)) == 0) {
                    token = strtok(NULL, delimiters);
                }
                if (token == NULL) {
                    log_err("Can't find values for organic BULK "
                            "DENSITY for layer %zu in soil file", layer);
                }
                sscanf(token, "%lf", &(temp->bulk_dens_org)[layer]);
                if (temp->bulk_dens_org[layer] <= 0 && temp->organic[layer] >
                    0) {
                    log_warn("layer %zu organic bulk density (%f) must "
                             "be > 0; setting to mineral bulk density "
                             "(%f)", layer, temp->bulk_dens_org[layer],
                             temp->bulk_dens_min[layer]);
                    temp->bulk_dens_org[layer] = temp->bulk_dens_min[layer];
                }
            }

            /* read layer soil density */
            for (layer = 0; layer < options.Nlayer; layer++) {
                token = strtok(NULL, delimiters);
                while (token != NULL && (length = strlen(token)) == 0) {
                    token = strtok(NULL, delimiters);
                }
                if (token == NULL) {
                    log_err("Can't find values for organic SOIL DENSITY for "
                            "layer %zu in soil file", layer);
                }
                sscanf(token, "%lf", &(temp->soil_dens_org)[layer]);
                if (temp->soil_dens_org[layer] <= 0 && temp->organic[layer] >
                    0) {
                    log_warn("layer %zu organic soil density (%f) must be "
                             "> 0; setting to mineral soil density (%f)",
                             layer, temp->soil_dens_org[layer],
                             temp->soil_dens_min[layer]);
                    temp->soil_dens_org[layer] = temp->soil_dens_min[layer];
                }
                if (temp->organic[layer] > 0 && temp->bulk_dens_org[layer] >=
                    temp->soil_dens_org[layer]) {
                    log_err("layer %zu organic bulk density (%f) "
                            "must be less than organic soil density (%f)",
                            layer, temp->bulk_dens_org[layer],
                            temp->soil_dens_org[layer]);
                }
            }
        }
        else {
            for (layer = 0; layer < options.Nlayer; layer++) {
                temp->organic[layer] = 0.0;
                temp->bulk_dens_org[layer] = MISSING;
                temp->soil_dens_org[layer] = MISSING;
            }
        }

        if (options.BULK_DENSITY_COMB) {
            /* read soil bulk density */
            for (layer = 0; layer < options.Nlayer; layer++) {
                token = strtok(NULL, delimiters);
                while (token != NULL && (length = strlen(token)) == 0) {
                    token = strtok(NULL, delimiters);
                }
                if (token == NULL) {
                    log_err("Can't find values for SOIL BULK DENSITY for "
                            "layer %zu in soil file", layer);
                }
                sscanf(token, "%lf", &(temp->bulk_density)[layer]);
            }
        }
        else {
            for (layer = 0; layer < options.Nlayer; layer++) {
                temp->bulk_density[layer] =
                    (1 -
                     temp->organic[layer]) * temp->bulk_dens_min[layer] +
                    temp->organic[layer] * temp->bulk_dens_org[layer];
            }
        }

        /* read cell gmt offset */
        token = strtok(NULL, delimiters);
        while (token != NULL && (length = strlen(token)) == 0) {
            token = strtok(NULL, delimiters);
        }
        if (token == NULL) {
            log_err("Can't find values for GMT OFFSET in soil file");
        }
        sscanf(token, "%lf", &off_gmt);

        /* read layer saturated point */
        for (layer = 0; layer < options.Nlayer; layer++) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for SATURATED POINT for layer %zu "
                        "in soil file", layer);
            }
            sscanf(token, "%lf", &(temp->Wsat[layer]));
        }

        /* read layer field capacity */
        for (layer = 0; layer < options.Nlayer; layer++) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for FIELD CAPACITY for layer %zu "
                        "in soil file", layer);
            }
            sscanf(token, "%lf", &(temp->Wfc[layer]));
        }

        /* read layer wilting point */
        for (layer = 0; layer < options.Nlayer; layer++) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for WILTING POINT for layer %zu "
                        "in soil file", layer);
            }
            sscanf(token, "%lf", &(temp->Wpwp[layer]));
        }

        /* read layer saturated matric potential */
        for (layer = 0; layer < options.Nlayer; layer++) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for SATURATED MATRIC for layer %zu "
                        "in soil file", layer);
            }
            sscanf(token, "%lf", &(temp->psi_sat[layer]));
        }

        /* read layer residual moisture content */
        for (layer = 0; layer < options.Nlayer; layer++) {
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
            if (token == NULL) {
                log_err("Can't find values for RESIDUAL MOISTURE CONTENT for "
                        "layer %zu in soil file", layer);
            }
            sscanf(token, "%lf", &(temp->resid_moist[layer]));
        }

        /* read frozen soil active flag */
        token = strtok(NULL, delimiters);
        while (token != NULL && (length = strlen(token)) == 0) {
            token = strtok(NULL, delimiters);
        }
        if (token == NULL) {
            log_err("Can't find values for FROZEN SOIL ACTIVE FLAG in "
                    "soil file");
        }
        sscanf(token, "%d", &tempint);
        temp->FS_ACTIVE = (char)tempint;

        /*******************************************
           End of soil parameters for this grid cell
        *******************************************/

        /*******************************************
           Compute Soil Layer Properties
        *******************************************/
        for (layer = 0; layer < options.Nlayer; layer++) {
            temp->soil_density[layer] =
                (1 -
                 temp->organic[layer]) * temp->soil_dens_min[layer] +
                temp->organic[layer] * temp->soil_dens_org[layer];
            if (temp->resid_moist[layer] == MISSING) {
                temp->resid_moist[layer] = param.SOIL_RESID_MOIST;
            }
            temp->porosity[layer] = 1.0 - temp->bulk_density[layer] /
                                    temp->soil_density[layer];
            temp->max_moist[layer] = temp->depth[layer] * 
                                     temp->porosity[layer];
        }

        /**********************************************
           Validate Soil Layer Thicknesses
        **********************************************/
        for (layer = 0; layer < options.Nlayer; layer++) {
            if (temp->depth[layer] < MINSOILDEPTH) {
                log_err("Model will not function with layer %zu "
                        "depth %f < %f m.", layer, temp->depth[layer],
                        MINSOILDEPTH);
            }
        }
        if (temp->depth[0] > temp->depth[1]) {
            log_err("Model will not function with layer %d depth"
                    "(%f m) > layer %d depth (%f m).", 0, temp->depth[0],
                    1, temp->depth[1]);
        }

        /**********************************************
           Compute Maximum Infiltration for Upper Layers
        **********************************************/
        if (options.Nlayer == 6) {
            temp->max_infil =
                (1.0 + temp->b_infilt) * (temp->max_moist[0] + 
                                 temp->max_moist[1] + temp->max_moist[2]);
        }
        else {
            log_err("The options.Nlayer (%zu) of this VIC version must be 6. ", 
                        options.Nlayer);
        }

        /****************************************************************
           Compute Soil Layer Critical and Wilting Point Moisture Contents
        ****************************************************************/
        for (layer = 0; layer < options.Nlayer; layer++) {
            if (temp->Wpwp[layer] < temp->resid_moist[layer]) {
                log_err("Calculated wilting point moisture (%f mm) is "
                        "less than calculated residual moisture (%f mm) "
                        "for layer %zu.\n\tIn the soil parameter file, "
                        "Wpwp_FRACT MUST be >= resid_moist / (1.0 - "
                        "bulk_density/soil_density).", temp->Wpwp[layer],
                        temp->resid_moist[layer], layer);
            }
        }

        /*******************************************************************
           Calculate grid cell area.
         ******************************************************************/
        compute_cell_area(temp);

        /*************************************************
           Allocate and Initialize Snow Band Parameters
        *************************************************/
        Nbands = options.SNOW_BAND;
        temp->AreaFract = calloc(Nbands, sizeof(*(temp->AreaFract)));
        check_alloc_status(temp->AreaFract, "Memory allocation error.");
        temp->BandElev = calloc(Nbands, sizeof(*(temp->BandElev)));
        check_alloc_status(temp->BandElev, "Memory allocation error.");
        temp->Tfactor = calloc(Nbands, sizeof(*(temp->Tfactor)));
        check_alloc_status(temp->Tfactor, "Memory allocation error.");
        temp->Pfactor = calloc(Nbands, sizeof(*(temp->Pfactor)));
        check_alloc_status(temp->Pfactor, "Memory allocation error.");

        /** Set default values for factors to use unmodified forcing data **/
        for (band = 0; band < Nbands; band++) {
            temp->AreaFract[band] = 0.;
            temp->BandElev[band] = temp->elevation;
            temp->Tfactor[band] = 0.;
            temp->Pfactor[band] = 1.;
        }
        temp->AreaFract[0] = 1.;
        
        size_t Nsoil = options.Nsoil;
        double dp = 0.; 
        for (layer = 0; layer < options.Nlayer; layer++) {
            dp += temp->depth[layer];
        }

        if (!options.EXP_TRANS) {
            temp->dz_soil[0] = 0.025;
            temp->dz_soil[1] = 0.05;
            temp->dz_soil[2] = 0.05;
            temp->dz_soil[3] = 0.1;
            temp->dz_soil[4] = 0.1;
            temp->dz_soil[5] = 0.1;
            temp->dz_soil[6] = 0.1;
            temp->dz_soil[7] = 0.3;
            temp->dz_soil[8] = 0.3;
            temp->dz_soil[9] = 0.3;
            temp->dz_soil[10] = 0.3;
            temp->dz_soil[11] = 0.275;
            /* Compute soil node depths */
            double sum_dz = 0.;
            for (k = 0; k < Nsoil; k++) {
                sum_dz += temp->dz_soil[k];
                temp->Zsum_soil[k] = sum_dz;
            }
            for (k = 0; k < Nsoil; k++) {
                temp->zc_soil[k] = temp->Zsum_soil[k] - temp->dz_soil[k] / 2.;
            }
        }
        else {
            // exponential grid transformation, EXP_TRANS = TRUE
            // calculate exponential function parameter
            // to force Zsum=dp at bottom node
        } // end if !EXP_TRANS

        /*************************************************
           Compute soil albedo in PAR range (400-700nm)
           following eqn 122 in Knorr 1997
        *************************************************/
        if (options.CARBON) {
            temp->AlbedoPar = 0.92 * param.ALBEDO_BARE_SOIL - 0.015;
        }
        // set soil albedo
        temp->AlbedoSat[0] = 0.09;  // saturated soil albedo at visible band
        temp->AlbedoSat[1] = 0.18;  // saturated soil albedo at NIR band
        temp->AlbedoDry[0] = 0.18;  // dry soil albedo at visible band
        temp->AlbedoDry[1] = 0.36;  // dry soil albedo at NIR band

        /*************************************************
           Miscellaneous terms for MTCLIM disaggregation
        *************************************************/
        /* Central Longitude of Current Time Zone */
        temp->time_zone_lng = off_gmt * 360. / HOURS_PER_DAY;
        /* Assume flat grid cell for radiation calculations */
        temp->slope = 0;
    } // end if(!(*MODEL_DONE) && (*RUN_MODEL))
}

