/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine reads the VIC model parameters file
 *****************************************************************************/

#include "vic_driver_shared_all.h"

/******************************************************************************
 * @brief    Read the VIC model parameters file
 *****************************************************************************/
void
get_parameters(FILE *paramfile)
{
    char                     cmdstr[MAXSTRING];
    char                     optstr[MAXSTRING];

    extern parameters_struct param;

    /** Read through parameter file to find parameters **/

    rewind(paramfile);
    fgets(cmdstr, MAXSTRING, paramfile);

    while (!feof(paramfile)) {
        if (cmdstr[0] != '#' && cmdstr[0] != '\n' && cmdstr[0] != '\0') {
            sscanf(cmdstr, "%s", optstr);

            /* Handle case of comment line in which '#' is indented */
            if (optstr[0] == '#') {
                fgets(cmdstr, MAXSTRING, paramfile);
                continue;
            }

            /*************************************
               Get Model Parameters
            *************************************/
            // Lapse Rate
            if (strcasecmp("LAPSE_RATE", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.LAPSE_RATE);
            }
            // Precipitation Guage Height
            else if (strcasecmp("GAUGE_HEIGHT", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.GAUGE_HEIGHT);
            }
            // Huge Resistance Term
            else if (strcasecmp("HUGE_RESIST", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.HUGE_RESIST);
            }
            // Surface Albedo Parameters
            else if (strcasecmp("ALBEDO_BARE_SOIL", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.ALBEDO_BARE_SOIL);
            }
            // Surface Emissivities
            else if (strcasecmp("EMISS_GRND", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.EMISS_GRND);
            }
            else if (strcasecmp("EMISS_ICE", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.EMISS_ICE);
            }
            else if (strcasecmp("EMISS_VEG", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.EMISS_VEG);
            }
            else if (strcasecmp("EMISS_SNOW", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.EMISS_SNOW);
            }
            else if (strcasecmp("EMISS_H2O", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.EMISS_H2O);
            }
            // Vegetation Parameters
            else if (strcasecmp("VEG_LAI_SNOW_MULTIPLIER", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.VEG_LAI_SNOW_MULTIPLIER);
            }
            else if (strcasecmp("VEG_LAI_WATER_FACTOR", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.VEG_LAI_WATER_FACTOR);
            }
            // Canopy Parameters
            else if (strcasecmp("CANOPY_CLOSURE", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.CANOPY_CLOSURE);
            }
            else if (strcasecmp("CANOPY_RSMAX", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.CANOPY_RSMAX);
            }
            // Saturation Vapor Pressure Parameters
            else if (strcasecmp("SVP_A0", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SVP_A0);
            }
            else if (strcasecmp("SVP_A1", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SVP_A1);
            }
            else if (strcasecmp("SVP_A2", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SVP_A2);
            }
            // Photosynthesis Parameters
            else if (strcasecmp("PHOTO_LRESC3", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_LRESC3);
            }
            else if (strcasecmp("PHOTO_LRESC4", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_LRESC4);
            }
            else if (strcasecmp("PHOTO_OX", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_OX);
            }
            else if (strcasecmp("PHOTO_KC", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_KC);
            }
            else if (strcasecmp("PHOTO_KO", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_KO);
            }
            else if (strcasecmp("PHOTO_EC", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_EC);
            }
            else if (strcasecmp("PHOTO_EO", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_EO);
            }
            else if (strcasecmp("PHOTO_EV", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_EV);
            }
            else if (strcasecmp("PHOTO_ER", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_ER);
            }
            // Snow Parameters
            else if (strcasecmp("SNOW_MAX_SURFACE_SWE", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_MAX_SURFACE_SWE);
            }
            else if (strcasecmp("SNOW_LIQUID_WATER_CAPACITY", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_LIQUID_WATER_CAPACITY);
            }
            else if (strcasecmp("SNOW_NEW_SNOW_DENSITY", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_NEW_SNOW_DENSITY);
            }
            else if (strcasecmp("SNOW_NEW_SNOW_DENS_MAX", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_NEW_SNOW_DENS_MAX);
            }
            else if (strcasecmp("SNOW_NEW_SNOW_DENS_MAX", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_NEW_SNOW_DENS_MAX);
            }
            else if (strcasecmp("SNOW_NEW_SNT_C1", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_NEW_SNT_C1);
            }
            else if (strcasecmp("SNOW_NEW_SNT_C2", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_NEW_SNT_C2);
            }
            else if (strcasecmp("SNOW_NEW_SNT_C3", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_NEW_SNT_C3);
            }
            else if (strcasecmp("SNOW_NEW_BRAS_DENOM", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_NEW_BRAS_DENOM);
            }
            else if (strcasecmp("SNOW_CONDUCT", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_CONDUCT);
            }
            else {
                log_warn("Unrecognized option in the parameter file:  %s "
                         "- check your spelling", optstr);
            }
        }
        fgets(cmdstr, MAXSTRING, paramfile);
    }
}

/******************************************************************************
 * @brief    Validate VIC model parameter values
 *****************************************************************************/
void
validate_parameters()
{
    extern parameters_struct param;

    // Validate Parameters
    // Lapse Rate
    if (!(param.LAPSE_RATE >= -1 && param.LAPSE_RATE <= 0)) {
        log_err("LAPSE_RATE must be defined on the interval [-1,0] (C/m)")
    }
    // Precipitation Guage Height
    if (!(param.GAUGE_HEIGHT >= 0 && param.GAUGE_HEIGHT <= 100)) {
        log_err("GAUGE_HEIGHT must be defined on the interval [0,100] (m)")
    }
    // Huge Resistance Term
    if (!(param.HUGE_RESIST >= 0.)) {
        log_err("HUGE_RESIST must be defined on the interval [0, inf) (s/m)");
    }
    // Surface Albedo Parameters
    if (!(param.ALBEDO_BARE_SOIL >= 0 && param.ALBEDO_BARE_SOIL <= 1)) {
        log_err("ALBEDO_BARE_SOIL must be defined on the interval [0,1] (-)")
    }
    // Surface Emissivities
    if (!(param.EMISS_GRND >= 0 && param.EMISS_GRND <= 1)) {
        log_err("EMISS_GRND must be defined on the interval [0,1] (-)")
    }
    if (!(param.EMISS_ICE >= 0 && param.EMISS_ICE <= 1)) {
        log_err("EMISS_ICE must be defined on the interval [0,1] (-)")
    }
    if (!(param.EMISS_VEG >= 0 && param.EMISS_VEG <= 1)) {
        log_err("EMISS_VEG must be defined on the interval [0,1] (-)")
    }
    if (!(param.EMISS_SNOW >= 0 && param.EMISS_SNOW <= 1)) {
        log_err("EMISS_SNOW must be defined on the interval [0,1] (-)")
    }
    if (!(param.EMISS_H2O >= 0 && param.EMISS_H2O <= 1)) {
        log_err("EMISS_H2O must be defined on the interval [0,1] (-)")
    }
    // Vegetation Parameters
    if (!(param.VEG_LAI_SNOW_MULTIPLIER >= 0.)) {
        log_err(
            "VEG_LAI_SNOW_MULTIPLIER must be defined on the interval [0, inf) (-)");
    }
    if (!(param.VEG_LAI_WATER_FACTOR >= 0.)) {
        log_err(
            "VEG_LAI_WATER_FACTOR must be defined on the interval [0, inf) (-)");
    }
    // Canopy Parameters
    if (!(param.CANOPY_CLOSURE >= 0.)) {
        log_err("CANOPY_CLOSURE must be defined on the interval [0, inf) (Pa)");
    }
    if (!(param.CANOPY_RSMAX >= 0.)) {
        log_err("CANOPY_RSMAX must be defined on the interval [0, inf) (s/m)");
    }
    // Photosynthesis Parameters
    // PHOTO_OMEGA - Currently, no constraints
    // PHOTO_FCI1C3 - Currently, no constraints
    // PHOTO_FCI1C4 - Currently, no constraints
    if (!(param.PHOTO_OX >= 0.)) {
        log_err(
            "PHOTO_OX must be defined on the interval [0, inf) (mol H2O/m2s)");
    }
    // PHOTO_KC - Currently, no constraints
    // PHOTO_KO - Currently, no constraints
    // PHOTO_EC - Currently, no constraints
    // PHOTO_EO - Currently, no constraints
    // PHOTO_EV - Currently, no constraints
    // PHOTO_ER - Currently, no constraints
    // PHOTO_ALC3 - Currently, no constraints
    // PHOTO_FRDC3 - Currently, no constraints
    // PHOTO_EK - Currently, no constraints
    // PHOTO_ALC4 - Currently, no constraints
    // PHOTO_FRDC4 - Currently, no constraints
    // PHOTO_THETA - Currently, no constraints
    // PHOTO_FRLEAF - Currently, no constraints
    // PHOTO_FRGROWTH - Currently, no constraints

    // Snow Parameters
    if (!(param.SNOW_MAX_SURFACE_SWE >= 0.)) {
        log_err(
            "SNOW_MAX_SURFACE_SWE must be defined on the interval [0, inf) (m)");
    }
    if (!(param.SNOW_LIQUID_WATER_CAPACITY >= 0 &&
          param.SNOW_LIQUID_WATER_CAPACITY <= 1)) {
        log_err(
            "SNOW_LIQUID_WATER_CAPACITY must be defined on the interval [0,1] (-)")
    }
    if (!(param.SNOW_NEW_SNOW_DENSITY >= 0.)) {
        log_err(
            "SNOW_NEW_SNOW_DENSITY must be defined on the interval [0, inf) (kg/m^3)");
    }
    if (!(param.SNOW_NEW_SNOW_DENS_MAX >= 0. &&
          param.SNOW_NEW_SNOW_DENS_MAX <= 700.)) {
        log_err(
            "SNOW_NEW_SNOW_DENS_MAX must be defined on the interval [0, 700) (kg/m^3)");
    }
    // SNOW_DENS_ETA0 - Currently, no constraints
    // SNOW_DENS_C1 - Currently, no constraints
    // SNOW_DENS_C2 - Currently, no constraints
    // SNOW_DENS_C5 - Currently, no constraints
    // SNOW_DENS_C6 - Currently, no constraints
    // SNOW_DENS_F - Currently, no constraints
    if (!(param.SNOW_CONDUCT >= 0.)) {
        log_err("SNOW_CONDUCT must be defined on the interval [0, inf) (W/mK)");
    }

    // Convergence Tolerances

}