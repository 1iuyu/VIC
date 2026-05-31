/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine reads the VIC model global control file, getting values for
 * global parameters, model options, and debugging controls.
 *****************************************************************************/

#include "vic_driver_image.h"

/******************************************************************************
 * @brief    Read the VIC model global control file, getting values for
 *           global parameters, model options, and debugging controls.
 *****************************************************************************/
void
get_global_param(FILE *gp)
{
    extern option_struct       options;
    extern global_param_struct global_param;
    extern param_set_struct    param_set;
    extern filenames_struct    filenames;
    extern size_t              NF, NR;

    char                       cmdstr[MAXSTRING];
    char                       optstr[MAXSTRING];
    char                       flgstr[MAXSTRING];
    char                       flgstr2[MAXSTRING];
    size_t                     file_num;
    int                        field;
    int                        status;
    unsigned int               tmpstartdate;
    unsigned int               tmpenddate;
    unsigned short int         lastday[MONTHS_PER_YEAR];


    /** Read through global control file to find parameters **/

    rewind(gp);
    fgets(cmdstr, MAXSTRING, gp);

    while (!feof(gp)) {
        if (cmdstr[0] != '#' && cmdstr[0] != '\n' && cmdstr[0] != '\0') {
            sscanf(cmdstr, "%s", optstr);

            /* Handle case of comment line in which '#' is indented */
            if (optstr[0] == '#') {
                fgets(cmdstr, MAXSTRING, gp);
                continue;
            }

            /*************************************
               Get Model Global Parameters
            *************************************/
            if (strcasecmp("MODEL_STEPS_PER_DAY", optstr) == 0) {
                sscanf(cmdstr, "%*s %zu", &global_param.model_steps_per_day);
            }
            else if (strcasecmp("STARTYEAR", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.startyear);
            }
            else if (strcasecmp("STARTMONTH", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.startmonth);
            }
            else if (strcasecmp("STARTDAY", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.startday);
            }
            else if (strcasecmp("STARTSEC", optstr) == 0) {
                sscanf(cmdstr, "%*s %u", &global_param.startsec);
            }
            else if (strcasecmp("NRECS", optstr) == 0) {
                sscanf(cmdstr, "%*s %zu", &global_param.nrecs);
            }
            else if (strcasecmp("ENDYEAR", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.endyear);
            }
            else if (strcasecmp("ENDMONTH", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.endmonth);
            }
            else if (strcasecmp("ENDDAY", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.endday);
            }
            else if (strcasecmp("CALENDAR", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                global_param.calendar = str_to_calendar(flgstr);
            }
            else if (strcasecmp("OUT_TIME_UNITS", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                global_param.time_units = str_to_timeunits(flgstr);
            }
            else if (strcasecmp("FROZEN_SOIL", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.FROZEN_SOIL = str_to_bool(flgstr);
            }
            else if ((strcasecmp("NOFLUX",
                                 optstr) == 0) ||
                     (strcasecmp("NO_FLUX", optstr) == 0)) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.NOFLUX = str_to_bool(flgstr);
            }
            else if (strcasecmp("SNOW_DENSITY", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("DENS_SNTHRM", flgstr) == 0) {
                    options.SNOW_DENSITY = DENS_SNTHRM;
                }
                else if (strcasecmp("DENS_BRAS", flgstr) == 0) {
                    options.SNOW_DENSITY = DENS_BRAS;
                }
                else {
                    log_err("Unknown SNOW_DENSITY option: %s", flgstr);
                }
            }
            else if (strcasecmp("CORRPREC", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.CORRPREC = str_to_bool(flgstr);
            }
            else if (strcasecmp("RESOLUTION", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &global_param.resolution);
            }
            else if (strcasecmp("AERO_RESIST", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("AR_ZENG", flgstr) == 0) {
                    options.AERO_RESIST = AR_ZENG;
                }
                else if (strcasecmp("AR_MEIER", flgstr) == 0) {
                    options.AERO_RESIST = AR_MEIER;
                }
                else {
                    log_err("Unknown AERO_RESIST option: %s", flgstr);
                }
            }
            else if (strcasecmp("TFALLBACK", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.TFALLBACK = str_to_bool(flgstr);
            }
            else if (strcasecmp("CANOPY_LAYERS", optstr) == 0) {
                sscanf(cmdstr, "%*s %zu", &options.Ncanopy);
            }
            else if (strcasecmp("CARBON", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.CARBON = str_to_bool(flgstr);
            }

            /*************************************
               Define log directory
            *************************************/
            else if (strcasecmp("LOG_DIR", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", filenames.log_path);
            }

            /*************************************
               Define state files
            *************************************/
            else if (strcasecmp("INIT_STATE", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("FALSE", flgstr) == 0) {
                    options.INIT_STATE = false;
                }
                else {
                    options.INIT_STATE = true;
                    strcpy(filenames.init_state.nc_filename, flgstr);
                }
            }
            else if (strcasecmp("STATENAME", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", filenames.statefile);
                options.SAVE_STATE = true;
            }
            else if (strcasecmp("STATEYEAR", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.stateyear);
            }
            else if (strcasecmp("STATEMONTH", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.statemonth);
            }
            else if (strcasecmp("STATEDAY", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.stateday);
            }
            else if (strcasecmp("STATESEC", optstr) == 0) {
                sscanf(cmdstr, "%*s %u", &global_param.statesec);
            }
            else if (strcasecmp("STATE_FORMAT", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("NETCDF3_CLASSIC", flgstr) == 0) {
                    options.STATE_FORMAT = NETCDF3_CLASSIC;
                }
                else if (strcasecmp("NETCDF3_64BIT_OFFSET", flgstr) == 0) {
                    options.STATE_FORMAT = NETCDF3_64BIT_OFFSET;
                }
                else if (strcasecmp("NETCDF4_CLASSIC", flgstr) == 0) {
                    options.STATE_FORMAT = NETCDF4_CLASSIC;
                }
                else if (strcasecmp("NETCDF4", flgstr) == 0) {
                    options.STATE_FORMAT = NETCDF4;
                }
                else {
                    log_err("STATE_FORMAT must be either NETCDF3_CLASSIC, "
                            "NETCDF3_64BIT_OFFSET, NETCDF4_CLASSIC, or NETCDF4.");
                }
            }

            /*************************************
               Define forcing files
            *************************************/
            else if (strcasecmp("FORCING1", optstr) == 0) {
                if (strcmp(filenames.f_path_pfx[0], "MISSING") != 0) {
                    log_err("Tried to define FORCING1 twice, if you want to "
                            "use two forcing files, the second must be "
                            "defined as FORCING2");
                }
                sscanf(cmdstr, "%*s %s", filenames.f_path_pfx[0]);
                file_num = 0;
                field = 0;
                // count the number of forcing variables in this file
                param_set.N_TYPES[file_num] = count_force_vars(gp);
            }
            else if (strcasecmp("FORCING2", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", filenames.f_path_pfx[1]);
                if (strcasecmp("FALSE", filenames.f_path_pfx[1]) == 0) {
                    strcpy(filenames.f_path_pfx[1], "MISSING");
                }
                file_num = 1;
                field = 0;
                // count the number of forcing variables in this file
                param_set.N_TYPES[file_num] = count_force_vars(gp);
            }
            else if (strcasecmp("FORCE_TYPE", optstr) == 0) {
                set_force_type(cmdstr, file_num, &field);
            }
            else if (strcasecmp("FORCE_STEPS_PER_DAY", optstr) == 0) {
                sscanf(cmdstr, "%*s %zu",
                       &param_set.force_steps_per_day[file_num]);
            }
            else if (strcasecmp("FORCEYEAR", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu",
                       &global_param.forceyear[file_num]);
            }
            else if (strcasecmp("FORCEMONTH", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu",
                       &global_param.forcemonth[file_num]);
            }
            else if (strcasecmp("FORCEDAY", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.forceday[file_num]);
            }
            else if (strcasecmp("FORCESEC", optstr) == 0) {
                sscanf(cmdstr, "%*s %u", &global_param.forcesec[file_num]);
            }

            /*************************************
               Define parameter files
            *************************************/
            else if (strcasecmp("CONSTANTS", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", filenames.constants);
            }
            else if (strcasecmp("DOMAIN", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", filenames.domain.nc_filename);
            }
            else if (strcasecmp("DOMAIN_TYPE", optstr) == 0) {
                get_domain_type(cmdstr);
            }
            else if (strcasecmp("PARAMETERS", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", filenames.params.nc_filename);
            }
            else if (strcasecmp("ROUT_PARAM", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("FALSE", flgstr) == 0) {
                    options.ROUT = false;
                }
                else {
                    options.ROUT = true;
                    strcpy(filenames.rout_params.nc_filename, flgstr);
                }
            }
            else if (strcasecmp("DENSITY_FROM_SOIL", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.PARAM_FROM_SOIL = str_to_bool(flgstr);
            }
            else if (strcasecmp("VEGLIB", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", filenames.veglib);
            }
            else if (strcasecmp("VEGLIB_FCAN", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.VEGLIB_FCAN = str_to_bool(flgstr);
            }
            else if (strcasecmp("VEGPARAM_LAI", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.VEGPARAM_LAI = str_to_bool(flgstr);
            }
            else if (strcasecmp("LAI_SRC", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("FROM_VEGHIST", flgstr) == 0) {
                    options.LAI_SRC = FROM_VEGHIST;
                }
                else if (strcasecmp("FROM_VEGPARAM", flgstr) == 0) {
                    options.LAI_SRC = FROM_VEGPARAM;
                }
                else if (strcasecmp("FROM_VEGLIB", flgstr) == 0) {
                    options.LAI_SRC = FROM_VEGLIB;
                }
                else {
                    log_err("Unrecognized value of LAI_SRC in the global "
                            "control file.");
                }
            }
            else if (strcasecmp("FCAN_SRC", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("FROM_VEGHIST", flgstr) == 0) {
                    options.FCAN_SRC = FROM_VEGHIST;
                }
                else if (strcasecmp("FROM_VEGPARAM", flgstr) == 0) {
                    options.FCAN_SRC = FROM_VEGPARAM;
                }
                else if (strcasecmp("FROM_VEGLIB", flgstr) == 0) {
                    options.FCAN_SRC = FROM_VEGLIB;
                }
                else if (strcasecmp("FROM_DEFAULT", flgstr) == 0) {
                    options.FCAN_SRC = FROM_DEFAULT;
                }
                else {
                    log_err("Unrecognized value of FCAN_SRC in the global "
                            "control file.");
                }
            }
            else if (strcasecmp("SNOW_BAND", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (str_to_bool(flgstr)) {
                    options.SNOW_BAND = SNOW_BAND_TRUE_BUT_UNSET;
                }
            }
            else if (strcasecmp("GLACIER_ID", optstr) == 0) {
                sscanf(cmdstr,"%*s %d", &options.GLACIER_ID);
            }
            /*************************************
               Define output files
            *************************************/
            else if (strcasecmp("RESULT_DIR", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", filenames.result_dir);
            }

            /*************************************
               Define output file contents
            *************************************/
            else if (strcasecmp("OUTFILE", optstr) == 0) {
                ; // do nothing
            }
            else if (strcasecmp("OUTVAR", optstr) == 0) {
                ; // do nothing
            }
            else if (strcasecmp("AGGFREQ", optstr) == 0) {
                ; // do nothing
            }
            else if (strcasecmp("HISTFREQ", optstr) == 0) {
                ; // do nothing
            }
            else if (strcasecmp("COMPRESS", optstr) == 0) {
                ; // do nothing
            }
            else if (strcasecmp("OUT_FORMAT", optstr) == 0) {
                ; // do nothing
            }

            /*************************************
               Fail when classic driver specific options are used
            *************************************/
            else if (strcasecmp("OUTPUT_FORCE", optstr) == 0) {
                log_err("OUTPUT_FORCE is not a valid option for this driver.  "
                        "Update your global parameter file accordingly.");
            }

            /***********************************
               Unrecognized Global Parameter Flag
            ***********************************/
            else {
                log_warn("Unrecognized option in the global parameter file: %s"
                         "\n - check your spelling", optstr);
            }
        }
        fgets(cmdstr, MAXSTRING, gp);
    }

    /******************************************
       Check for undefined required parameters
    ******************************************/

    // Validate model time step
    if (global_param.model_steps_per_day == 0) {
        log_err("Model time steps per day has not been defined.  Make sure "
                "that the global file defines MODEL_STEPS_PER_DAY.");
    }
    else if (global_param.model_steps_per_day != 1 &&
             global_param.model_steps_per_day <
             MIN_SUBDAILY_STEPS_PER_DAY) {
        log_err("The specified number of model steps per day (%zu) > 1 and < "
                "the minimum number of subdaily steps per day (%d).  Make "
                "sure that the global file defines MODEL_STEPS_PER_DAY of at "
                "least (%d).", global_param.model_steps_per_day,
                MIN_SUBDAILY_STEPS_PER_DAY,
                MIN_SUBDAILY_STEPS_PER_DAY);
    }
    else if (global_param.model_steps_per_day >
             MAX_SUBDAILY_STEPS_PER_DAY) {
        log_err("The specified number of model steps per day (%zu) > the "
                "the maximum number of subdaily steps per day (%d).  Make "
                "sure that the global file defines MODEL_STEPS_PER_DAY of at "
                "most (%d).", global_param.model_steps_per_day,
                MAX_SUBDAILY_STEPS_PER_DAY,
                MAX_SUBDAILY_STEPS_PER_DAY);
    }
    else if ((global_param.model_steps_per_day > HOURS_PER_DAY) &&
             (global_param.model_steps_per_day % HOURS_PER_DAY) != 0) {
        log_err("The specified number of model steps per day (%zu) is > 24 "
                "and is not evenly divided by 24.",
                global_param.model_steps_per_day);
    }
    else {
        global_param.step_dt = SEC_PER_DAY /
                          (double) global_param.model_steps_per_day;
    }

    // set NR and NF
    NF = param_set.force_steps_per_day[0] / global_param.model_steps_per_day;
    if (NF == 1) {
        NR = 0;
    }
    else {
        NR = NF;
    }

    // Validate simulation start date
    if (global_param.startyear == 0) {
        log_err("Simulation start year has not been defined.  Make sure that "
                "the global file defines STARTYEAR.");
    }
    if (global_param.startmonth == 0) {
        log_err("Simulation start month has not been defined.  Make sure that "
                "the global file defines STARTMONTH.");
    }
    else if (global_param.startmonth > MONTHS_PER_YEAR) {
        log_err("The specified simulation start month (%hu) > 12. Make "
                "sure that the global file defines a positive integer for "
                "STARTMONTH.", global_param.startmonth);
    }
    if (global_param.startday == 0) {
        log_err("Simulation start day has not been defined.  Make sure that "
                "the global file defines STARTDAY.");
    }
    if (global_param.model_steps_per_day == 1) {
        global_param.startsec = 0;
    }
    else if (global_param.startsec > SEC_PER_DAY) {
        log_err("The specified simulation start second (%u) > 86400  Make sure "
                "that the global file defines a positive integer "
                "for STARTSEC.",
                global_param.startsec);
    }

    // Validate simulation end date and/or number of timesteps
    make_lastday(global_param.calendar, global_param.endyear, lastday);

    if (global_param.nrecs == 0 && global_param.endyear == 0 &&
        global_param.endmonth == 0 && global_param.endday == 0) {
        log_err("The model global file MUST define EITHER the number of "
                "records to simulate (NRECS), or the year (ENDYEAR), month "
                "(ENDMONTH), and day (ENDDAY) of the last full simulation day");
    }
    else if (global_param.nrecs == 0) {
        if (global_param.endyear == 0) {
            log_err("Simulation end year has not been defined.  Make sure "
                    "that the global file defines ENDYEAR.");
        }
        if (global_param.endmonth == 0) {
            log_err("Simulation end month has not been defined.  Make sure "
                    "that the global file defines ENDMONTH.");
        }
        else if (global_param.endmonth > MONTHS_PER_YEAR) {
            log_err("The specified simulation end month (%hu) < 0.  Make sure "
                    "that the global file defines a positive integer for "
                    "ENDMONTH.", global_param.endmonth);
        }
        if (global_param.endday == 0) {
            log_err("Simulation end day has not been defined.  Make sure "
                    "that the global file defines ENDDAY.");
        }
        else if (global_param.endday > lastday[global_param.endmonth - 1]) {
            log_err("The specified simulation end day (%hu) > the number of "
                    "days in the ENDMONTH (%hu).  Make sure that the global "
                    "file defines a positive integer for ENDDAY.",
                    global_param.endday, global_param.endmonth);
        }
        tmpstartdate = global_param.startyear * 10000 +
                       global_param.startmonth * 100 +
                       global_param.startday;
        tmpenddate = global_param.endyear * 10000 +
                     global_param.endmonth * 100 +
                     global_param.endday;
        if (tmpenddate < tmpstartdate) {
            log_err("The specified simulation end date (%04d-%02d-%02d) is "
                    "EARLIER than the specified start date (%04d-%02d-%02d).",
                    global_param.endyear, global_param.endmonth,
                    global_param.endday,
                    global_param.startyear, global_param.startmonth,
                    global_param.startday);
        }
    }
    else if (global_param.nrecs < 1) {
        log_err("The specified duration of simulation (%zu) < 1 time step. "
                "Make sure that the global file defines a positive integer "
                "for NRECS.", global_param.nrecs);
    }

    // Validate forcing files and variables
    if (strcmp(filenames.f_path_pfx[0], "MISSING") == 0) {
        log_err("No forcing file has been defined.  Make sure that the global "
                "file defines FORCING1.");
    }

    // Get information from the forcing file(s)
    // Open first-year forcing files and get info
    snprintf(filenames.forcing[0].nc_filename, 
             sizeof(filenames.forcing[0].nc_filename), "%s%04d.nc",
             filenames.f_path_pfx[0], global_param.startyear);
    status = nc_open(filenames.forcing[0].nc_filename, NC_NOWRITE,
                     &(filenames.forcing[0].nc_id));
    check_nc_status(status, "Error opening %s",
                    filenames.forcing[0].nc_filename);
    get_forcing_file_info(&param_set, 0);
    if (param_set.N_TYPES[1] != 0) {
        sprintf(filenames.forcing[1].nc_filename, "%s%04d.nc",
                filenames.f_path_pfx[1], global_param.startyear);
        status = nc_open(filenames.forcing[1].nc_filename, NC_NOWRITE,
                         &(filenames.forcing[1].nc_id));
        check_nc_status(status, "Error opening %s",
                        filenames.forcing[1].nc_filename);
        get_forcing_file_info(&param_set, 1);
    }

    if (param_set.N_TYPES[1] != 0 && global_param.forceyear[1] == 0) {
        global_param.forceyear[1] = global_param.forceyear[0];
        global_param.forcemonth[1] = global_param.forcemonth[0];
        global_param.forceday[1] = global_param.forceday[0];
        global_param.forcesec[1] = global_param.forcesec[0];
        global_param.forceskip[1] = 0;
        global_param.forceoffset[1] = global_param.forceskip[1];
    }
    if (param_set.force_steps_per_day[0] == 0) {
        log_err("Forcing file time steps per day has not been "
                "defined.  Make sure that the global file defines "
                "FORCE_STEPS_PER_DAY.");
    }
    else {
        param_set.FORCE_DT[0] = SEC_PER_DAY /
                                (double) param_set.force_steps_per_day[0];
    }
    if (param_set.force_steps_per_day[1] > 0) {
        param_set.FORCE_DT[1] = SEC_PER_DAY /
                                (double) param_set.force_steps_per_day[1];
    }
    else {
        param_set.FORCE_DT[1] = param_set.FORCE_DT[0];
    }

    // Validate result directory
    if (strcmp(filenames.result_dir, "MISSING") == 0) {
        log_err("No results directory has been defined.  Make sure that the "
                "global file defines the result directory on the line that "
                "begins with \"RESULT_DIR\".");
    }

    // Validate parameter file information
    if (strcmp(filenames.params.nc_filename, "MISSING") == 0) {
        log_err("A parameters file has not been defined.  Make sure that the "
                "global file defines the parameters parameter file on the line "
                "that begins with \"PARAMETERS\".");
    }

    // Validate the input state file information
    if (options.INIT_STATE) {
        if (strcmp(filenames.init_state.nc_filename, "MISSING") == 0) {
            log_err("\"INIT_STATE\" was specified, but no input state file "
                    "has been defined.  Make sure that the global file "
                    "defines the inputstate file on the line that begins "
                    "with \"INIT_STATE\".");
        }
    }

    // Validate the output state file information
    if (options.SAVE_STATE) {
        if (strcmp(filenames.statefile, "MISSING") == 0) {
            log_err("\"SAVE_STATE\" was specified, but no output state "
                    "file has been defined.  Make sure that the global "
                    "file defines the output state file on the line that "
                    "begins with \"SAVE_STATE\".");
        }
        if (global_param.stateyear == 0 || global_param.statemonth == 0 ||
            global_param.stateday == 0) {
            log_err("Incomplete specification of the date to save state "
                    "for state file (%s).\nSpecified date (yyyy-mm-dd-sssss): "
                    "%04d-%02d-%02d-%05u\nMake sure STATEYEAR, STATEMONTH, "
                    "and STATEDAY are set correctly in your global parameter "
                    "file.", filenames.statefile, global_param.stateyear,
                    global_param.statemonth, global_param.stateday,
                    global_param.statesec);
        }
        // Check for month, day in range
        make_lastday(global_param.calendar, global_param.stateyear,
                     lastday);
        if (global_param.stateday > lastday[global_param.statemonth - 1] ||
            global_param.statemonth < 1 ||
            global_param.statemonth > MONTHS_PER_YEAR ||
            global_param.stateday < 1 || global_param.stateday > 31 ||
            global_param.statesec > SEC_PER_DAY) {
            log_err("Unusual specification of the date to save state "
                    "for state file (%s).\nSpecified date (yyyy-mm-dd-sssss): "
                    "%04d-%02d-%02d-%05u\nMake sure STATEYEAR, STATEMONTH, "
                    "STATEDAY and STATESEC are set correctly in your global "
                    "parameter file.", filenames.statefile,
                    global_param.stateyear, global_param.statemonth,
                    global_param.stateday, global_param.statesec);
        }
    }
    // Set the statename here temporarily to compare with INIT_STATE name
    if (options.SAVE_STATE) {
        snprintf(flgstr2, sizeof(flgstr2),
                "%s.%04i%02i%02i_%05u.nc",
                filenames.statefile, global_param.stateyear,
                global_param.statemonth, global_param.stateday,
                global_param.statesec);
    }
    if (options.INIT_STATE && options.SAVE_STATE &&
        (strcmp(filenames.init_state.nc_filename, flgstr2) == 0)) {
        log_err("The save state file (%s) has the same name as the "
                "initialize state file (%s).  The initialize state file "
                "will be destroyed when the save state file is opened.",
                filenames.statefile, filenames.init_state.nc_filename);
    }

    // Validate soil parameter/simulation mode combinations
    if (options.Nlayer > MAX_LAYERS) {
        log_err("Global file wants more soil moisture layers (%zu) than "
                "are defined by MAX_LAYERS (%d).  Edit vic_run/include/vic_def.h "
                "and recompile.", options.Nlayer,
                MAX_LAYERS);
    }
    // Default file formats (if unset)
    if (options.SAVE_STATE && options.STATE_FORMAT == UNSET_FILE_FORMAT) {
        options.STATE_FORMAT = NETCDF4_CLASSIC;
    }


    /*********************************
       Output major options
    *********************************/
    display_current_settings(DISP_VERSION);
}
