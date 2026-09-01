/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine reads the VIC model global control file, getting information
 * for output variables list (if any).
 *****************************************************************************/

#include "vic_driver_shared_image.h"

/******************************************************************************
 * @brief    Get output info from global parameter file.
 *****************************************************************************/
void
parse_output_info(FILE           *gp,
                  stream_struct **streams,
                  size_t         *nstreams,
                  dmy_struct     *dmy_current)
{
    extern domain_struct local_domain;

    char                 cmdstr[MAXSTRING];
    char                 optstr[MAXSTRING];
    char                 flgstr[MAXSTRING];
    int                  streamnum;
    char                 varname[MAXSTRING];
    int                  outvarnum;
    char                 freq_type_str[MAXSTRING];
    char                 freq_value_str[MAXSTRING];
    char                 format[MAXSTRING];
    char                 typestr[MAXSTRING];
    int                  type;
    char                 multstr[MAXSTRING];
    char                 aggstr[MAXSTRING];
    char                 domainstr[MAXSTRING];
    char                 nvarsstr[MAXSTRING];
    double               mult;
    unsigned short int   freq;
    int                  freq_n;
    dmy_struct           freq_dmy;
    unsigned short int   agg_type;
    unsigned short int   domain;
    int                  found;
    size_t              *ivar = NULL;

    // Initialize
    streamnum = -1;
    *nstreams = 0;

    ivar = malloc(local_domain.ncells_active * sizeof(*ivar));
    check_alloc_status(ivar, "Memory allocation error.");

    // 构造一个hru数组，因为alloc_out_data定义在shared_all.h下，无法传递local_domain.locations
    for (size_t i = 0; i < local_domain.ncells_active; i++) {
        ivar[i] = local_domain.locations[i].nveg;
    }

    // rewind the global parameter file to the begining and parse only the
    // output file info.
    rewind(gp);
    fgets(cmdstr, MAXSTRING, gp);

    /** Read through global control file to find output info **/
    while (fgets(cmdstr, MAXSTRING, gp) != NULL) {
        if (cmdstr[0] != '#' && cmdstr[0] != '\n' && cmdstr[0] != '\0') {
            sscanf(cmdstr, "%s", optstr);

            if (strcasecmp("OUTFILE", optstr) == 0) {
                if (streamnum >= 0 &&
                    outvarnum != (int)(*streams)[streamnum].nvars) {
                    log_err("Output stream %s specifies %zu variables, but %d OUTVAR entries were found.", 
                            (*streams)[streamnum].prefix, (*streams)[streamnum].nvars, outvarnum);
                }
                streamnum++;
                if (streamnum >= MAX_OUTPUT_STREAMS) {
                    log_err("Found too many output files, was expecting "
                            "%d but found %d", MAX_OUTPUT_STREAMS,
                            streamnum);
                }

                found = sscanf(cmdstr, "%*s %s %s", (*streams)[streamnum].prefix, nvarsstr);

                if (found != 2) { 
                    log_err("Invalid specification for OUTFILE. " 
                            "Expected: OUTFILE <prefix> <nvars>"); 
                }

                size_t nvars = (size_t) atoi(nvarsstr);

                if (nvars == 0) { 
                    log_err("Number of output variables for OUTFILE %s " 
                            "must be greater than zero.", (*streams)[streamnum].prefix); 
                }

                setup_stream(&(*streams)[streamnum], nvars, local_domain.ncells_active, ivar);

                // set default file format
                (*streams)[streamnum].file_format = NETCDF4_CLASSIC;

                outvarnum = 0;
                (*nstreams)++;

            }
            else if (strcasecmp("AGGFREQ", optstr) == 0) {
                if (streamnum < 0) {
                    log_err("Error in global param file: \"OUTFILE\" must be "
                            "specified before you can specify \"AGGFREQ\".");
                }
                found = sscanf(cmdstr, "%*s %s %s", freq_type_str,
                               freq_value_str);

                if (!found) {
                    log_err("No arguments found after AGGFREQ");
                }
                // parse the frequency string to an enum value
                freq = str_to_freq_flag(freq_type_str);

                if (freq == FREQ_DATE) {
                    // Make sure we have a datestring
                    if (found != 2) {
                        log_err(
                            "AGGFREQ was set to DATE but no date string was found");
                    }
                    // parse date from freq_value_str
                    strpdmy(freq_value_str, "%Y-%m-%d", &freq_dmy);
                    // set the alarm
                    set_alarm(dmy_current, freq, &freq_dmy,
                              (&(*streams)[streamnum].agg_alarm));
                }
                else {
                    if (found != 2) {
                        // Default frequency is 1
                        freq_n = 1;
                    }
                    else {
                        // get the frequency value as an integer
                        freq_n = atoi(freq_value_str);
                    }
                    // set the alarm
                    set_alarm(dmy_current, freq, &freq_n,
                              (&(*streams)[streamnum].agg_alarm));
                }
            }
            else if (strcasecmp("HISTFREQ", optstr) == 0) {
                if (streamnum < 0) {
                    log_err("Error in global param file: \"OUTFILE\" must be "
                            "specified before you can specify \"HISTFREQ\".");
                }
                found = sscanf(cmdstr, "%*s %s %s", freq_type_str,
                               freq_value_str);

                if (!found) {
                    log_err("No arguments found after HISTFREQ");
                }
                // parse the frequency string to an enum value
                freq = str_to_freq_flag(freq_type_str);

                if (freq == FREQ_DATE) {
                    // Make sure we have a datestring
                    if (found != 2) {
                        log_err(
                            "HISTFREQ was set to DATE but no date string was found");
                    }
                    // parse date from freq_value_str
                    strpdmy(freq_value_str, "%Y-%m-%d", &freq_dmy);
                    // set the alarm
                    set_alarm(dmy_current, freq, &freq_dmy,
                              (&(*streams)[streamnum].write_alarm));
                }
                else {
                    if (found != 2) {
                        // Default frequency is 1
                        freq_n = 1;
                    }
                    else {
                        // get the frequency value as an integer
                        freq_n = atoi(freq_value_str);
                    }
                    // set the alarm
                    set_alarm(dmy_current, freq, &freq_n,
                              (&(*streams)[streamnum].write_alarm));
                }
            }
            else if (strcasecmp("COMPRESS", optstr) == 0) {
                if (streamnum < 0) {
                    log_err("Error in global param file: \"OUTFILE\" must be "
                            "specified before you can specify \"COMPRESS\".");
                }
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    (*streams)[streamnum].compress =
                        COMPRESSION_LVL_DEFAULT;
                }
                else if (strcasecmp("FALSE", flgstr) == 0) {
                    (*streams)[streamnum].compress = 0;
                }
                else {
                    (*streams)[streamnum].compress = atoi(flgstr);
                }
            }
            else if (strcasecmp("OUT_FORMAT", optstr) == 0) {
                if (streamnum < 0) {
                    log_err("Error in global param file: \"OUTFILE\" must be "
                            "specified before you can specify \"OUT_FORMAT\".");
                }
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("NETCDF3_CLASSIC", flgstr) == 0) {
                    (*streams)[streamnum].file_format = NETCDF3_CLASSIC;
                }
                else if (strcasecmp("NETCDF3_64BIT_OFFSET", flgstr) == 0) {
                    (*streams)[streamnum].file_format =
                        NETCDF3_64BIT_OFFSET;
                }
                else if (strcasecmp("NETCDF4_CLASSIC", flgstr) == 0) {
                    (*streams)[streamnum].file_format = NETCDF4_CLASSIC;
                }
                else if (strcasecmp("NETCDF4", flgstr) == 0) {
                    (*streams)[streamnum].file_format = NETCDF4;
                }
                else {
                    log_err(
                        "Image driver file format must be a valid NETCDF format");
                }
            }
            else if (strcasecmp("OUTVAR", optstr) == 0) {
                if (streamnum < 0) {
                    log_err("Error in global param file: \"OUTFILE\" must be "
                            "specified before you can specify \"OUTVAR\".");
                }

                if (outvarnum >= (int)(*streams)[streamnum].nvars) { 
                    log_err("Too many OUTVAR entries for output stream %s. Expected %zu variables.", 
                            (*streams)[streamnum].prefix, (*streams)[streamnum].nvars); 
                }
                // parse outvar options
                strcpy(varname, "");
                strcpy(format, "");
                strcpy(typestr, "");
                strcpy(multstr, "");
                strcpy(aggstr, "");
                strcpy(domainstr, "");
                found = sscanf(cmdstr, "%*s %s %s %s %s %s %s", varname,
                               format, typestr, multstr, aggstr, domainstr);
                if (!found) {
                    log_err("OUTVAR specified but no variable was listed");
                }
                // interpret string options, set defaults if necessary
                str_to_ascii_format(format);
                agg_type = str_to_agg_type(aggstr);
                type = str_to_out_type(typestr);
                mult = str_to_out_mult(multstr);
                domain = str_to_out_domain(domainstr);

                // Add OUTVAR to stream
                set_output_var(&((*streams)[streamnum]), varname, outvarnum,
                               format, type, mult, agg_type, domain);
                outvarnum++;
            }
        }      
    }
    // Check that the number of OUTVAR entries matches the number specified in OUTFILE.
    if (streamnum >= 0 && outvarnum != (int)(*streams)[streamnum].nvars) { 
        log_err("Output stream %s specifies %zu variables, but %d OUTVAR entries were found.", 
                (*streams)[streamnum].prefix, (*streams)[streamnum].nvars, outvarnum); 
    }
    
    free(ivar);
}
