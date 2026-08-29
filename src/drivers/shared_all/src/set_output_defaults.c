/******************************************************************************
 * @section DESCRIPTION
 *
 * Set the output_stream and out_data structures to default values. These can
 * be overridden by the user in the global control file.
 *****************************************************************************/

#include "vic_driver_shared_all.h"

/******************************************************************************
 * @brief    Set the output_stream and out_data structures to default values.
             These can be overridden by the user in the global control file.
 *****************************************************************************/
void
set_output_defaults(stream_struct **streams,
                    dmy_struct     *dmy_current,
                    unsigned short  default_file_format)
{
    extern option_struct options;

    size_t               streamnum;
    size_t               varnum;
    alarm_struct         default_alarm;
    int                  default_freq_n = 1;


    set_alarm(dmy_current, FREQ_NDAYS, &default_freq_n, &default_alarm);

    for (streamnum = 0; streamnum < options.Noutstreams; streamnum++) {
        (*streams)[streamnum].agg_alarm = default_alarm;
        (*streams)[streamnum].file_format = default_file_format;
    }

    // Variables in first file
    streamnum = 0;
    varnum = 0;
    strcpy((*streams)[streamnum].prefix, "fluxes");
    set_output_var(&((*streams)[streamnum]), "OUT_PREC", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT, OUT_DOMAIN_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_EVAP", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT, OUT_DOMAIN_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_RUNOFF", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT, OUT_DOMAIN_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_BASEFLOW", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT, OUT_DOMAIN_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_SOIL_TEMP", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT, OUT_DOMAIN_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_SOIL_LIQ", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT, OUT_DOMAIN_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_SWNET", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT, OUT_DOMAIN_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_LATENT", varnum++, "%.4f",
                    OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT, OUT_DOMAIN_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_LWNET", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT, OUT_DOMAIN_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_SURF_TEMP", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT, OUT_DOMAIN_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_EVAP_BARE", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT, OUT_DOMAIN_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_SENSIBLE", varnum++,
                    "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT, OUT_DOMAIN_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_GRND_FLUX", varnum++,
                    "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT, OUT_DOMAIN_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_VEGT", varnum++,
                   "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT, OUT_DOMAIN_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_CANOPY_SWQ", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT, OUT_DOMAIN_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_AIR_TEMP", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT, OUT_DOMAIN_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_WIND", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT, OUT_DOMAIN_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_SWE", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT, OUT_DOMAIN_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_SNOW_DEPTH", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT, OUT_DOMAIN_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_SNOW_PACK_ICE", varnum++,
                   "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT, OUT_DOMAIN_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_SNOW_PACK_LIQ", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT, OUT_DOMAIN_DEFAULT);
}
