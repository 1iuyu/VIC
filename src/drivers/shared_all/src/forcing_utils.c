/******************************************************************************
 * @section DESCRIPTION
 *
 * Library of utilities for processing meteorological forcings.
 *****************************************************************************/

#include "vic_driver_shared_all.h"

/******************************************************************************
 * @brief    calculate 1d average
 *****************************************************************************/
double
average(double *ar,
        size_t  n)
{
    size_t i;
    double sum = 0.;

    if (n <= 0) {
        log_err("Divide by zero or negative");
    }
    else if (n == 1) {
        return ar[0];
    }
    else {
        for (i = 0; i < n; i++) {
            sum += ar[i];
        }
    }

    return sum / n;
}

/******************************************************************************
 * @brief   convert specific humidity (q) to vapor pressure (vp) based on
 *          pressure (p)
 *
 * @param q specific humidity
 * @param p pressure
 *
 * @return vp vapor pressure (units are the same as p)
 *****************************************************************************/
double
q_to_vp(double q,
        double p)
{
    double vp;

    // full equation
     vp = (q * p) / (0.622 + 0.378 * q);

    return vp;
}
/******************************************************************************
 * @brief   convert surface pressure (Pa) to density (kg/m3) based on
 *          pressure (p), specific humidity (kg/kg), and temperature
 *
 * @param t temperature
 * @param p pressure
 * @param q specific humidity
 *
 * @return rho surface pressure
 *****************************************************************************/
double
air_density(double t,
            double p,
            double q)
{
    // full equation
    double egcm = q * p / (CONST_DEPS + (1.0 - CONST_DEPS) * q);

    double rho = (p - (1.0 - CONST_DEPS) * egcm) / (CONST_RDAIR * t);

    return rho;
}

/******************************************************************************
 * @brief   compute atmospheric potential temperature (K) based on
 *          pressure (p), and temperature
 *
 * @param t temperature
 * @param p pressure
 *
 * @return theta atmospheric potential temperature (K)
 *****************************************************************************/
double
compute_theta(double t,
              double p)
{
    double theta;
    double P0;
    P0 = CONST_PSREF * 100;
    // full equation
    theta = t * pow((P0 / p), CONST_RDAIR / CONST_CPDAIR);

    return theta;
}

/******************************************************************************
 * @brief   compute atmospheric virtual potential temperature (K) based on
 *          specific humidity (kg/kg), and potential temperature (K)
 *
 * @param q specific humidity (kg/kg)
 * @param theta potential temperature (K)
 *
 * @return theta_v atmospheric virtual potential temperature (K)
 *****************************************************************************/
double
compute_theta_v(double q,
                double theta)
{
    double theta_v;

    // full equation
    theta_v = theta * (1.0 + 0.61 * q);

    return theta_v;
}
/******************************************************************************
 * @brief    This routine downscales incoming shortwave radiation from the
 *           grid-cell scale to the elevation-band scale using the difference
 *           in atmospheric pressure between the grid cell and elevation band.
 *****************************************************************************/
double
downscale_shortwave(double sw_cell,
                    double p_cell,
                    double p_band,
                    double coszen,
                    double julian_day)
{
    double rt_R;
    double toa_swrad;
    double clr_idx;
    double diff_wgt;
    double diff_swrad;
    double beam_swrad;
    double beam_swrad_band;
    double diff_swrad_band;
    double k_band;
    double opt_factor;
    double sw_band;

    /* No incoming solar radiation */
    if (coszen <= 0.0 || sw_cell <= 0.0) {
        return 0.0;
    }

    /* Earth-Sun distance correction */
    rt_R = 1.0 - 0.01672 * cos(0.9856 * (julian_day - 4));

    /* Top-of-atmosphere shortwave radiation */
    toa_swrad = CONST_SOLAR * rt_R * rt_R * coszen;

    /* Clearness index */
    if (toa_swrad < 1.0e-7) {
        clr_idx = 0.0;
    }
    else {
        clr_idx = sw_cell / toa_swrad;
    }

    clr_idx = min(max(clr_idx, 0.0), 1.0);

    /* Diffuse radiation weight */
    diff_wgt = 0.952 - 1.041 * exp(-exp(min(2.3 - 4.702 * clr_idx, 3.5)));

    diff_wgt = min(max(diff_wgt, 0.0), 1.0);

    /* Separate diffuse and beam radiation */
    diff_swrad = sw_cell * diff_wgt;
    beam_swrad = sw_cell * (1.0 - diff_wgt);

    /* Broadband attenuation coefficient */
    if (clr_idx <= 0.0) {
        k_band = 0.0;
    }
    else {
        k_band = log(clr_idx) / p_band;
    }

    // Optical path correction caused by pressure difference.
    opt_factor = exp(k_band * (p_cell - p_band));

    if (!isfinite(opt_factor) ||
        opt_factor > 10000.0 ||
        opt_factor < -10000.0) {
        opt_factor = 0.0;
    }

    /* Flat terrain: illumination factor = 1 */
    beam_swrad_band = beam_swrad * opt_factor;

    /* Diffuse radiation is not downscaled */
    diff_swrad_band = diff_swrad;

    sw_band = beam_swrad_band + diff_swrad_band;

    /* Limit the downscaled radiation to ±50% */
    sw_band = min(sw_band, sw_cell * 1.5);
    sw_band = max(sw_band, sw_cell * 0.5);

    if (sw_band < 1.0e-4) {
        sw_band = 0.0;
    }

    return sw_band;
}
/******************************************************************************
 * @brief    This routine downscales incoming longwave radiation from the
 *           grid-cell scale to the elevation-band scale using atmospheric
 *           temperature, specific humidity, and pressure.
 *****************************************************************************/
double
downscale_longwave(double t_cell,
                   double q_cell,
                   double p_cell,
                   double lw_cell,
                   double t_band,
                   double q_band,
                   double p_band)
{
    double es_cell, es_band;
    double vp_cell, vp_band;
    double emissivity_clear_cell;
    double emissivity_clear_band;
    double emissivity_allsky_cell;
    double lw_band;

    /* Saturation vapor pressure: Pa */
    svp_flags(t_cell, p_cell,
              &es_cell, NULL, NULL, NULL, ESAT);

    svp_flags(t_band, p_band,
              &es_band, NULL, NULL, NULL, ESAT);

    /* Water vapor pressure: Pa */
    vp_cell = q_cell * es_cell / 100.0;
    vp_band = q_band * es_band / 100.0;

    /* Clear-sky emissivity */
    emissivity_clear_cell = 0.23 + 0.43 * pow(vp_cell / t_cell, 1.0 / 5.7);

    emissivity_clear_band = 0.23 + 0.43 * pow(vp_band / t_band, 1.0 / 5.7);

    /* All-sky emissivity at cell scale */
    emissivity_allsky_cell = lw_cell / (5.67e-8 * pow(t_cell, 4));

    /* Downscale longwave radiation */
    lw_band = (emissivity_clear_band + (emissivity_allsky_cell - emissivity_clear_cell)) *
                5.67e-8 * pow(t_band, 4);

    /* Limit the downscaling range to ±50% */
    lw_band = min(lw_band, lw_cell * 1.5);
    lw_band = max(lw_band, lw_cell * 0.5);

    return max(lw_band, 0.0);
}

/******************************************************************************
 * @brief    This routine determines the counts the number of forcing variables
             in each forcing file specified in the global parameter file.
 *****************************************************************************/
size_t
count_force_vars(FILE *gp)
{
    size_t        nvars;
    unsigned long start_position;
    char          cmdstr[MAXSTRING];
    char          optstr[MAXSTRING];

    // Figure out where we are in the input file
    fflush(gp);
    start_position = ftell(gp);

    // read the first line
    fgets(cmdstr, MAXSTRING, gp);

    // initalize nvars
    nvars = 0;

    // Loop through the lines
    while (!feof(gp)) {
        if (cmdstr[0] != '#' && cmdstr[0] != '\n' && cmdstr[0] != '\0') {
            // line is not blank or a comment
            sscanf(cmdstr, "%s", optstr);

            // if the line starts with FORCE_TYPE
            if (strcasecmp("FORCE_TYPE", optstr) == 0) {
                nvars++;
            }
            // else if we arive at another forcing file break out of loop
            else if (strcasecmp("FORCING1", optstr) == 0 ||
                     strcasecmp("FORCING2", optstr) == 0) {
                break;
            }
        }
        fgets(cmdstr, MAXSTRING, gp);
    }

    // put the position in the file back to where we started
    fseek(gp, start_position, SEEK_SET);

    return nvars;
}
