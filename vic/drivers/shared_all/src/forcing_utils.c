/******************************************************************************
 * @section DESCRIPTION
 *
 * Library of utilities for processing meteorological forcings.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

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
    // vp = q/(q+CONST_EPS*(1-q))*p;

    // approximation used in VIC
    vp = q * p / CONST_EPS;

    return vp;
}

/******************************************************************************
 * @brief   convert surface pressure (Pa) to density (kg/m3) based on
 *          pressure (p), vapor pressure (vp), and temperature
 *
 * @param t temperature
 * @param p pressure
 *
 * @return rho surface pressure
 *****************************************************************************/
double
air_density(double t,
            double p)
{
    double rho;

    // full equation
    // rho = (p*1000)/(Rd * *t+CONST_TKFRZ) + (pv*1000)/(Rv * *t+CONST_TKFRZ);

    // approximation used in VIC
    rho = p / (CONST_RDAIR * (CONST_TKFRZ + t));

    return rho;
}

/******************************************************************************
 * @brief   return 1 if it will snow, otherwise return 0
 *****************************************************************************/
char
will_it_snow(double *snowf,
             size_t  n)
{
    size_t i;

    for (i = 0; i < n; i++) {
        if (snowf[i] > 0.) {
            return 1;
        }
    }

    return 0;
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
