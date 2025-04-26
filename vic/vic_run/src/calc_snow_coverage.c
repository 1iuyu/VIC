/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine computes the current fraction of the vegetation band that is
 * covered with snow. Note: This new snow cover fraction based on Noah-MP, see
 * Niu and Yang (2007, JGR) scheme for more details.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    Compute the current fraction of the vegetation band that is
 *           covered with snow.
 *****************************************************************************/
double
calc_snow_coverage(double  depth,
                   double  swq)
{
    double SnowDensBulk; // bulk density of snow [Kg/m3] 
    double MeltFac; // melting factor for snow cover frac
    double coverage; // snow cover fraction
    double SNOW_MeltFac; // snowmelt m parameter
    double SNOW_CoverFac; // snow cover factor [m]

    SNOW_MeltFac = 2.5; // Original Noah-MP hard-coded
    SNOW_CoverFac = 0.005; // Original Noah-MP hard-coded
    coverage = 0.;

    if (depth > 0.0) {
        SnowDensBulk = swq * MM_PER_M / depth;
        MeltFac = pow(SnowDensBulk / 100.0, SNOW_MeltFac);
        coverage = tanh(depth / (SNOW_CoverFac * MeltFac));
    }

    return(coverage);
}
