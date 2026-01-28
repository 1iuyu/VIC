/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine computes the current fraction of the vegetation band that is
 * covered with snow.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    Compute the current fraction of the vegetation band that is
 *           covered with snow.
 *****************************************************************************/
void
calc_snow_coverage(double             Cv,
                   unsigned short     veg_class,
                   snow_data_struct  *snow,
                   soil_con_struct   *soil_con)

{
    extern option_struct     options;

    double MeltFac;
    double coverage;
    double GridSize;
    double density;
    double SNOW_MeltFac;
    double SNOW_CoverFac;
    double gridScalePara;  // GridSize: grid spacing (500-36000)[m].

    double snow_depth = snow->snow_depth;   // snow depth [m]

    /* initialization */
    coverage = 0.;
    /* glacier snow cover fraction */
    if (veg_class == options.GLACIER_ID) {
        if (snow->swq > 0.) {
            coverage = 1.0;
        }
    }
    /* ground snow cover fraction */
    else {
        if (snow_depth > 0.0) {
            GridSize = sqrt(soil_con->cell_area * Cv);
            gridScalePara = min((max(GridSize, 500.0) / M_PER_KM), 36.0);
            SNOW_MeltFac = 0.9713 + tanh(0.7436 * gridScalePara);
            SNOW_CoverFac = 0.0062 * sinh(0.0555 * gridScalePara) + 0.0555;
            density = snow->swq / snow_depth;
            MeltFac = pow(density / 100, SNOW_MeltFac);
            coverage = tanh(snow_depth / (SNOW_CoverFac * MeltFac));
        }
    }
    snow->coverage = coverage;
}
