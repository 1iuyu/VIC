/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine computes the current fraction of the vegetation band that is
 * covered with snow.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    Compute the current fraction of the vegetation band that is
 *           covered with snow.
 *****************************************************************************/
void
calc_snow_coverage(double             Cv,
                   bool               IS_GLAC,
                   double             snowfall,
                   snow_data_struct  *snow,
                   soil_con_struct   *soil_con)

{
    extern parameters_struct param;
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
    if (IS_GLAC == true) {
        double temp_intsnow = 0.0;
        double int_snow = param.TOL_A;
        double smr = 0.0;
        if (snow->swq == 0.) {
            if (snowfall > 0.0) {
                coverage = tanh(0.1 * snowfall);
            }
            else {
                coverage = 0.0;
            }
            snow->coverage = coverage;
        }
        else {
            if (fabs(snow->swq - snow->last_swq) > 0.0) {  // 存在积雪融化
                if (snowfall > 0.0) {
                    temp_intsnow = (snow->swq + snowfall) /
                        (0.5 * (cos(CONST_PI * (1.0 - pow(max(snow->coverage, param.TOL_A), 1.0 / 10.0))) + 1.0));
                    int_snow = min(1.e8, temp_intsnow);
                }
                smr = min(1.0, snow->swq / min(int_snow, 2000.0));

                snow->coverage = 1.0 - pow(acos(min(1.0, (2.0 * smr - 1.0))) / CONST_PI, 10.0); // n_melt = 200.0
            }
            if (snowfall > 0.0) {
                snow->coverage = snow->coverage + tanh(0.1 * snowfall) * (1.0 - snow->coverage);
            }
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
        snow->coverage = coverage;
    }
}
