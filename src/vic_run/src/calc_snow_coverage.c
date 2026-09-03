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
calc_snow_coverage(double            step_dt,
                   double            snowfall,
                   cell_data_struct *cell,
                   snow_data_struct *snow,
                   soil_con_struct  *soil_con,
                   int               flag)
{
    extern parameters_struct param;

    /* initialization */
    double coverage = 0.0;
    double topo_std = soil_con->topo_std;
    //double facter = 1.0e-4 * (topo_std * 3.3 - 9.5);
    double n_melt = 200.0 / max(10.0, topo_std);
    /* glacier snow cover fraction */
    if (cell->IS_GLAC) {
        n_melt = 10.0;
    }

    double new_snow = snowfall * step_dt;
    if (flag == SNOW_ACCU) {
        double k_accum = 1.15 * pow(topo_std, -0.55);
        if (snow->swq == 0.0) {
            if (snowfall > 0.0) {
                coverage = tanh(k_accum * new_snow);
            }
            else {
                coverage = 0.0;
            }
            snow->coverage = coverage;
        }
        else {
            if (snowfall > 0.0) {
                snow->coverage = snow->coverage + tanh(k_accum * new_snow) *
                                (1.0 - snow->coverage);
            }
            snow->coverage =
                min(1.0, max(0.0, snow->coverage));            
        }
        if (snow->coverage > param.TOL_A) {
            // Characteristic SWE corresponding to the snow-covered fraction before snowmelt.
            double temp_intsnow = (snow->swq + new_snow) / 
                        (0.5 * (cos(CONST_PI * pow(1.0 - snow->coverage, 1.0 / n_melt)) + 1.0));

            snow->ref_swq = min(1.e8, temp_intsnow);
        }
    }
    else if (flag == SNOW_MELT) {
        if (snow->swq <= 0.0) {
            snow->coverage = 0.0;
            snow->ref_swq = 0.0;
        }
        else if (snow->swq < snow->last_swq) {  // 存在积雪融化

            double smr = min(1.0, snow->swq / min(snow->ref_swq, 2000.0));
            snow->coverage = 1.0 - pow(acos(min(1.0, (2.0 * smr - 1.0))) / CONST_PI, n_melt); // n_melt = 200.0
            snow->coverage = min(1.0, max(0.0, snow->coverage));
        }
    }
}
