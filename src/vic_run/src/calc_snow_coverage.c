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
                   double            air_temp,
                   double            snowfall,
                   double            rainfall,
                   cell_data_struct *cell,
                   snow_data_struct *snow)
{
    extern parameters_struct param;

    /* initialization */
    double coverage = 0.0;
    double n_melt = 200.0;
    /* glacier snow cover fraction */
    if (cell->IS_GLAC) {
        n_melt = 10.0;
    }

    double temp_intsnow = 0.0;
    double int_snow = param.TOL_A;
    double smr = 0.0;
    double new_snow = snowfall * step_dt;
    if (snow->swq == 0.0) {
        if (snowfall > 0.0) {
            coverage = tanh(0.1 * new_snow);
        }
        else {
            coverage = 0.0;
        }
        snow->coverage = coverage;
    }
    else {
        if (fabs(snow->swq - snow->last_swq) < 0.0) {  // 存在积雪融化
            if (snowfall > 0.0) {
                temp_intsnow = (snow->swq + new_snow) /
                    (0.5 * (cos(CONST_PI * (1.0 - pow(max(snow->coverage, param.TOL_A), 1.0 / n_melt))) + 1.0));
                int_snow = min(1.e8, temp_intsnow);
            }
            smr = min(1.0, snow->swq / min(int_snow, 2000.0));

            snow->coverage = 1.0 - pow(acos(min(1.0, (2.0 * smr - 1.0))) / CONST_PI, n_melt); // n_melt = 200.0
        }
        if (snowfall > 0.0) {
            snow->coverage = snow->coverage + tanh(0.1 * new_snow) * (1.0 - snow->coverage);
        }
    }

    /* snowpack water processs */
    update_snow(step_dt, air_temp,
                snowfall, rainfall, snow);

}
