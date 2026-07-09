/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes the lake variables for each new grid cell.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    This routine Calculate fraction of land surfaces which are 
 *  submerged based on surface microtopography and surface water storage.
 *****************************************************************************/
int
calc_surf_water(double            step_dt,
                snow_data_struct *snow,
                cell_data_struct *cell,
                soil_con_struct  *soil_con)
{
    extern parameters_struct param;
    size_t iter = 0;
    size_t MAX_ITER = 10;
    double fd;
    double dfdd;
    double beta = -3.0;
    double surf2soil = 0.0;
    double frac_h2o = cell->frac_h2o;
    double h2osfc = cell->h2osfc;
    double coverage = snow->coverage;
    double slope0 = pow(0.4, 1.0 / beta);
    double micro_sigma = pow(soil_con->slope + slope0, beta);
    if (frac_h2o > param.TOL_B) {
        if (cell->IS_WET) {
            double d = 0.0;
            double old_d = 0.0;
            double sigma = MM_PER_M * micro_sigma;
            do {
                old_d = d;
                fd = 0.5 * d * (1.0 + erf(d / (sigma * sqrt(2.0)))) +
                    sigma / sqrt(2.0 * CONST_PI) * exp(-d * d / 
                            (2.0 * sigma * sigma)) - h2osfc;
                dfdd = 0.5 * (1.0 + erf(d / (sigma * sqrt(2.0))));
                d = d - fd / dfdd;
                
                if (fabs(d - old_d) > param.TOL_B) break;
                iter++;
            } 
            while (iter < MAX_ITER);
            frac_h2o = 0.5 * (1.0 + erf(d / (sigma * sqrt(2.0))));
        }
        else if (cell->IS_GLAC) {
            double gamma = 1.375;
            double rho_ratio = CONST_RHOICE / CONST_RHOFW;
            double k = 0.0365 * pow(soil_con->cell_area, gamma - 1.0) * rho_ratio * MM_PER_M;
            if (h2osfc > 0.0 && k > 0.0) {
                double x = pow(h2osfc / k, 1.0 / gamma);  // 初始猜测
                double fx, dfdx;
                do {
                    fx = k * pow(x, gamma) - h2osfc;
                    dfdx = k * gamma * pow(x, gamma - 1.0);
                    
                    if (fabs(dfdx) < 1e-12) break;
                    
                    double dx = fx / dfdx;
                    x = x - dx;
                    
                    if (fabs(dx) < param.TOL_B) break;
                    iter++;
                }
                while (iter < MAX_ITER);
                
                frac_h2o = x;
                
                // 物理约束
                frac_h2o = max(min(0.0, frac_h2o), 1.0);

            } else {
                frac_h2o = 0.0;
            }
        }
    }
    else {
        frac_h2o = 0.0;
        surf2soil = h2osfc / step_dt;
    }
    if (coverage > (1.0 - frac_h2o) && snow->swq > 0.0) {
        if (frac_h2o > 0.01) {
            frac_h2o = max(0.01, 1.0 - coverage);
            coverage = 1.0 - frac_h2o;
        }
        else {
            coverage = 1.0 - frac_h2o;
        }
    }
    h2osfc -= surf2soil;
    cell->liq[0] += surf2soil * step_dt;
    snow->coverage = coverage;
    cell->h2osfc = h2osfc;
    cell->frac_h2o = frac_h2o;
    // Update the number of nodes in the soil column based on the presence of surface water
    if (cell->h2osfc > 0.0) {
        cell->Nnode = soil_con->Nbedrock + snow->Nsnow + 1;
    }
    else {
        cell->Nnode = soil_con->Nbedrock + snow->Nsnow;
    }
    
    return (0);
}