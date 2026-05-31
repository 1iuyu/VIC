/******************************************************************************
 * @section DESCRIPTION
 *
 * Initialize Obukhov length scale and wind speed.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    Initialize Obukhov length scale and wind speed for the 
 *           Monin-Obukhov (M-O) Similarity Theory (MOST).
 *****************************************************************************/
double
initialize_MOST(double  wind,
                double  dthv,    
                double  zLdisp,  // reference height "minus" zero displacement height [m]
                double  theta_v, // 虚位温
                double  roughness,
                double *corr_wind)
{
    double Ri_bulk = 0.;    // bulk Richardson number
    double zeta = 0.;
    double L_disp = 0.;     // Obukhov length scale (m)
    double wstar = 0.5;     // convective velocity scale [m/s]   
    if (dthv >= 0.0) {  // stable conditions
        *corr_wind = max(wind, 0.1);
    }
    else {  // unstable conditions
        *corr_wind = sqrt(wind * wind + wstar * wstar);
    }
    Ri_bulk = CONST_G * dthv * zLdisp /
                      (theta_v * (*corr_wind) * (*corr_wind));
    if (Ri_bulk >= 0.0) {   // neutral or stable
        zeta = Ri_bulk * log(zLdisp / roughness) /
                                    (1.0 - 5.0 * min(Ri_bulk, 0.19));
        zeta = min(0.5, max(zeta, 0.01));
    }
    else {  // unstable
        zeta = Ri_bulk * log(zLdisp / roughness);
        zeta = max(-100.0, min(zeta, -0.01));
    }
    L_disp = zLdisp / zeta;

    return(L_disp);
}
