/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate the aerodynamic resistances.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    Calculate the aerodynamic resistance for each vegetation layer,
 *           based on Monin-Obukhov (M-O) Similarity Theory (MOST).
 *****************************************************************************/
double
initialize_MOST(double  wind,
                double  dthv,    
                double  zLdisp,  // reference height "minus" zero displacement height [m]
                double  theta_v,
                double  roughness,
                double *corr_wind)
{
    double Ri_bulk = 0.;    // bulk Richardson number
    double zeta = 0.;
    double L_disp = 0.;     // Obukhov length scale (m)
    double wstar = 0.5;     // convective velocity scale [m/s]   
    if (dthv >= 0.0) {
        *corr_wind = max(wind, 0.1);
    }
    else {
        *corr_wind = sqrt(wind * wind + wstar * wstar);
    }
    Ri_bulk = CONST_G * dthv * zLdisp /
                      (theta_v * (*corr_wind) * (*corr_wind));
    if (Ri_bulk >= 0.0) {   // neutral or stable
        zeta = Ri_bulk * log(zLdisp / roughness) /
                                    (1.0 - 5.0 * min(Ri_bulk, 0.19));
        zeta = max(zeta, 0.01);
    }
    else {  // unstable
        zeta = Ri_bulk * log(zLdisp / roughness);
        zeta = max(-100.0, min(zeta, -0.01));
    }
    L_disp = zLdisp / zeta;

    return(L_disp);
}
