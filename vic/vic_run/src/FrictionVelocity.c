/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate the aerodynamic resistances.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    Calculate the aerodynamic resistance for each vegetation layer,
 *           based on Monin-Obukhov (M-O) Similarity Theory (MOST).
 *****************************************************************************/
int
FrictionVelocity(double  L_disp,
                 double  corr_wind,
                 double *ref_height,
                 double *ustar,
                 double *temp_profile,
                 double *Qair_profile,
                 double *Z0m_grnd,
                 double  displacement)
{
    double zeta = 0.;
    double zetam = 1.574;
    double zetat = 0.465;
    double zLdisp = 0.;

    // Wind profile
    zLdisp = ref_height[0] - displacement;

    zeta = zLdisp / L_disp;
    if (zeta < -zetam) {
        *ustar = CONST_KARMAN * corr_wind / (log(-zetam * L_disp / Z0m_grnd[0]) -
                    StabilityFunc1(-zetam) + 
                    StabilityFunc1(Z0m_grnd[0] / L_disp) +
                    1.14 * (pow(-zeta, 0.333) - pow(zetam, 0.333)));
    }
    else if (zeta < 0.) {
        *ustar = CONST_KARMAN * corr_wind / (log(zLdisp / Z0m_grnd[0]) -
                    StabilityFunc1(zeta) +
                    StabilityFunc1(Z0m_grnd[0] / L_disp));
    }
    else if (zeta <= 1.) {
        *ustar = CONST_KARMAN * corr_wind / (log(zLdisp / Z0m_grnd[0]) +
                    5.0 * zeta - 5.0 * Z0m_grnd[0] / L_disp);
    }
    else {
        *ustar = CONST_KARMAN * corr_wind / (log(L_disp / Z0m_grnd[0]) +
                    5.0 - 5.0 * Z0m_grnd[0] / 
                    L_disp + (5.0 * log(zeta) + zeta - 1.0));
    }

    // Temperature profile
    zLdisp = ref_height[1] - displacement;
    zeta = zLdisp / L_disp;

    if (zeta < -zetam) {
        // Unstable
        *temp_profile = CONST_KARMAN / (log(-zetat * L_disp/ Z0m_grnd[1]) -
                            StabilityFunc2(-zetat) +
                            StabilityFunc2(Z0m_grnd[1] / L_disp) +
                            0.8 * (pow(zetat, -0.333) - pow(-zeta, -0.333)));
    }
    else if (zeta < 0.) {
        // Unstable
         *temp_profile = CONST_KARMAN / (log(zLdisp / Z0m_grnd[1]) -
                            StabilityFunc2(zeta) +
                            StabilityFunc2(Z0m_grnd[1] / L_disp));
    }
    else if (zeta <= 1.) {
        // Stable
        *temp_profile = CONST_KARMAN / (log(zLdisp / Z0m_grnd[1]) +
                            5.0 * zeta - 5.0 * Z0m_grnd[1] / L_disp);
    }
    else {
        // Very stable
        *temp_profile = CONST_KARMAN / (log(L_disp / Z0m_grnd[1]) +
                            5.0 - 5.0 * Z0m_grnd[1] / L_disp +
                            (5.0 * log(zeta) + zeta - 1.0));
    }

    // Humidity profile
    if (Z0m_grnd[1] == Z0m_grnd[2]) {
        *Qair_profile = *temp_profile;
    }
    else {
        zLdisp = ref_height[2] - displacement;
        zeta = zLdisp / L_disp;
        if (zeta < -zetam) {
            // Unstable
            *Qair_profile = CONST_KARMAN / (log(-zetat * L_disp/ Z0m_grnd[2]) -
                                StabilityFunc2(-zetat) +
                                StabilityFunc2(Z0m_grnd[2] / L_disp) +
                                0.8 * (pow(zetat, -0.333) - pow(-zeta, -0.333)));
        }
        else if (zeta < 0.) {
            // Unstable
            *Qair_profile = CONST_KARMAN / (log(zLdisp / Z0m_grnd[2]) -
                                StabilityFunc2(zeta) +
                                StabilityFunc2(Z0m_grnd[2] / L_disp));
        }
        else if (zeta <= 1.) {
            // Stable
            *Qair_profile = CONST_KARMAN / (log(zLdisp / Z0m_grnd[2]) +
                                5.0 * zeta - 5.0 * Z0m_grnd[2] / L_disp);
        }
        else {
            // Very stable
            *Qair_profile = CONST_KARMAN / (log(L_disp / Z0m_grnd[2]) +
                                5.0 - 5.0 * Z0m_grnd[2] / L_disp +
                                (5.0 * log(zeta) + zeta - 1.0));
        }
    }

    return 0;

}
/******************************************************************************
 * @brief    Calculate the aerodynamic resistance for each vegetation layer,
 *           based on Monin-Obukhov (M-O) Similarity Theory (MOST).
 *****************************************************************************/
double
StabilityFunc1(double zeta)
{
    double chik, chik2;
    
    chik2 = sqrt(1.0 - 16.0 * zeta);
    chik = sqrt(chik2);
    
    return 2.0 * log((1.0 + chik) * 0.5) +
           log((1.0 + chik2) * 0.5) -
           2.0 * atan(chik) + CONST_PI * 0.5;
}
/******************************************************************************
 * @brief    Calculate the aerodynamic resistance for each vegetation layer,
 *           based on Monin-Obukhov (M-O) Similarity Theory (MOST).
 *****************************************************************************/
double
StabilityFunc2(double zeta)
{
    double chik2;
    
    chik2 = sqrt(1.0 - 16.0 * zeta);
    
    return 2.0 * log((1.0 + chik2) * 0.5);
}