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
int
CalcAerodynamic(size_t  iter,
                double  Qair,
                double  air_temp,    
                double *zL_disp,       /* M-O stability (z/L), above displacement */
                double *ustar,         /* friction velocity [m/s] */
                double *aero_resist,   /* aerodynamic resistances */
                double *tmp_sensible,
                double  wind,          /* adjusted wind speed */
                double  air_density,
                double  displacement,  /* vegetation displacement */
                double  ref_height,    /* vegetation reference height */
                double  roughness)
{
    extern parameters_struct param;

    size_t      moz_signchg_count;
    double      tmp_CM;
    double      tmp_CH;
    double      tmp_CM2;
    double      tmp_CH2; 
    double      L_disp;
    double      zL_disp_2m;
    double      tmp_Tvir;
    double      old_Z;
    double      TMP1;
    double      TMP2;
    double      TMP3;    
    double      psi_m;
    double      psi_h;
    double      psi_m_2m;
    double      psi_h_2m;
    double      phi_m;
    double      phi_h;
    double      phi_m_2m;
    double      phi_h_2m;

    /* initialization */
    old_Z = (*zL_disp);
    moz_signchg_count = 0;
    if (ref_height < displacement) {
        log_err("Warning: ref_height < displacement.");
    }

    /* temporary drag coefficients */
    tmp_CM = log((ref_height - displacement) / roughness);
    tmp_CH = log((ref_height - displacement) / roughness);
    tmp_CM2 = log((2.0 + roughness) / roughness);
    tmp_CH2 = log((2.0 + roughness) / roughness);

    /* compute M-O stability parameter */
    if (iter == 0) {
        (*ustar) = 0.;
        L_disp = 0.;
        (*zL_disp) = 0.;
        zL_disp_2m = 0.;
    }
    else {
        tmp_Tvir = (1.0 + 0.61 * Qair) * air_temp;
        TMP1 = CONST_KARMAN * (CONST_G / tmp_Tvir) * 
                    (*tmp_sensible) / (air_density * CONST_CPDAIR);
        if (fabs(TMP1) < param.TOL_A) {
            TMP1 = param.TOL_A;
        }
        L_disp = -1. * pow(*ustar, 3) / TMP1;
        (*zL_disp) = min((ref_height - displacement) / L_disp, 1.0);
        zL_disp_2m = min((2.0 + roughness) / L_disp, 1.0);
    }

    /* accumulate number of times moz changes sign. */
    if (old_Z * (*zL_disp) < 0.) {
        moz_signchg_count++;
    }
    if (moz_signchg_count >= 2) {
        (*zL_disp) = 0.;
        psi_m = 0.;
        psi_h = 0.;
        zL_disp_2m = 0.;
        psi_m_2m = 0.;
        psi_h_2m = 0.;
    }

    /* evaluate stability-dependent variables using moz from prior iteration */
    if (*zL_disp < 0.) {
        TMP1 = pow(1.0 - 16.0 * (*zL_disp), 0.25);
        TMP2 = log((1.0 + TMP1 * TMP1) / 2.0);
        TMP3 = log((1.0 + TMP1) / 2.0);
        phi_m = 2.0 * TMP3 + TMP2 - 2.0 * atan(TMP1) + 1.5707963;
        phi_h = 2 * TMP2;
        // 2-meter quantities
        double TMP1_2m = pow((1.0 - 16.0 * zL_disp_2m), 0.25);
        double TMP2_2m = log((1.0 + TMP1_2m * TMP1_2m) / 2.0);
        double TMP3_2m = log((1.0 + TMP1_2m) / 2.0);
        phi_m_2m = 2.0 * TMP3_2m + TMP2_2m - 2.0 * atan(TMP1_2m) + 1.5707963;
        phi_h_2m = 2 * TMP2_2m;
    }
    else {
        phi_m = -5 * (*zL_disp);
        phi_h = phi_m;
        phi_m_2m = -5 * zL_disp_2m;
        phi_h_2m = phi_m_2m;
    }

    /* except for first iteration, weight stability factors for previous */
    if (iter == 0) {
        psi_m = phi_m;
        psi_h = phi_h;
        psi_m_2m = phi_m_2m;
        psi_h_2m = phi_h_2m;
    }
    else {
        psi_m = 0.5 * (psi_m + phi_m);
        psi_h = 0.5 * (psi_h + phi_h);
        psi_m_2m = 0.5 * (psi_m_2m + phi_m_2m);
        psi_h_2m = 0.5 * (psi_h_2m + phi_h_2m); 
    }

    /* exchange coefficients */
    psi_m = min(psi_m, 0.9 * tmp_CM);
    psi_h = min(psi_h, 0.9 * tmp_CH);
    psi_m_2m = min(psi_m_2m, 0.9 * tmp_CM2);
    psi_h_2m = min(psi_h_2m, 0.9 * tmp_CH2);

    double CMFM   = tmp_CM  - psi_m;
    double CHFH   = tmp_CH  - psi_h;
    double CM2FM2 = tmp_CM2 - psi_m_2m;
    double CH2FH2 = tmp_CH2 - psi_h_2m;

    if (fabs(CMFM) < param.TOL_A) {
        CMFM = param.TOL_A;
    }  
    if (fabs(CHFH) < param.TOL_A) {
        CHFH = param.TOL_A;
    }
    if (fabs(CM2FM2) < param.TOL_A) {
        CM2FM2 = param.TOL_A;
    }
    if (fabs(CH2FH2) < param.TOL_A) {
        CH2FH2 = param.TOL_A;
    }

    double CM = CONST_KARMAN * CONST_KARMAN / (CMFM * CMFM);
    double CH = CONST_KARMAN * CONST_KARMAN / (CMFM * CHFH);

    /* friction velocity */
    (*ustar) = wind * sqrt(CM);

    double CH_2m = CONST_KARMAN * (*ustar) / CH2FH2;

    /* aerodynamic resistance [0]Ra, [1]SH ,[2]LH */
    aero_resist[0] = max(1.0, 1.0 / (CM * wind));
    aero_resist[1] = max(1.0, 1.0 / (CH * wind));
    aero_resist[2] = aero_resist[1];

    return (0);
}
