/******************************************************************************
* @section DESCRIPTION
*
* This routine was written to Calculate fluxes absorbed by vegetation,
* reflected by vegetation, and transmitted through vegetation for unit 
* incoming direct or diffuse flux given an underlying ground with known albedo.
******************************************************************************/

#include "vic_run.h"

/******************************************************************************
* @brief   This routine was written to Calculate fluxes absorbed by vegetation,
*          reflected by vegetation, and transmitted through vegetation for unit
*          incoming direct or diffuse flux given an underlying ground with known albedo.
******************************************************************************/
void
canopy_two_stream(double             coszen,
                  energy_bal_struct *energy,
                  cell_data_struct  *cell,
                  veg_var_struct    *veg_var,
                  veg_lib_struct    *veg_lib)
{
	extern parameters_struct param;
    extern option_struct     options;
	double gap_fracdir;
	double gap_fracdfs;
	double veg_density;
	double can_depth;
	double tmp_angle;
	double can_radius_z;
	double gap_bc;
	double fvd;
	double tmp_veg;
	double gap_wc;
	double tmp_coszen;
	double LOI;
	double Phi1, Phi2;
	double tau_leafdir;
	double tau_leafdfs;
	double scat_leaf;
	double tmp_0, tmp_1, tmp_2;
    double tmp_3, tmp_4, tmp_5;
	double tmp_6, tmp_7, tmp_8;
	double tmp_9;
	double omega_0;
    double proj_area;
	double coef_leafdir;
	double coef_leafdfs;
	double scat_canopy;
	double coef_candir;
	double coef_candfs;
	double B, C, D, F, H;
	double P1, P2, P3, P4;
	double S1, S2;
	double sigma;
	double U1, U2, U3;
	double D1, D2;
	double H1, H2, H3, H4, H5;
	double H6, H7, H8, H9, H10;
    double LAIcanopy;
    size_t i, j, Ncanopy;
    double A1, A2;
    double fabd_sun, fabd_sha;
    double fabi_sun, fabi_sha;
    double V, U, DV, DU;
    double dH2, dH3, dH5, dH6;
    double dH7, dH8, dH9, dH10;
    double dA1, dA2;
    double d_ftid, d_fabd;
    double d_ftii, d_fabi;
    double d_fabd_sun, d_fabd_sha;
    double d_fabi_sun, d_fabi_sha;
    double *AbsDirSun = energy->AbsDirSun;
    double *AbsDirSha = energy->AbsDirSha;
    double *AbsDfsSun = energy->AbsDfsSun;
    double *AbsDfsSha = energy->AbsDfsSha;
    // set pointer
    double *ReflectVeg = energy->ReflectVeg;
    double *TransmitVeg = energy->TransmitVeg;
    double *ReflGrndDir = energy->ReflGrndDir;
    double *ReflGrndDfs = energy->ReflGrndDfs;
    double *AbsSubDir = energy->AbsSubDir;
    double *AbsSubDfs = energy->AbsSubDfs;
    double *ReflSubDir = energy->ReflSubDir;
    double *ReflSubDfs = energy->ReflSubDfs;
    double *ShortDir2Dir = energy->ShortDir2Dir;
    double *ShortDfs2Dir = energy->ShortDfs2Dir;
    double *ShortDir2Dfs = energy->ShortDir2Dfs;
    double *ShortDfs2Dfs = energy->ShortDfs2Dfs;
    double *AlbedoGrndDir = energy->AlbedoGrndDir;
    double *AlbedoGrndDfs = energy->AlbedoGrndDfs;
    double *AlbedoSurfDir = energy->AlbedoSurfDir;
    double *AlbedoSurfDfs = energy->AlbedoSurfDfs;
    double *LAI_z = veg_var->LAI_z;
    double *SAI_z = veg_var->SAI_z;
    double *fsun_z = veg_var->fsun_z;
    double *LAIsun_z = veg_var->LAIsun_z;
    double *LAIsha_z = veg_var->LAIsha_z;
    double COI = veg_lib->COI;
    double NetLAI = veg_var->NetLAI;
    double NetSAI = veg_var->NetSAI;
    double fcanopy = veg_var->fcanopy;
    double wetFrac = veg_var->wetFrac;
    double Tfoliage = energy->Tfoliage;
    double Canopy_Upper = veg_lib->Canopy_Upper;
    double Canopy_Lower = veg_lib->Canopy_Lower;
    double Canopy_Radius = veg_lib->Canopy_Radius;
    double NetVEG = NetLAI + NetSAI;
    Ncanopy = cell->Ncanopy;
	/* compute within and between gaps */
    if (NetVEG == 0.0) {
        gap_fracdir = 1.0;
        gap_fracdfs = 1.0;
    }
    else {
        veg_density = -log(max(1.0 - fcanopy, 0.01)) / 
                            (CONST_PI * pow(Canopy_Radius, 2.0));
        can_depth = Canopy_Upper - Canopy_Lower;
        can_radius_z = 0.5 * can_depth;
        tmp_angle = atan(can_radius_z / Canopy_Radius * 
                                        tan(acos(max(0.01, coszen))));
        gap_bc = exp(-veg_density * CONST_PI * Canopy_Radius * 
                                        Canopy_Radius / cos(tmp_angle));
        fvd = NetVEG / (1.33 * CONST_PI * pow(Canopy_Radius, 3.0) *
                                (can_radius_z / Canopy_Radius) * veg_density);
        tmp_veg = can_depth * fvd;
        gap_wc = (1.0 - gap_bc) * exp(-0.5 * tmp_veg / coszen);
        gap_fracdir = min(1.0 - fcanopy, gap_bc + gap_wc);
        gap_fracdfs = 0.05;
    }
    tmp_coszen = max(0.001, coszen);
    LOI = min(max(COI, -0.4), 0.6);
    if (fabs(LOI) < 0.01) {
        LOI = 0.01;
    }
    Phi1 = 0.5 - 0.633 * LOI - 0.330 * LOI * LOI;
    Phi2 = 0.877 * (1.0 - 2.0 * Phi1);
    proj_area = Phi1 + Phi2 * tmp_coszen;
    tau_leafdir = proj_area / tmp_coszen;
    tau_leafdfs = (1.0 - Phi1 / Phi2 * log((Phi1 + Phi2) / Phi1)) / Phi2;
    tmp_0 = max(proj_area + Phi2 * tmp_coszen, param.TOL_A);
    tmp_1 = Phi1 * tmp_coszen;
    tmp_2 = 1.0 - tmp_1 / tmp_0 * log((tmp_1 + tmp_0) / tmp_1);
    for (i = 0; i < options.Nswband; i++) {
        scat_leaf = ReflectVeg[i] + TransmitVeg[i];
        omega_0 = 0.5 * scat_leaf * proj_area / tmp_0 * tmp_2;
        coef_leafdir = (1.0 + tau_leafdfs * tau_leafdir) / 
                    (scat_leaf * tau_leafdfs * tau_leafdir) * omega_0;
        coef_leafdfs = 0.5 * (ReflectVeg[i] + TransmitVeg[i] + (ReflectVeg[i] - 
                        TransmitVeg[i]) * pow((1.0 + LOI) / 2.0, 2.0)) / scat_leaf;
        /* adjust omega, betad, and betai for intercepted snow */
        if (Tfoliage > CONST_TKFRZ) {
            tmp_0 = scat_leaf;
            tmp_1 = coef_leafdir;
            tmp_2 = coef_leafdfs;
        }
        else {
            tmp_0 = (1.0 - wetFrac) * scat_leaf + wetFrac * param.SNOW_OMEGAS[i];
            tmp_1 = ((1.0 - wetFrac) * scat_leaf * coef_leafdir + wetFrac * 
                            param.SNOW_OMEGAS[i] * param.SNOW_BETADS) / tmp_0;
            tmp_2 = ((1.0 - wetFrac) * scat_leaf * coef_leafdfs + wetFrac * 
                            param.SNOW_OMEGAS[i] * param.SNOW_BETAIS) / tmp_0;
        }
        scat_canopy = tmp_0;
        coef_candir = tmp_1;
        coef_candfs = tmp_2;
        /* absorbed, reflected, transmitted fluxes per unit incoming radiation */
        B = 1.0 - scat_canopy + scat_canopy * coef_candfs;
        C = scat_canopy * coef_candfs;
        tmp_0 = tau_leafdfs * tau_leafdir;
        D = tmp_0 * scat_canopy * coef_candir;
        F = tmp_0 * scat_canopy * (1.0 - coef_candir);
        tmp_1 = B * B - C * C;
        H = sqrt(tmp_1) / tau_leafdfs;
        sigma = tmp_0 * tmp_0 - tmp_1;
        if (fabs(sigma) < param.TOL_A) {
            sigma = sign(param.TOL_A, sigma);
        }
        P1 = B + tau_leafdfs * H;
        P2 = B - tau_leafdfs * H;
        P3 = B + tmp_0;
        P4 = B - tmp_0;
        S1 = exp(-1.0 * min(H * NetVEG, 40.0));
        S2 = exp(-1.0 * min(tau_leafdir * NetVEG, 40.0));
        // direct
        U1 = B - C / AlbedoGrndDfs[i];
        U2 = B - C * AlbedoGrndDfs[i];
        U3 = F + C * AlbedoGrndDfs[i];                        
        tmp_2 = U1 - tau_leafdfs * H;
        tmp_3 = U1 + tau_leafdfs * H;
        D1 = P1 * tmp_2 / S1 - P2 * tmp_3 * S1;
        if (fabs(D1) < param.TOL_A) {
            D1 = sign(param.TOL_A, D1);
        }
        tmp_4 = U2 + tau_leafdfs * H;
        tmp_5 = U2 - tau_leafdfs * H;
        D2 = tmp_4 / S1 - tmp_5 * S1;
        if (fabs(D2) < param.TOL_A) {
            D2 = sign(param.TOL_A, D2);
        }
        H1 = -D * P4 - C * F;
        tmp_6 = D - H1 * P3 / sigma;
        tmp_7 = (D - C - H1 / sigma * (U1 + tmp_0)) * S2;
        H2 = (tmp_6 * tmp_2 / S1 - P2 * tmp_7) / D1;
        H3 = -(tmp_6 * tmp_3 * S1 - P1 * tmp_7) / D1;
        H4 = -F * P3 - C * D;
        tmp_8 = H4 / sigma;
        tmp_9 = (U3 - tmp_8 * (U2 - tmp_0)) * S2;
        H5 = -(tmp_8 * tmp_4 / S1 + tmp_9) / D2;
        H6 = (tmp_8 * tmp_5 * S1 + tmp_9) / D2;
        AlbedoSurfDir[i] = (H1 / sigma + H2 + H3) * (1.0 - gap_fracdir) + 
                            AlbedoGrndDir[i] * gap_fracdir;
        ReflSubDir[i] = (H1 / sigma + H2 + H3) * (1.0 - gap_fracdir);
        ReflGrndDir[i] = AlbedoGrndDir[i] * gap_fracdir;
        ShortDfs2Dir[i] = (H4 * S2 / sigma + H5 * S1 + H6 / S1) * (1.0 - gap_fracdir);
        ShortDir2Dir[i] = S2 * (1.0 - gap_fracdir) + gap_fracdir;
        AbsSubDir[i] = 1.0 - AlbedoSurfDir[i] - 
                 (1.0 - AlbedoGrndDir[i]) * ShortDir2Dir[i] -
                 (1.0 - AlbedoGrndDfs[i]) * ShortDfs2Dir[i];
        A1 = H1 / sigma * (1.0 - S2 * S2) /
            (2.0 * tau_leafdir) + H2 * (1.0 - S2 * S1) / (tau_leafdir + H) +
            H3 * (1.0 - S2 / S1) / (tau_leafdir - H);

        A2 = H4 / sigma * (1.0 - S2 * S2) / (2.0 * tau_leafdir) + H5 *
            (1.0 - S2 * S1) / (tau_leafdir + H) + H6 * (1.0 - S2 / S1) /
            (tau_leafdir - H);
        fabd_sun = (1.0 - scat_canopy) * (1.0 - S2 + (A1 + A2) / tau_leafdfs);
        fabd_sha = AbsSubDir[i] - fabd_sun;
        // Diffuse
        U1 = B - C / AlbedoGrndDfs[i];
        U2 = B - C * AlbedoGrndDfs[i];
        tmp_2 = U1 - tau_leafdfs * H;
        tmp_3 = U1 + tau_leafdfs * H;
        D1 = P1 * tmp_2 / S1 - P2 * tmp_3 * S1;
        tmp_4 = U2 + tau_leafdfs * H;
        tmp_5 = U2 - tau_leafdfs * H;
        D2 = tmp_4 / S1 - tmp_5 * S1;
        H7 = C * tmp_2 / (D1 * S1);
        H8 = -C * tmp_3 * S1 / D1;
        H9 = tmp_4 / (D2 * S1);
        H10 = -tmp_5 * S1 / D2;
        AlbedoSurfDfs[i] = (H7 + H8) * (1.0 - gap_fracdfs) + AlbedoGrndDfs[i] * gap_fracdfs;
        ReflSubDfs[i] = AlbedoSurfDfs[i];
        ReflGrndDfs[i] = 0.0;
        ShortDfs2Dfs[i] = (H9 * S1 + H10 / S1) * (1.0 - gap_fracdfs) + gap_fracdfs;
        ShortDir2Dfs[i] = 0.0;
        AbsSubDfs[i] = 1.0 - AlbedoSurfDfs[i] - (1.0 - AlbedoGrndDir[i]) *
                    ShortDir2Dfs[i] - (1.0 - AlbedoGrndDfs[i]) * ShortDfs2Dfs[i];
        A1 = H7 * (1.0 - S2 * S1) / (tau_leafdir + H) + H8 * (1.0 - S2 / S1) / (tau_leafdir - H);
        A2 = H9 * (1.0 - S2 * S1) / (tau_leafdir + H) + H10 * (1.0 - S2 / S1) / (tau_leafdir - H);
        fabi_sun = (1.0 - scat_canopy) * (A1 + A2) / tau_leafdfs;
        fabi_sha = AbsSubDfs[i] - fabi_sun;
        if (i == 0) {
            if (Ncanopy == 1) {
                fsun_z[0] = (1.0 - S2) / min(tau_leafdir * NetVEG, 40.0);
                LAIsun_z[0] = LAI_z[0] * fsun_z[0];
                LAIsha_z[0] = LAI_z[0] * (1.0 - fsun_z[0]);
                AbsDirSun[0] = fabd_sun / (fsun_z[0] * NetVEG);
                AbsDfsSun[0] = fabi_sun / (fsun_z[0] * NetVEG);
                AbsDirSha[0] = fabd_sha / ((1.0 - fsun_z[0]) * NetVEG);
                AbsDfsSha[0] = fabi_sha / ((1.0 - fsun_z[0]) * NetVEG);
                double ksun_vcmax = (1.0 - exp(-(tau_leafdir + 0.30) * NetLAI)) / (tau_leafdir + 0.30);
                double ksha_vcmax = (1.0 - exp(-0.30 * NetLAI)) / 0.30 - ksun_vcmax;
                if (NetLAI > 0.0) {
                    ksun_vcmax /= (fsun_z[0] * NetLAI);
                    ksha_vcmax /= ((1.0 - fsun_z[0]) * NetLAI);
                }
                else {
                    ksun_vcmax = 0.0;
                    ksha_vcmax = 0.0;
                }
                veg_var->ksun_vcmax = ksun_vcmax;
                veg_var->ksha_vcmax = ksha_vcmax;
            }
            else {
                for (j = 0; j < Ncanopy; j++) {
                    if (j == 0) {
                        LAIcanopy = 0.5 * (LAI_z[j] + SAI_z[j]);
                    }
                    else {
                        LAIcanopy += 0.5 * (LAI_z[j-1] + SAI_z[j-1] + LAI_z[j] + SAI_z[j]);
                    }
                    S1 = exp(-min(H * LAIcanopy, 40.0));
                    S2 = exp(-min(tau_leafdir * LAIcanopy, 40.0));
                    fsun_z[j] = S2;
                    LAIsun_z[j] = LAI_z[j] * fsun_z[j];
                    LAIsha_z[j] = LAI_z[j] * (1.0 - fsun_z[j]);
                    V = D1;
                    DV = H * P1 * tmp_2 / S1 + H * P2 * tmp_3 * S1;
                    U = tmp_6 * tmp_2 / S1 - P2 * tmp_7;
                    DU = H * tmp_6 * tmp_2 / S1 + tau_leafdir * P2 * tmp_7;
                    dH2 = (V * DU - U * DV) / (V * V);

                    U = -tmp_6 * tmp_3 * S1 + P1 * tmp_7;
                    DU = H * tmp_6 * tmp_3 * S1 - tau_leafdir * P1 * tmp_7;
                    dH3 = (V * DU - U * DV) / (V * V);

                    V = D2;
                    DV = H * tmp_4 / S1 + H * tmp_5 * S1;

                    U = -H4 / sigma * tmp_4 / S1 - tmp_9;
                    DU = -H * H4 / sigma * tmp_4 / S1 + tau_leafdir * tmp_9;
                    dH5 = (V * DU - U * DV) / (V * V);

                    U = H4 / sigma * tmp_5 * S1 + tmp_9;
                    DU = -H * H4 / sigma * tmp_5 * S1 - tau_leafdir * tmp_9;
                    dH6 = (V * DU - U * DV) / (V * V);

                    dA1 = H1 / sigma * S2 * S2
                        + H2 * S2 * S1
                        + H3 * S2 / S1
                        + (1.0 - S2 * S1) / (tau_leafdir + H) * dH2
                        + (1.0 - S2 / S1) / (tau_leafdir - H) * dH3;

                    dA2 = H4 / sigma * S2 * S2
                        + H5 * S2 * S1
                        + H6 * S2 / S1
                        + (1.0 - S2 * S1) / (tau_leafdir + H) * dH5
                        + (1.0 - S2 / S1) / (tau_leafdir - H) * dH6;
                    
                    d_ftid = -tau_leafdir * H4 / sigma * S2
                            - H * H5 * S1
                            + H * H6 / S1
                            + dH5 * S1
                            + dH6 / S1;

                    d_fabd = -(dH2 + dH3)
                            + (1.0 - AlbedoGrndDir[i]) * tau_leafdir * S2
                            - (1.0 - AlbedoGrndDfs[i]) * d_ftid;

                    d_fabd_sun = (1.0 - scat_canopy) *
                                (tau_leafdir * S2 +
                                (dA1 + dA2) / tau_leafdfs);
                    d_fabd_sha = d_fabd - d_fabd_sun;
                    AbsDirSun[j] = max(d_fabd_sun, 0.0) / max(fsun_z[j], param.TOL_A);
                    AbsDirSha[j] = max(d_fabd_sha, 0.0)
                                    / max(1.0 - fsun_z[j], param.TOL_A);
                    // Diffuse
                    V = D1;
                    DV = H * P1 * tmp_2 / S1 + H * P2 * tmp_3 * S1;
                    U = C * tmp_2 / S1;
                    DU = H * C * tmp_2 / S1;
                    dH7 = (V * DU - U * DV) / (V * V);
                    U = -C * tmp_3 * S1;
                    DU = H * C * tmp_3 * S1;
                    dH8 = (V * DU - U * DV) / (V * V);
                    V = D2;
                    DV = H * tmp_4 / S1 + H * tmp_5 * S1;
                    U = tmp_4 / S1;
                    DU = H * tmp_4 / S1;
                    dH9 = (V * DU - U * DV) / (V * V);
                    U = -tmp_5 * S1;
                    DU = H * tmp_5 * S1;
                    dH10 = (V * DU - U * DV) / (V * V);
                    dA1 = H7 * S2 * S1
                        + H8 * S2 / S1
                        + (1.0 - S2 * S1) / (tau_leafdir + H) * dH7
                        + (1.0 - S2 / S1) / (tau_leafdir - H) * dH8;
                    dA2 = H9 * S2 * S1
                        + H10 * S2 / S1
                        + (1.0 - S2 * S1) / (tau_leafdir + H) * dH9
                        + (1.0 - S2 / S1) / (tau_leafdir - H) * dH10;
                    d_ftii = -H * H9 * S1 + H * H10 / S1 + dH9 * S1 + dH10 / S1;
                    d_fabi = -(dH7 + dH8) - (1.0 - AlbedoGrndDfs[i]) * d_ftii;
                    d_fabi_sun = (1.0 - scat_canopy) * (dA1 + dA2) / tau_leafdfs;
                    d_fabi_sha = d_fabi - d_fabi_sun;
                    AbsDfsSun[j] = max(d_fabi_sun, 0.0) / max(fsun_z[j], param.TOL_A);
                    AbsDfsSha[j] = max(d_fabi_sha, 0.0) / max(1.0 - fsun_z[j], param.TOL_A);
                }
            }
        }
    }
}

