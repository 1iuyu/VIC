/******************************************************************************
* @section DESCRIPTION
*
* This routine was written to Calculate fluxes absorbed by vegetation,
* reflected by vegetation, and transmitted through vegetation for unit 
* incoming direct or diffuse flux given an underlying ground with known albedo.
******************************************************************************/

#include <vic_run.h>

/******************************************************************************
* @brief   This routine was written to Calculate fluxes absorbed by vegetation,
*          reflected by vegetation, and transmitted through vegetation for unit
*          incoming direct or diffuse flux given an underlying ground with known albedo.
******************************************************************************/
void
canopy_two_stream(size_t	         lidx,
				  size_t	         index,
				  double   	         Tcanopy,
				  double	         NetLAI,
				  double	         NetSAI,
				  double	         fcanopy,
                  double             coszen,
				  double	        *proj_area,
				  double	         wetFrac,
                  energy_bal_struct *energy,
                  veg_lib_struct    *veg_lib)
{
	extern parameters_struct param;

	double		gap_frac_direct;
	double		gap_frac_diffuse;
	double		veg_density;
	double		can_depth;
	double		tmp_angle;
	double		can_radius_z;
	double		gap_bc;
	double		fvd;
	double		tmp_veg;
	double		gap_wc;
	double		tmp_coszen;
	double		LOI;
	double		Phi1, Phi2;
	double		tau_leaf_direct;
	double		tau_leaf_diffuse;
	double		scat_leaf;
	double		tmp_0, tmp_1, tmp_2;
    double      tmp_3, tmp_4, tmp_5;
	double		tmp_6, tmp_7, tmp_8;
	double		tmp_9;
	double		omega_0;
	double		coef_leaf_direct;
	double		coef_leaf_diffuse;
	double		scat_canopy;
	double		coef_can_direct;
	double		coef_can_diffuse;
	double		B, C, D, F, H;
	double		P1, P2, P3, P4;
	double		S1, S2;
	double		sigma;
	double		U1, U2, U3;
	double		D1, D2;
	double		H1, H2, H3, H4, H5;
	double		H6, H7, H8, H9, H10;
	double		transmit_sw_direct;
	double		transmit_sw_diffuse;
	double		refl_total;
	double		refl_over;
	double		refl_grnd;

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

    double  COI = veg_lib->COI;
    double  Canopy_Upper = veg_lib->Canopy_Upper;
    double  Canopy_Lower = veg_lib->Canopy_Lower;
    double  Canopy_Radius = veg_lib->Canopy_Radius;

	/* compute within and between gaps */
    if (NetLAI + NetSAI == 0.0) {
        gap_frac_direct = 1.0;
        gap_frac_diffuse = 1.0;
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
        fvd = (NetLAI + NetSAI) / (1.33 * CONST_PI * pow(Canopy_Radius, 3.0) *
                                (can_radius_z / Canopy_Radius) * veg_density);
        tmp_veg = can_depth * fvd;
        gap_wc = (1.0 - gap_bc) * exp(-0.5 * tmp_veg / coszen);
        gap_frac_direct = min(1.0 - fcanopy, gap_bc + gap_wc);
        gap_frac_diffuse = 0.05;
    }
    tmp_coszen = max(0.001, coszen);
    LOI = min(max(COI, -0.4), 0.6);
    if (LOI < 0.01) {
        LOI = 0.01;
    }
    Phi1 = 0.5 - 0.633 * LOI - 0.330 * LOI * LOI;
    Phi2 = 0.877 * (1.0 - 2.0 * Phi1);
    (*proj_area) = Phi1 + Phi2 * tmp_coszen;
    tau_leaf_direct = (*proj_area) / tmp_coszen;
    tau_leaf_diffuse = (1.0 - Phi1 / Phi2 * log((Phi1 + Phi2) / Phi1)) / Phi2;
    scat_leaf = ReflectVeg[lidx] + TransmitVeg[lidx];
    tmp_0 = (*proj_area) + Phi2 * tmp_coszen;
    tmp_1 = Phi1 * tmp_coszen;
    omega_0 = 0.5 * scat_leaf * (*proj_area) / tmp_0 * (1.0 - tmp_1 / tmp_0 * 
                                                log((tmp_1 + tmp_0) / tmp_1));
    coef_leaf_direct = (1.0 + tau_leaf_diffuse * tau_leaf_direct) / 
                    (scat_leaf * tau_leaf_diffuse * tau_leaf_direct) * omega_0;
    coef_leaf_diffuse = 0.5 * (ReflectVeg[lidx] + TransmitVeg[lidx] + (ReflectVeg[lidx] - 
                        TransmitVeg[lidx]) * pow((1.0 + LOI) / 2.0, 2.0)) / scat_leaf;

    /* adjust omega, betad, and betai for intercepted snow */
    if (Tcanopy > 0.0) {
        tmp_0 = scat_leaf;
        tmp_1 = coef_leaf_direct;
        tmp_2 = coef_leaf_diffuse;
    }
    else {
        tmp_0 = (1.0 - wetFrac) * scat_leaf + wetFrac * param.SNOW_OMEGAS[lidx];
        tmp_1 = ((1.0 - wetFrac) * scat_leaf * coef_leaf_direct + wetFrac * 
                        param.SNOW_OMEGAS[lidx] * coef_leaf_direct) / tmp_0;
        tmp_2 = ((1.0 - wetFrac) * scat_leaf * coef_leaf_diffuse + wetFrac * 
                        param.SNOW_OMEGAS[lidx] * coef_leaf_diffuse) / tmp_0;
    }
    scat_canopy = tmp_0;
    coef_can_direct = tmp_1;
    coef_can_diffuse = tmp_2;

    /* absorbed, reflected, transmitted fluxes per unit incoming radiation */
    B = 1.0 - scat_canopy + scat_canopy * coef_can_diffuse;
    C = scat_canopy * coef_can_diffuse;
    tmp_0 = tau_leaf_diffuse * tau_leaf_direct;
    D = tmp_0 * scat_canopy * coef_can_direct;
    F = tmp_0 * scat_canopy * (1.0 - coef_can_direct);
    tmp_1 = B * B - C * C;
    H = sqrt(tmp_1) / tau_leaf_diffuse;
    sigma = tmp_0 * tmp_0 - tmp_1;
    if (fabs(sigma) < param.TOL_A) {
        sigma = sign(param.TOL_A, sigma);
    }
    P1 = B + tau_leaf_diffuse * H;
    P2 = B - tau_leaf_diffuse * H;
    P3 = B + tmp_0;
    P4 = B - tmp_0;
    S1 = exp(-H * (NetLAI + NetSAI));
    S2 = exp(-tau_leaf_direct * (NetLAI + NetSAI));
    // direct
    if (index == 0) {
        U1 = B - C / AlbedoGrndDir[lidx];
        U2 = B - C * AlbedoGrndDir[lidx];
        U3 = F + C * AlbedoGrndDir[lidx];
    }
    else {  // diffuse
        U1 = B - C / AlbedoGrndDfs[lidx];
        U2 = B - C * AlbedoGrndDfs[lidx];
        U3 = F + C * AlbedoGrndDfs[lidx];                        
    }
    tmp_2 = U1 - tau_leaf_diffuse * H;
    tmp_3 = U1 + tau_leaf_diffuse * H;
    D1 = P1 * tmp_2 / S1 - P2 * tmp_3 * S1;
    tmp_4 = U2 + tau_leaf_diffuse * H;
    tmp_5 = U2 - tau_leaf_diffuse * H;
    D2 = tmp_4 / S1 - tmp_5 * S1;
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
    H7 = (C * tmp_2) / (D1 * S1);
    H8 = (-C * tmp_3 * S1) / D1;
    H9 = tmp_4 / (D2 * S1);
    H10 = (-tmp_5 * S1) / D2;

    /* direct and diffuse fluxes below vegetation Niu and Yang (2004), JGR. */
    if (index == 0) {
        transmit_sw_direct = S2 * (1.0 - gap_frac_direct) + gap_frac_direct;
        transmit_sw_diffuse = (H4 * S2 / sigma + H5 * S1 + H6 / S1) * (1.0 - gap_frac_direct);
    }
    else {
        transmit_sw_direct = 0.0;
        transmit_sw_diffuse = (H9 * S1 + H10 / S1) * (1.0 - gap_frac_diffuse) + gap_frac_diffuse;
    }
    if (index == 0) {
        ShortDir2Dir[lidx] = transmit_sw_direct;
        ShortDfs2Dir[lidx] = transmit_sw_diffuse;
    }
    else {
        ShortDir2Dfs[lidx] = transmit_sw_direct;
        ShortDfs2Dfs[lidx] = transmit_sw_diffuse;
    }

    /* flux reflected by the surface (veg and grnd) */
    if (index == 0) {
        refl_total = (H1 / sigma + H2 + H3) * (1.0 - gap_frac_direct) + 
                                    AlbedoGrndDir[lidx] * gap_frac_direct;
        refl_over = (H1 / sigma + H2 + H3) * (1.0 - gap_frac_direct);
        refl_grnd = AlbedoGrndDir[lidx] * gap_frac_direct;
    }
    else {
        refl_total = (H7 + H8) * (1.0 - gap_frac_diffuse) + AlbedoGrndDfs[lidx] * gap_frac_diffuse;
        refl_over = (H7 + H8) * (1.0 - gap_frac_diffuse) + AlbedoGrndDfs[lidx] * gap_frac_diffuse;
        refl_grnd = 0.0;                        
    }
    if (index == 0) {
        AlbedoSurfDir[lidx] = refl_total;
        ReflSubDir[lidx] = refl_over;
        ReflGrndDir[lidx] = refl_grnd;
    }
    else {
        AlbedoSurfDfs[lidx] = refl_total;
        ReflSubDfs[lidx] = refl_over;
        ReflGrndDfs[lidx] = refl_grnd;                        
    }

    /* flux absorbed by vegetation */
    if (index == 0) {
        AbsSubDir[lidx] = 1.0 - AlbedoSurfDir[lidx] - (1.0 - AlbedoGrndDir[lidx]) * 
                    ShortDir2Dir[lidx] - (1.0 - AlbedoGrndDfs[lidx]) * ShortDfs2Dir[lidx];
    }
    else {
        AbsSubDfs[lidx] = 1.0 - AlbedoSurfDfs[lidx] - (1.0 - AlbedoGrndDir[lidx]) *
                    ShortDir2Dfs[lidx] - (1.0 - AlbedoGrndDfs[lidx]) * ShortDfs2Dfs[lidx];
    }

}

