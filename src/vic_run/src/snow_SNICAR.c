/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine calculate albedo of snow containing impurities and the 
 * evolution of snow effective radius.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief  
 * Compute the albedo of snow containing impurities and the evolution of 
 * snow effective radius.
 *****************************************************************************/
int
snow_SNICAR(size_t             comp_type,
            double             step_dt,
            double             coszen,
            double             snowfall,
            double            *mass_cnc_aer,
            energy_bal_struct *energy,
            cell_data_struct  *cell,
            snow_data_struct  *snow,
			soil_con_struct   *soil_con)
{
    extern parameters_struct param;
    double zc_wave[5] = {0.5, 0.85, 1.1, 1.35, 3.25};
    // 高斯积分点（弧度）
    static const double difgaus_pt[8] = {
        0.9894009, 0.9445750, 0.8656312, 0.7554044,
        0.6178762, 0.4580168, 0.2816036, 0.0950125
    };
    // 高斯权重
    static const double difgaus_wt[8] = {
        0.0271525, 0.0622535, 0.0951585, 0.1246290,
        0.1495960, 0.1691565, 0.1826034, 0.1894506
    };
    // 非球形雪粒不对称因子参数化系数
    static const double g_wvl[8] = {
        0.25, 0.70, 1.41, 1.90, 2.50, 3.50, 4.00, 5.00
    };
    static const double g_b0[7] = {
        9.76029E-1, 9.67798E-1, 1.00111, 1.00224,
        9.64295E-1, 9.97475E-1, 9.97475E-1
    };
    static const double g_b1[7] = {
        5.21042E-1, 4.96181E-1, 1.83711E-1, 1.37082E-1,
        5.50598E-2, 8.48743E-2, 8.48743E-2
    };
    static const double g_b2[7] = {
        -2.66792E-4, 1.14088E-3, 2.37011E-4, -2.35905E-4,
        8.40449E-4, -4.71484E-4, -4.71484E-4
    };
    double *ss_alb_dir;
    double *ss_alb_dfs;
    double *asym_snow_dir;
    double *asym_snow_dfs;
    double *mass_ext_dir;
    double *mass_ext_dfs;
    static const double ss_alb_bcphil_dir[5] = {0.758058, 0.708629, 0.654803, 0.599902, 0.415976};
    static const double ss_alb_bcphil_dfs[5] = {0.758067, 0.709739, 0.655509, 0.605827, 0.482284};
    static const double ss_alb_bcphob_dfs[5] = {0.366232, 0.301432, 0.250557, 0.214858, 0.151357};
    static const double ss_alb_bcphob_dir[5] = {0.366239, 0.300250, 0.250005, 0.211183, 0.132030};
    static const double ss_alb_dust1_dir[5] = {0.982610, 0.993620, 0.993858, 0.993084, 0.983789};
    static const double ss_alb_dust1_dfs[5] = {0.982676, 0.993606, 0.993863, 0.993282, 0.988931};
    static const double ss_alb_dust2_dir[5] = {0.951550, 0.981651, 0.989356, 0.992101, 0.992631};
    static const double ss_alb_dust2_dfs[5] = {0.951619, 0.981442, 0.989298, 0.992030, 0.993126};
    static const double ss_alb_dust3_dir[5] = {0.916266, 0.960250, 0.969149, 0.973407, 0.981918};
    static const double ss_alb_dust3_dfs[5] = {0.916129, 0.960537, 0.969218, 0.973622, 0.984143};
    static const double ss_alb_dust4_dir[5] = {0.865081, 0.934038, 0.950694, 0.955140, 0.966210};
    static const double ss_alb_dust4_dfs[5] = {0.865239, 0.933621, 0.950604, 0.954943, 0.963697};
    static const double ss_alb_dust5_dir[5] = {0.785926, 0.886459, 0.913632, 0.925360, 0.940251};
    static const double ss_alb_dust5_dfs[5] = {0.786042, 0.885755, 0.913356, 0.925003, 0.936343};
    static const double asm_prm_bcphob_dir[5] = {0.440809, 0.333372, 0.274836, 0.235804, 0.160655};
    static const double asm_prm_bcphob_dfs[5] = {0.440763, 0.334865, 0.275412, 0.239339, 0.179417};
    static const double ss_alb_ocphil_dir[5] = {0.963861, 0.999168, 0.999052, 0.998544, 0.854790};
    static const double ss_alb_ocphil_dfs[5] = {0.963963, 0.999169, 0.999054, 0.998668, 0.977566};
    static const double asm_prm_ocphil_dir[5] = {0.694844, 0.643041, 0.590617, 0.545134, 0.419638};
    static const double asm_prm_ocphil_dfs[5] = {0.694939, 0.644230, 0.591232, 0.549758, 0.459431};
    static const double ext_cff_mss_ocphil_dir[5] = {50107.07, 24353.02, 13745.27, 8713.34, 3646.99};
    static const double ext_cff_mss_ocphil_dfs[5] = {50099.28, 24654.86, 13829.43, 9088.47, 4048.98};
    static const double ss_alb_ocphob_dir[5] = {0.773160, 0.990013, 0.986363, 0.982205, 0.964227};
    static const double ss_alb_ocphob_dfs[5] = {0.773243, 0.990086, 0.986413, 0.982688, 0.970958};
    static const double asm_prm_ocphob_dir[5] = {0.579585, 0.461883, 0.387327, 0.333355, 0.223369};
    static const double asm_prm_ocphob_dfs[5] = {0.579570, 0.463703, 0.388107, 0.338380, 0.251348};
    static const double ext_cff_mss_ocphob_dir[5] = {4928.72, 1350.66, 611.17, 344.96, 111.45};
    static const double ext_cff_mss_ocphob_dfs[5] = {4922.89, 1374.15, 615.92, 362.70, 139.92};
    static const double asm_prm_bcphil_dir[5] = {0.645843, 0.576849, 0.516039, 0.464377, 0.337056};
    static const double asm_prm_bcphil_dfs[5] = {0.645911, 0.578235, 0.516745, 0.469552, 0.372236};
    static const double ext_cff_mss_bcphob_dir[5] = {12376.36, 7867.52, 5678.42, 4494.98, 2780.40};
    static const double ext_cff_mss_bcphob_dfs[5] = {12376.82, 7927.52, 5697.23, 4591.12, 3133.56};
    static const double ext_cff_mss_bcphil_dir[5] = {49947.83, 27984.92, 18265.60, 13249.65, 7306.51};
    static const double ext_cff_mss_bcphil_dfs[5] = {49937.52, 28253.70, 18347.13, 13642.25, 8042.38};
    static const double asm_prm_dust1_dir[5] = {0.677934, 0.707388, 0.662214, 0.601745, 0.410256};
    static const double asm_prm_dust1_dfs[5] = {0.677869, 0.708229, 0.663041, 0.606302, 0.480786};
    static const double asm_prm_dust2_dir[5] = {0.695228, 0.637981, 0.672995, 0.705675, 0.681674};
    static const double asm_prm_dust2_dfs[5] = {0.695183, 0.637701, 0.672411, 0.704000, 0.706122};
    static const double asm_prm_dust3_dir[5] = {0.778833, 0.746215, 0.683282, 0.620030, 0.643931};
    static const double asm_prm_dust3_dfs[5] = {0.778799, 0.746797, 0.684328, 0.623973, 0.623736};
    static const double asm_prm_dust4_dir[5] = {0.821899, 0.781644, 0.764708, 0.746967, 0.706380};
    static const double asm_prm_dust4_dfs[5] = {0.821856, 0.782220, 0.765058, 0.746750, 0.728913};
    static const double asm_prm_dust5_dir[5] = {0.861434, 0.821163, 0.804414, 0.793842, 0.767681};
    static const double asm_prm_dust5_dfs[5] = {0.861395, 0.821530, 0.804506, 0.793945, 0.775716};
    static const double ext_cff_mss_dust1_dir[5] = {2583.01, 2422.90, 1665.01, 1140.80, 447.84};
    static const double ext_cff_mss_dust1_dfs[5] = {2583.66, 2440.83, 1673.37, 1181.23, 556.93};
    static const double ext_cff_mss_dust2_dir[5] = {814.37, 947.72, 1131.89, 1231.69, 968.79};
    static const double ext_cff_mss_dust2_dfs[5] = {814.22, 944.10, 1129.72, 1229.83, 1101.81};
    static const double ext_cff_mss_dust3_dir[5] = {375.88, 405.16, 393.46, 386.27, 497.39};
    static const double ext_cff_mss_dust3_dfs[5] = {375.91, 404.54, 393.81, 384.89, 467.81};
    static const double ext_cff_mss_dust4_dir[5] = {190.59, 196.58, 203.70, 199.13, 220.91};
    static const double ext_cff_mss_dust4_dfs[5] = {190.59, 196.49, 203.83, 198.14, 216.84};
    static const double ext_cff_mss_dust5_dir[5] = {92.79, 94.66, 95.94, 96.87, 100.66};
    static const double ext_cff_mss_dust5_dfs[5] = {92.79, 94.62, 95.88, 96.53, 98.94};
    // initialize
    size_t i, j, k, nidx;
    size_t min_mie_radius = 30;
    bool virtual_flag;
    bool DELTA_flag = true;
    bool solver_flag = true;
    double flux_dir;
    double flux_dfs;
    double tmp_pack_liq;
    double tmp_pack_ice;
    double tmp_radius;
    double tmp_albedo;
    double refk;
    double F_sfc_pls;
    double band_wgt[5];
    double ss_alb_aer[SNOW_NUM_AER];
    double asm_prm_aer[SNOW_NUM_AER];
    double ext_cff_mss_aer[SNOW_NUM_AER];
    double trndir[MAX_SNOWS+1];
    double trntdr[MAX_SNOWS+1];
    double trndfs[MAX_SNOWS+1];
    double rupdir[MAX_SNOWS+1];
    double rupdfs[MAX_SNOWS+1];
    double rdndfs[MAX_SNOWS+1];
    double dfdir[MAX_SNOWS+1];
    double dfdfs[MAX_SNOWS+1];
    double rdir[MAX_SNOWS+1];
    double rdfs_a[MAX_SNOWS+1];
    double rdfs_b[MAX_SNOWS+1];
    double tdir[MAX_SNOWS+1];
    double tdfs_a[MAX_SNOWS+1];
    double tdfs_b[MAX_SNOWS+1];
    double trnlay[MAX_SNOWS+1];
    double dftmp[MAX_SNOWS+1];
    double albedo[MAX_SWBANDS];
    double abs_flux[MAX_SWBANDS];
    double abs_flux_dir[SNICAR_BANDS];
    double abs_flux_dfs[SNICAR_BANDS];
    double albsfc_layer[SNICAR_BANDS];
    double *radius = snow->radius;
    double *pack_ice = snow->pack_ice;
    double *pack_liq = snow->pack_liq;
    double *AlbedoSoilDir = energy->AlbedoSoilDir;
    double *AlbedoSoilDfs = energy->AlbedoSoilDfs;
    double *AlbedoSnowDir = energy->AlbedoSnowDir;
    double *AlbedoSnowDfs = energy->AlbedoSnowDfs;
    static const double DIR_wgt[5] = {1.0, 0.471, 0.178, 0.139, 0.212};
    static const double DFS_wgt[5] = {1.0, 0.625, 0.190, 0.101, 0.084};

    // 检查是否需要进行辐射传输计算
    if (coszen > 0.0 && snow->swq > MIN_SNOWMASS) {
        size_t tmp_Nsnow = snow->Nsnow;
        if (tmp_Nsnow == 0) {
            virtual_flag = true;
            tmp_Nsnow = 1;
            tmp_pack_ice = snow->swq;
            tmp_pack_liq = 0.0;
            tmp_radius = param.SNOW_RADIUS_MIN;
        }
        else {
            virtual_flag = false;
            tmp_pack_ice = pack_ice[i];
        }
        // Set spectral underlying surface albedos to their corresponding VIS or NIR albedos
        if (cell->IS_VEG) { // soil tile
            if (comp_type == BAND_DIR) {
                for (i = 0; i < SNICAR_BANDS; i++) {
                    albsfc_layer[i] = AlbedoSoilDir[1];
                }
                albsfc_layer[0] = AlbedoSoilDir[0];
            }
            else if (comp_type == BAND_DFS) {
                for (i = 0; i < SNICAR_BANDS; i++) {
                    albsfc_layer[i] = AlbedoSoilDfs[1];
                }
                albsfc_layer[0] = AlbedoSoilDfs[0];              
            }
        }
        else if (cell->IS_GLAC) {   // glacier tile
            for (i = 0; i < SNICAR_BANDS; i++) {
                albsfc_layer[i] = param.LAKE_ALBEDO[1];
            }
            albsfc_layer[0] = param.LAKE_ALBEDO[0];      
        }
        else if (cell->IS_WET) {   // lake or wetland tile
            for (i = 0; i < SNICAR_BANDS; i++) {
                albsfc_layer[i] = param.GLAC_ALBEDO[1];
            }
            albsfc_layer[0] = param.GLAC_ALBEDO[0];             
        }

        // The following weights are appropriate for 
        // surface-incident flux in a mid-latitude winter atmosphere
        if (comp_type == BAND_VIS) {
            for (i = 0; i < SNICAR_BANDS; i++) {
                band_wgt[i] = DIR_wgt[i];
            }
        }
        else if (comp_type == BAND_NIR) {
            for (i = 0; i < SNICAR_BANDS; i++) {
                band_wgt[i] = DFS_wgt[i];
            }
        }
        double exp_min = exp(-10.0);
        double ss_alb_snow = 0.0;
        double asm_prm_snow = 0.0;
        double ext_cff_mss_snow = 0.0;
        double mu_not = max(coszen, 0.01);
        // Loop over snow spectral bands
        for (i = 0; i < SNICAR_BANDS; i++) {
            while (solver_flag) {
                if (comp_type == BAND_VIS) {
                    flux_dir = 1.0 / (mu_not * CONST_PI);
                    flux_dfs = 0.0;
                }
                else if (comp_type == BAND_NIR) {
                    flux_dir = 0.0;
                    flux_dfs = 1.0;
                }
                // Pre-emptive error handling: aerosols can reap havoc on these absorptive bands.
                if (i >= 3) {
                    for (k = 0; k < tmp_Nsnow; k++) {
                        for (j = 0; j < SNOW_NUM_AER; j++) {
                            nidx = k * j;
                            mass_cnc_aer[nidx] = 0.0;
                        }
                    }
                }
                // 设置光学性质
                for (j = 0; j < tmp_Nsnow; j++) {
                    nidx = radius[j] - min_mie_radius + 1;
                    if (comp_type == BAND_DIR) {                    
                        // snow optical properties (direct radiation)
                        ss_alb_snow = ss_alb_dir[nidx];
                        asm_prm_snow = asym_snow_dir[nidx];
                        ext_cff_mss_snow = mass_ext_dir[nidx];
                    }
                    else if (comp_type == BAND_DFS) {
                        // snow optical properties (diffuse radiation)
                        ss_alb_snow = ss_alb_dfs[nidx];
                        asm_prm_snow = asym_snow_dfs[nidx];
                        ext_cff_mss_snow = mass_ext_dfs[nidx];
                    }
                    asm_prm_snow = max(0.99, asm_prm_snow);
                }
                // aerosol species 2 optical properties, hydrophobic BC
                ss_alb_aer[1]      = ss_alb_bcphob_dfs[i];
                asm_prm_aer[1]     = asm_prm_bcphob_dfs[i];
                ext_cff_mss_aer[1] = ext_cff_mss_bcphob_dfs[i];
                // aerosol species 3 optical properties, hydrophilic OC
                ss_alb_aer[2]      = ss_alb_ocphil_dfs[i];
                asm_prm_aer[2]     = asm_prm_ocphil_dfs[i];
                ext_cff_mss_aer[2] = ext_cff_mss_ocphil_dfs[i];
                // aerosol species 4 optical properties, hydrophobic OC
                ss_alb_aer[3]      = ss_alb_ocphob_dfs[i];
                asm_prm_aer[3]     = asm_prm_ocphob_dfs[i];
                ext_cff_mss_aer[3] = ext_cff_mss_ocphob_dfs[i];
                // Optics for BC/dust-snow external mixing:
                // aerosol species 1 optical properties, hydrophilic BC
                ss_alb_aer[0]      = ss_alb_bcphil_dfs[i];
                asm_prm_aer[0]     = asm_prm_bcphil_dfs[i];
                ext_cff_mss_aer[0] = ext_cff_mss_bcphil_dfs[i];
                // aerosol species 5 optical properties, dust size1
                ss_alb_aer[4]      = ss_alb_dust1_dfs[i];
                asm_prm_aer[4]     = asm_prm_dust1_dfs[i];
                ext_cff_mss_aer[4] = ext_cff_mss_dust1_dfs[i];
                // aerosol species 6 optical properties, dust size2
                ss_alb_aer[5]      = ss_alb_dust2_dfs[i];
                asm_prm_aer[5]     = asm_prm_dust2_dfs[i];
                ext_cff_mss_aer[5] = ext_cff_mss_dust2_dfs[i];
                // aerosol species 7 optical properties, dust size3
                ss_alb_aer[6]      = ss_alb_dust3_dfs[i];
                asm_prm_aer[6]     = asm_prm_dust3_dfs[i];
                ext_cff_mss_aer[6] = ext_cff_mss_dust3_dfs[i];
                // aerosol species 8 optical properties, dust size4
                ss_alb_aer[7]      = ss_alb_dust4_dfs[i];
                asm_prm_aer[7]     = asm_prm_dust4_dfs[i];
                ext_cff_mss_aer[7] = ext_cff_mss_dust4_dfs[i];
                double tau = 0.0, omega = 0.0, g = 0.0;
                // 计算雪的光学厚度
                for (j = 0; j < tmp_Nsnow; j++) {
                    double L_aer = 0.0, tau_aer = 0.0;
                    double tau_sum = 0.0, omega_sum = 0.0, g_sum = 0.0;
                    double L_snow = pack_ice[j] + pack_liq[j];
                    double tau_snow = L_snow * ext_cff_mss_snow;
                    // 计算气溶胶的光学贡献
                    for (j = 0; j < SNOW_NUM_AER; j++) {
                        L_aer = L_snow * mass_cnc_aer[i];
                        tau_aer = L_aer * ext_cff_mss_aer[i];
                        tau_sum += tau_aer;
                        omega_sum += tau_aer * ss_alb_aer[i];
                        g_sum += tau_aer * ss_alb_aer[i] * asm_prm_aer[i];
                    }
                    tau = tau_sum + tau_snow;
                    omega = (1 / tau) * (omega_sum + (ss_alb_snow * tau_snow));
                    g = (1 / (tau * omega)) * (g_sum + (asm_prm_snow * ss_alb_snow * tau_snow));
                }
                double g_star = 0.0, omega_star = 0.0, tau_star = 0.0;
                if (DELTA_flag) {
                    for (j = 0; j < tmp_Nsnow; j++) {
                        g_star = g / (1 + g);
                        omega_star = (1.0 - g * g) * omega / (1.0 - omega * (g * g));
                        tau_star = (1.0 - omega * (g * g)) * tau;
                    }
                }
                else {
                    for (j = 0; j < tmp_Nsnow; j++) {
                        g_star = g;
                        omega_star = omega;
                        tau_star = tau;
                    }
                }
                // Start Adding-doubling RT solver
                for (j = 0; j <= tmp_Nsnow; j++) {
                    trndir[j] = 0.0;
                    trntdr[j] = 0.0;
                    trndfs[j] = 0.0;
                    rupdir[j] = 0.0;
                    rupdfs[j] = 0.0;
                    rdndfs[j] = 0.0;
                }
                // initialize top interface of top layer
                trndir[0] = 1.0;
                trntdr[0] = 1.0;
                trndfs[0] = 1.0;
                rdndfs[0] = 0.0;
                // 1.snow and aerosol layer column mass (L_snw, L_aer [kg/m^2])
                // 2.optical Depths (tau_snw, tau_aer)
                // 3.weighted Mie properties (tau, omega, g)
                for (j = 0; j < tmp_Nsnow; j++) {
                    if (zc_wave[i] <= 1.2) {
                        // BC-snow internal mixing applied to hydrophilic BC if activated
                        // BC-snow internal mixing primarily affect snow single-scattering albedo

                    }
                }

                // begin main level loop for snow layer interfaces except for the very bottom
                for (j = 0; j < tmp_Nsnow; j++) {
                    rdir[j] = 0.0;
                    rdfs_a[j] = 0.0;
                    rdfs_b[j] = 0.0;
                    tdir[j] = 0.0;
                    tdfs_a[j] = 0.0;
                    tdfs_b[j] = 0.0;
                    trnlay[j] = 0.0;
                    if (trntdr[j] > 0.001) {
                        // Delta变换后的单次散射性质
                        double ts = tau_star;
                        double ws = omega_star;
                        double gs = g_star;
                        // Delta-Eddington解
                        double lm = sqrt(3.0 * (1.0 - ws) * (1.0 - ws * gs));
                        double ue = 1.5 * (1.0 - ws * gs) / lm;
                        double extins = max(exp_min, exp(-lm * ts));
                        double ne = ((ue + 1.0) * (ue + 1.0) / extins) - 
                                    ((ue - 1.0) * (ue - 1.0) * extins);
                        // 漫射反射率和透射率
                        rdfs_a[j] = (ue * ue - 1.0) * (1.0 / extins - extins) / ne;
                        tdfs_a[j] = 4.0 * ue / ne;
                        // 直接辐射处理
                        trnlay[j] = max(exp_min, exp(-ts / mu_not));
                        
                        double alp = 0.75 * ws * mu_not * 
                                    ((1.0 + gs * (1.0 - ws)) / (1.0 - lm * lm * mu_not * mu_not));
                        double gam = 0.5 * ws * 
                                    ((1.0 + 3.0 * gs * (1.0 - ws) * mu_not * mu_not) / 
                                    (1.0 - lm * lm * mu_not * mu_not));
                        double apg = alp + gam;
                        double amg = alp - gam;
                        
                        rdir[j] = apg * rdfs_a[j] + amg * (tdfs_a[j] * trnlay[j] - 1.0);
                        tdir[j] = apg * tdfs_a[j] + 
                                (amg * rdfs_a[j] - apg + 1.0) * trnlay[j];
                        
                        // 高斯积分修正漫射性质
                        double R1 = rdfs_a[j];
                        double T1 = tdfs_a[j];
                        double swt = 0.0, smr = 0.0, smt = 0.0;
                        for (k = 0; k < MAX_GAUSSIAN; k++) {
                            double mu = difgaus_pt[k];
                            double gwt = difgaus_wt[k];
                            swt += mu * gwt;
                            
                            double trn = max(exp_min, exp(-ts / mu));
                            alp = 0.75 * ws * mu * 
                                ((1.0 + gs * (1.0 - ws)) / (1.0 - lm * lm * mu * mu));
                            gam = 0.5 * ws * 
                                ((1.0 + 3.0 * gs * (1.0 - ws) * mu * mu) / 
                                (1.0 - lm * lm * mu * mu));
                            apg = alp + gam;
                            amg = alp - gam;
                            
                            double rdr_temp = apg * R1 + amg * T1 * trn - amg;
                            double tdr_temp = apg * T1 + amg * R1 * trn - apg * trn + trn;
                            
                            smr += mu * rdr_temp * gwt;
                            smt += mu * tdr_temp * gwt;
                        }
                        rdfs_a[j] = smr / swt;
                        tdfs_a[j] = smt / swt;
                        rdfs_b[j] = rdfs_a[j];
                        tdfs_b[j] = tdfs_a[j];
                    } // trntdr(i) > trmin
                    // 计算下一界面的通量
                    trndir[j+1] = trndir[j] * trnlay[j];
                    
                    double refkm1 = 1.0 / (1.0 - rdndfs[j] * rdfs_a[j]);
                    double tdrrdir = trndir[j] * rdir[j];
                    double tdndif = trntdr[j] - trndir[j];
                    
                    trntdr[j+1] = trndir[j] * tdir[j] + 
                                (tdndif + tdrrdir * rdndfs[j]) * refkm1 * tdfs_a[j];
                    rdndfs[j+1] = rdfs_b[j] + 
                                (tdfs_b[j] * rdndfs[j] * refkm1 * tdfs_a[j]);
                    trndfs[j+1] = trndfs[j] * refkm1 * tdfs_a[j];
                }
                // 从下到上计算反射率
                if (cell->IS_VEG) {
                    if (comp_type == BAND_VIS) {
                        if (i == 0) {
                            rupdir[tmp_Nsnow] = AlbedoSoilDir[0];
                            rupdfs[tmp_Nsnow] = AlbedoSoilDir[0];    
                        }
                        rupdir[tmp_Nsnow] = AlbedoSoilDir[1];
                        rupdfs[tmp_Nsnow] = AlbedoSoilDir[1];                   
                    }
                    else if (comp_type == BAND_NIR) {
                        if (i == 0) {
                            rupdir[tmp_Nsnow] = AlbedoSoilDfs[0];
                            rupdfs[tmp_Nsnow] = AlbedoSoilDfs[0];                          
                        }
                        rupdir[tmp_Nsnow] = AlbedoSoilDfs[1];
                        rupdfs[tmp_Nsnow] = AlbedoSoilDfs[1];
                    }
                }
                else if (cell->IS_GLAC) {
                    if (i == 0) {
                        rupdir[tmp_Nsnow] = param.GLAC_ALBEDO[0];
                        rupdfs[tmp_Nsnow] = param.GLAC_ALBEDO[0];                          
                    }
                    rupdir[tmp_Nsnow] = param.GLAC_ALBEDO[1];
                    rupdfs[tmp_Nsnow] = param.GLAC_ALBEDO[1];                   
                }
                else if (cell->IS_WET) {
                    if (i == 0) {
                        rupdir[tmp_Nsnow] = param.LAKE_ALBEDO[0];
                        rupdfs[tmp_Nsnow] = param.LAKE_ALBEDO[0];                          
                    }
                    rupdir[tmp_Nsnow] = param.LAKE_ALBEDO[1];
                    rupdfs[tmp_Nsnow] = param.LAKE_ALBEDO[1];                    
                }

                for (j = tmp_Nsnow - 1; j >= 0; j--) {
                    double refkp1 = 1.0 / (1.0 - rdfs_b[j] * rupdfs[j+1]);
                    
                    rupdir[j] = rdir[j] + 
                                (trnlay[j] * rupdir[j+1] + 
                                (tdir[j] - trnlay[j]) * rupdfs[j+1]) * refkp1 * tdfs_b[j];
                    rupdfs[j] = rdfs_a[j] + 
                                tdfs_a[j] * rupdfs[j+1] * refkp1 * tdfs_b[j];
                }
                // 计算界面净通量
                for (j = 0; j <= tmp_Nsnow; j++) {
                    refk = 1.0 / (1.0 - rdndfs[j] * rupdfs[j]);
                    // netflux, down - up
                    dfdir[j] = trndir[j] + 
                            (trntdr[j] - trndir[j]) * (1.0 - rupdfs[j]) * refk -
                            trndir[j] * rupdir[j] * (1.0 - rdndfs[j]) * refk;
                    
                    if (dfdir[j] < 1e-11) {
                        dfdir[j] = 0.0;
                    }
                    // dfdif = fdifdn - fdifup
                    dfdfs[j] = trndfs[j] * (1.0 - rupdfs[j]) * refk;
                    if (dfdfs[j] < 1e-11) {
                        dfdfs[j] = 0.0;
                    }
                }
                // SNICAR_AD_RT is called twice for direct and diffuse incident fluxes
                if (comp_type == BAND_DIR) {
                    tmp_albedo = rupdir[0];
                    refk = 1.0 / (1.0 - rdndfs[0] * rupdfs[0]);
                    F_sfc_pls = (trndir[0] * rupdir[0] + (trntdr[0] - trndir[0]) * rupdfs[0]) * refk;
                    for (j = 0; j <= tmp_Nsnow; j++) {
                        dftmp[j] = dfdir[j];
                    }
                }
                else if (comp_type == BAND_DFS) {
                    tmp_albedo = rupdfs[0];
                    refk = 1.0 / (1.0 - rdndfs[0] * rupdfs[0]);
                    F_sfc_pls = trndfs[0] * rupdfs[0] * refk;
                    for (j = 0; j <= tmp_Nsnow; j++) {
                        dftmp[j] = dfdfs[j];
                    }
                }
                // Absorbed flux in each layer
                for (j = 0; j < tmp_Nsnow; j++) {
                    abs_flux[nidx] = dftmp[j] - dftmp[j+1];
                    if (abs_flux[nidx] > 0.0001) {
                        log_err("Error in snow SNICAR: negative absoption(%.4f)", abs_flux[nidx]);
                    }
                }
                // If there are no snow layers (but still snow)
                if (virtual_flag) {
                    abs_flux[0] = 0.0;
                    abs_flux[1] = 0.0;
                }
                for (j = 0; j < tmp_Nsnow; j++) {
                    nidx = j * SNICAR_BANDS; 
                    if (abs_flux[nidx] < 0.0) {
                        abs_flux[nidx] = 0.0;
                    }
                }
                // no need to repeat calculations for adding-doubling solver
                solver_flag = false;
            }
            albedo[i] = tmp_albedo;
            // Check that albedo is less than 1
            if (albedo[i] > 1.0) {
                log_err("Error in snow SNICAR: albedo[%d](%.4f) > 1.0", albedo[i], i);
            }
        }
        // 加权平均得到VIS和NIR反照率
        double albedo_sum = 0.0;
        if (comp_type == BAND_DIR) {
            AlbedoSnowDir[0] = albedo[0];
        }
        else if (comp_type == BAND_DFS) {
            AlbedoSnowDfs[0] = albedo[0];
        }
        for (k = 1; k < SNICAR_BANDS; k++) {
            albedo_sum += band_wgt[k] * albedo[k];
        }
        if (comp_type == BAND_DIR) {
            AlbedoSnowDir[1] = albedo_sum;
        }
        else if (comp_type == BAND_DFS) {
            AlbedoSnowDfs[1] = albedo_sum;
        }
        // 加权平均得到VIS和NIR吸收的能量
        double flux_sum = 0.0;
        if (comp_type == BAND_DIR) {
            for (j = 0; j < tmp_Nsnow; j++) {
                nidx = j * 1;
                abs_flux_dir[nidx] = abs_flux[nidx];
            }
        }
        else if (comp_type == BAND_DFS) {
            for (j = 0; j < tmp_Nsnow; j++) {
                nidx = j * 1;
                abs_flux_dfs[nidx] = abs_flux[nidx];
            }           
        }
        for (j = 0; j < tmp_Nsnow; j++) {
            for (k = 1; k < SNICAR_BANDS; k++) {
                nidx = j * SNICAR_BANDS + k;
                flux_sum += band_wgt[k] * abs_flux[nidx];
            }
            if (comp_type == BAND_DIR) {
                abs_flux_dir[nidx] = flux_sum;
            }
            else if (comp_type == BAND_DFS) {
                abs_flux_dfs[nidx] = flux_sum;
            }
        }
        // 太阳天顶角调整（高天顶角时的修正）
        if (mu_not < 0.2588 && comp_type == BAND_DIR) {  // cos(75°) = 0.2588
            double sza_c1 = 0.085730 - 0.630883 * mu_not + 1.303723 * mu_not * mu_not;
            double sza_c0 = 1.467291 - 3.338043 * mu_not + 6.807489 * mu_not * mu_not;
            double sza_factor = sza_c1 * (log10(radius[0]) - 6.0) + sza_c0;
            sza_factor = max(sza_factor, 0.0);
            
            double flux_sza_adjust = AlbedoSoilDir[1] * (sza_factor - 1.0);
            AlbedoSoilDir[1] *= sza_factor;
            
            // 调整顶层吸收通量
            abs_flux_dir[tmp_Nsnow] -= flux_sza_adjust;
        }
    }
    else if (coszen > 0.0 && snow->swq > 0.0) {
        if (cell->IS_VEG) {
            if (comp_type == BAND_DIR) {
                AlbedoSnowDir[0] = AlbedoSoilDir[0];
                AlbedoSnowDir[1] = AlbedoSoilDir[1];
            }
            else if (comp_type == BAND_DFS) {
                AlbedoSnowDfs[0] = AlbedoSoilDfs[0];
                AlbedoSnowDfs[1] = AlbedoSoilDfs[1];
            }
        }
        else if (cell->IS_GLAC) {
            AlbedoSnowDfs[0] = param.GLAC_ALBEDO[0];
            AlbedoSnowDfs[1] = param.GLAC_ALBEDO[1];
        }
        else if (cell->IS_WET) {
            AlbedoSnowDfs[0] = param.LAKE_ALBEDO[0];
            AlbedoSnowDfs[1] = param.LAKE_ALBEDO[1];            
        }
    }
    else {
        if (comp_type == BAND_DIR) {
            AlbedoSnowDir[0] = 0.0;
            AlbedoSnowDir[1] = 0.0;         
        }
        else if (comp_type == BAND_DFS) {
            AlbedoSnowDir[0] = 0.0;
            AlbedoSnowDfs[1] = 0.0;
        }
    }

    return (0);
}