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
            double             coszen,
            energy_bal_struct *energy,
            cell_data_struct  *cell,
            snow_data_struct  *snow)
{
    extern parameters_struct param;
    extern option_struct options;
    extern optical_struct optical;
    static const double zc_wave[5] = {
                        0.5, 0.85, 1.1, 1.35, 3.25};
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
    static const double g_F07_c2[7] = {
        1.349959E-1,  1.115697E-1,  9.853958E-2,  5.557793E-2,
        -1.233493E-1,  0.0,          0.0
    };
    static const double g_F07_c1[7] = {
        -3.987320E-1, -3.723287E-1, -3.924784E-1, -3.259404E-1,
        4.429054E-2, -1.726586E-1, -1.726586E-1
    };
    static const double g_F07_c0[7] = {
        7.938904E-1, 8.030084E-1, 8.513932E-1, 8.692241E-1,
        7.085850E-1, 6.412701E-1, 6.412701E-1
    };
    static const double g_F07_p2[7] = {
        3.165543E-3, 2.014810E-3, 1.780838E-3, 6.987734E-4,
        -1.882932E-2, -2.277872E-2, -2.277872E-2
    };
    static const double g_F07_p1[7] = {
        1.140557E-1, 1.143152E-1, 1.143814E-1, 1.071238E-1,
        1.353873E-1, 1.914431E-1, 1.914431E-1
    };
    static const double g_F07_p0[7] = {
        5.292852E-1, 5.425909E-1, 5.601598E-1, 6.023407E-1,
        6.473899E-1, 4.634944E-1, 4.634944E-1
    };

    // initialize
    size_t i, j, k, igb, nidx;
    size_t min_mie_radius = 30;
    bool virtual_flag;
    bool DELTA_flag = true;
    bool solver_flag = true;
    bool UseAerosol_flag = false;
    double flux_dir;
    double flux_dfs;
    double tmp_pack_liq[MAX_SNOWS] = {0};
    double tmp_pack_ice[MAX_SNOWS] = {0};
    double tmp_radius[MAX_SNOWS] = {0};
    double tmp_albedo;
    double refk;
    double F_sfc_pls;
    double band_wgt[SNICAR_BANDS] = {0};
    double ss_alb_aer[SNOW_NUM_AER] = {0};
    double asm_prm_aer[SNOW_NUM_AER] = {0};
    double ext_cff_mss_aer[SNOW_NUM_AER] = {0};
    double trndir[MAX_SNOWS+1] = {0};
    double trntdr[MAX_SNOWS+1] = {0};
    double trndfs[MAX_SNOWS+1] = {0};
    double rupdir[MAX_SNOWS+1] = {0};
    double rupdfs[MAX_SNOWS+1] = {0};
    double rdndfs[MAX_SNOWS+1] = {0};
    double dfdir[MAX_SNOWS+1] = {0};
    double dfdfs[MAX_SNOWS+1] = {0};
    double rdir[MAX_SNOWS+1] = {0};
    double rdfs_a[MAX_SNOWS+1] = {0};
    double rdfs_b[MAX_SNOWS+1] = {0};
    double tdir[MAX_SNOWS+1] = {0};
    double tdfs_a[MAX_SNOWS+1] = {0};
    double tdfs_b[MAX_SNOWS+1] = {0};
    double trnlay[MAX_SNOWS+1] = {0};
    double dftmp[MAX_SNOWS+1] = {0};
    double albedo[SNICAR_BANDS] = {0};
    double g_ice_Cg_tmp[7] = {0};
    double g_wvl_ct[7] = {0};
    double gg_ice_F07_tmp[7] = {0};
    double abs_flux[MAX_SNOWS+1][MAX_SWBANDS] = {0};
    double albsfc_layer[SNICAR_BANDS] = {0};
    double mass_cnc_aer[MAX_SNOWS][SNOW_NUM_AER] = {0}; // 初始化为0
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
            tmp_pack_ice[0] = snow->swq;
            tmp_pack_liq[0] = 0.0;
            tmp_radius[0] = param.SNOW_RADIUS_MIN;
        }
        else {
            virtual_flag = false;
            tmp_Nsnow = snow->Nsnow;
            for (j = 0; j < tmp_Nsnow; j++) {
                tmp_pack_ice[j] = pack_ice[j];
                tmp_pack_liq[j] = pack_liq[j];
                tmp_radius[j] = radius[j];
            }
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
        // initialize for nonspherical snow grains
        for (i = 0; i < 7; i++) {
            g_wvl_ct[i] = g_wvl[i+1] * 0.5 + g_wvl[i] * 0.5;
        }
        double exp_min = exp(-10.0);
        double diam_ice = 0.0;
        double fs_shape = 0.0;
        double AR_tmp = 0.0;
        double mu_not = max(coszen, 0.01);
        double g_star[MAX_SNOWS] = {0};
        double omega_star[MAX_SNOWS] = {0};
        double tau_star[MAX_SNOWS] = {0};
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
                if (i >= 3 && UseAerosol_flag) {
                    memset(mass_cnc_aer, 0, sizeof(mass_cnc_aer));
                }
                // aerosol species 2 optical properties, hydrophobic BC
                ss_alb_aer[1]      = optical.ss_alb_bcphob_dfs[i];
                asm_prm_aer[1]     = optical.asm_prm_bcphob_dfs[i];
                ext_cff_mss_aer[1] = optical.ext_cff_mss_bcphob_dfs[i];
                // aerosol species 3 optical properties, hydrophilic OC
                ss_alb_aer[2]      = optical.ss_alb_ocphil_dfs[i];
                asm_prm_aer[2]     = optical.asm_prm_ocphil_dfs[i];
                ext_cff_mss_aer[2] = optical.ext_cff_mss_ocphil_dfs[i];
                // aerosol species 4 optical properties, hydrophobic OC
                ss_alb_aer[3]      = optical.ss_alb_ocphob_dfs[i];
                asm_prm_aer[3]     = optical.asm_prm_ocphob_dfs[i];
                ext_cff_mss_aer[3] = optical.ext_cff_mss_ocphob_dfs[i];
                // Optics for BC/dust-snow external mixing:
                // aerosol species 1 optical properties, hydrophilic BC
                ss_alb_aer[0]      = optical.ss_alb_bcphil_dfs[i];
                asm_prm_aer[0]     = optical.asm_prm_bcphil_dfs[i];
                ext_cff_mss_aer[0] = optical.ext_cff_mss_bcphil_dfs[i];
                // aerosol species 5 optical properties, dust size1
                ss_alb_aer[4]      = optical.ss_alb_dust1_dfs[i];
                asm_prm_aer[4]     = optical.asm_prm_dust1_dfs[i];
                ext_cff_mss_aer[4] = optical.ext_cff_mss_dust1_dfs[i];
                // aerosol species 6 optical properties, dust size2
                ss_alb_aer[5]      = optical.ss_alb_dust2_dfs[i];
                asm_prm_aer[5]     = optical.asm_prm_dust2_dfs[i];
                ext_cff_mss_aer[5] = optical.ext_cff_mss_dust2_dfs[i];
                // aerosol species 7 optical properties, dust size3
                ss_alb_aer[6]      = optical.ss_alb_dust3_dfs[i];
                asm_prm_aer[6]     = optical.asm_prm_dust3_dfs[i];
                ext_cff_mss_aer[6] = optical.ext_cff_mss_dust3_dfs[i];
                // aerosol species 8 optical properties, dust size4
                ss_alb_aer[7]      = optical.ss_alb_dust4_dfs[i];
                asm_prm_aer[7]     = optical.asm_prm_dust4_dfs[i];
                ext_cff_mss_aer[7] = optical.ext_cff_mss_dust4_dfs[i];
                // aerosol species 9 optical properties, dust size4
                ss_alb_aer[8]      = optical.ss_alb_dust5_dfs[i];
                asm_prm_aer[8]     = optical.asm_prm_dust5_dfs[i];
                ext_cff_mss_aer[8] = optical.ext_cff_mss_dust5_dfs[i];
                // 设置光学性质
                for (j = 0; j < tmp_Nsnow; j++) {
                    double ss_alb_snow = 0.0;
                    double asm_prm_snow = 0.0;
                    double ext_cff_mss_snow = 0.0;
                    nidx = (size_t) round(tmp_radius[j] - min_mie_radius);
                    nidx = min(max(nidx, 0), SNICAR_RADII - 1);
                    if (comp_type == BAND_DIR) {
                        // snow optical properties (direct radiation)
                        ss_alb_snow = optical.ss_alb_dir[i][nidx];
                        ext_cff_mss_snow = optical.mass_ext_dir[i][nidx];
                    }
                    else if (comp_type == BAND_DFS) {
                        // snow optical properties (diffuse radiation)
                        ss_alb_snow = optical.ss_alb_dfs[i][nidx];
                        ext_cff_mss_snow = optical.mass_ext_dfs[i][nidx];
                    }
                    
                    // 计算不对称因子（球形雪）
                    if (options.SNOW_SHAPE == SPHERE) {
                        if (comp_type == BAND_DIR) {
                            asm_prm_snow = optical.asym_snow_dir[i][nidx];
                        }
                        else if (comp_type == BAND_DFS) {
                            asm_prm_snow = optical.asym_snow_dfs[i][nidx];
                        }
                    }
                    else {
                        diam_ice = 2.0 * tmp_radius[j];
                        // 根据形状设置参数
                        if (options.SNOW_SHAPE == SPHEROID) {
                            fs_shape = 0.929;
                            AR_tmp = 0.5;
                            
                            for (igb = 0; igb < 7; igb++) {
                                g_ice_Cg_tmp[igb] = g_b0[igb] * 
                                    pow(fs_shape / 0.788, g_b1[igb]) * 
                                    pow(diam_ice, g_b2[igb]);
                                gg_ice_F07_tmp[igb] = g_F07_c0[igb] + 
                                    g_F07_c1[igb] * AR_tmp + 
                                    g_F07_c2[igb] * AR_tmp * AR_tmp;
                            }
                        }
                        else if (options.SNOW_SHAPE == HEXAGONAL) {
                            fs_shape = 0.788;
                            AR_tmp = 2.5;
                            
                            for (igb = 0; igb < 7; igb++) {
                                g_ice_Cg_tmp[igb] = g_b0[igb] * 
                                    pow(fs_shape / 0.788, g_b1[igb]) * 
                                    pow(diam_ice, g_b2[igb]);
                                gg_ice_F07_tmp[igb] = g_F07_p0[igb] + 
                                    g_F07_p1[igb] * log(AR_tmp) + 
                                    g_F07_p2[igb] * log(AR_tmp) * log(AR_tmp);
                            }
                        }
                        else if (options.SNOW_SHAPE == KOCH) {
                            diam_ice = 2.0 * tmp_radius[j] / 0.544;
                            fs_shape = 0.712;
                            AR_tmp = 2.5;
                            
                            for (igb = 0; igb < 7; igb++) {
                                g_ice_Cg_tmp[igb] = g_b0[igb] * 
                                    pow(fs_shape / 0.788, g_b1[igb]) * 
                                    pow(diam_ice, g_b2[igb]);
                                gg_ice_F07_tmp[igb] = g_F07_p0[igb] + 
                                    g_F07_p1[igb] * log(AR_tmp) + 
                                    g_F07_p2[igb] * log(AR_tmp) * log(AR_tmp);
                            }
                        }
                        else {
                            log_err("Unknown SNOW_SHAPE option");
                        }
                        // 插值到目标波长
                        double g_Cg_intp = piecewise_linear_interp(7, g_wvl_ct, g_ice_Cg_tmp, zc_wave[i]);
                        double gg_F07_intp = piecewise_linear_interp(7, g_wvl_ct, gg_ice_F07_tmp, zc_wave[i]);
                        
                        // Fu(2007) eq.2.2
                        double g_ice_F07 = gg_F07_intp + 0.5 * (1.0 - gg_F07_intp) / ss_alb_snow;
                        
                        // He et al. (2017) eq.6
                        asm_prm_snow = g_ice_F07 * g_Cg_intp;
                    }
                    asm_prm_snow = max(0.99, asm_prm_snow);

                    // 计算雪的光学厚度
                    double L_aer = 0.0, tau_aer = 0.0;
                    double tau_sum = 0.0, omega_sum = 0.0, g_sum = 0.0;
                    double L_snow = tmp_pack_ice[j] + tmp_pack_liq[j];
                    double tau_snow = L_snow * ext_cff_mss_snow;
                    // 计算气溶胶的光学贡献
                    for (k = 0; k < SNOW_NUM_AER; k++) {
                        L_aer = L_snow * mass_cnc_aer[i][k];
                        tau_aer = L_aer * ext_cff_mss_aer[k];
                        tau_sum += tau_aer;
                        omega_sum += tau_aer * ss_alb_aer[k];
                        g_sum += tau_aer * ss_alb_aer[k] * asm_prm_aer[k];
                    }
                    double tau = tau_sum + tau_snow;
                    double omega = (1 / tau) * (omega_sum + (ss_alb_snow * tau_snow));
                    double g = (1 / (tau * omega)) * (g_sum + (asm_prm_snow * ss_alb_snow * tau_snow));
                    if (DELTA_flag) {
                        for (j = 0; j < tmp_Nsnow; j++) {
                            g_star[j] = g / (1 + g);
                            omega_star[j] = (1.0 - g * g) * omega / (1.0 - omega * (g * g));
                            tau_star[j] = (1.0 - omega * (g * g)) * tau;
                        }
                    }
                    else {
                        for (j = 0; j < tmp_Nsnow; j++) {
                            g_star[j] = g;
                            omega_star[j] = omega;
                            tau_star[j] = tau;
                        }
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
                        double ts = tau_star[j];
                        double ws = omega_star[j];
                        double gs = g_star[j];
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

                for (k = tmp_Nsnow; k > 0; k--) {
                    j = k - 1; // // 防止j下溢到SIZE_MAX
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
                    abs_flux[j][i] = dftmp[j] - dftmp[j+1];
                    if (abs_flux[j][i] < -0.0001) {
                        log_err("Error in snow SNICAR: negative absoption(%.4f)", abs_flux[j][i]);
                    }
                }
                // If there are no snow layers (but still snow)
                if (virtual_flag) {
                    abs_flux[0][i] = dftmp[0] - dftmp[1];
                    abs_flux[tmp_Nsnow][i] = dftmp[tmp_Nsnow];
                }
                for (j = 0; j < tmp_Nsnow; j++) { 
                    if (abs_flux[j][i] < 0.0) {
                        abs_flux[j][i] = 0.0;
                    }
                }
                // no need to repeat calculations for adding-doubling solver
                solver_flag = false;
            }
            // Energy conservation check:
            double F_abs_sum = 0.0;
            for (j = 0; j < tmp_Nsnow; j++) {
                F_abs_sum += abs_flux[j][i];
            }
            double energy_err = (mu_not * CONST_PI * flux_dir) + flux_dfs - 
                                    F_abs_sum - dftmp[tmp_Nsnow] - F_sfc_pls;
            if (fabs(energy_err) > 0.0001) {
                log_err("Error in snow SNICAR: energy not conserved(%.4f)", energy_err);
            }
            albedo[i] = tmp_albedo;
            // Check that albedo is less than 1
            if (albedo[i] > 1.0) {
                log_err("Error in snow SNICAR: albedo[%zu](%.4f) > 1.0", i, albedo[i]);
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
                energy->AbsShortDir[j][i] = abs_flux[j][i];
            }
        }
        else if (comp_type == BAND_DFS) {
            for (j = 0; j < tmp_Nsnow; j++) {
                energy->AbsShortDfs[j][i] = abs_flux[j][i];
            }           
        }
        for (j = 0; j < tmp_Nsnow; j++) {
            for (k = 1; k < SNICAR_BANDS; k++) {
                flux_sum += band_wgt[k] * abs_flux[j][k];
            }
            if (comp_type == BAND_DIR) {
                energy->AbsShortDir[j][i] = flux_sum;
            }
            else if (comp_type == BAND_DFS) {
                energy->AbsShortDfs[j][i] = flux_sum;
            }
        }
        // 太阳天顶角调整（高天顶角时的修正）
        if (mu_not < 0.2588 && comp_type == BAND_DIR) {  // cos(75°) = 0.2588
            double sza_c1 = 0.085730 - 0.630883 * mu_not + 1.303723 * mu_not * mu_not;
            double sza_c0 = 1.467291 - 3.338043 * mu_not + 6.807489 * mu_not * mu_not;
            double sza_factor = sza_c1 * (log10(tmp_radius[0]) - 6.0) + sza_c0;
            sza_factor = max(sza_factor, 0.0);
            
            double flux_sza_adjust = AlbedoSoilDir[1] * (sza_factor - 1.0);
            AlbedoSoilDir[1] *= sza_factor;
            
            // 调整顶层吸收通量
            energy->AbsShortDir[0][BAND_NIR] -= flux_sza_adjust;
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
            AlbedoSnowDfs[0] = 0.0;
            AlbedoSnowDfs[1] = 0.0;
        }
    }

    return (0);
}