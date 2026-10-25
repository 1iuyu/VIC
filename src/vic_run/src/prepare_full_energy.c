/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine returns the soil thermal properties, moisture and ice
 * contents for the top two layers for use with the QUICK_FLUX ground heat flux
 * solution.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    This subroutine returns the soil thermal properties, moisture and
 *           ice contents for the layers.
 *****************************************************************************/
void
prepare_full_energy(double             pressure,
                    cell_data_struct  *cell,
                    energy_bal_struct *energy,
                    snow_data_struct  *snow,
                    soil_con_struct   *soil_con)
{
    extern parameters_struct param;
    size_t i, lidx;
    // 初始化
    size_t Nsnow = snow->Nsnow;
    size_t Nnode = cell->Nnode;
    size_t Nsoil = cell->Nsoil;
    double dzp, zdn, zdp, denom, k_int;
    double enthliq, enthice, enthair;
    double tmp_density = 0.0;
    double qsaT = 0.0;
    double qsdT = 0.0;
    double h2osfc = cell->h2osfc;
    double coverage = snow->coverage;
    double zc_node[MAX_SNOWS+1] = {0};
    double Zsum_node[MAX_SNOWS+1] = {0};
    // 指针赋值
    double *T = energy->T;
    double *liq = energy->liq;
    double *ice = energy->ice;
    double *clay_node = soil_con->clay_node;
    double *sand_node = soil_con->sand_node;
    double *vol_clay = soil_con->vol_clay;
    double *vol_sand = soil_con->vol_sand;
    double *vol_silt = soil_con->vol_silt;
    double *vol_gravel = soil_con->vol_gravel;
    double *vol_organic = soil_con->vol_organic;
    double *Cs_node = energy->Cs_node;
    double *Cm_node = energy->Cm_node;
    double *pack_ice = snow->pack_ice;
    double *pack_liq = snow->pack_liq;
    double *Wsat_node = soil_con->Wsat_node;
    double *pack_T = snow->pack_T;
    double *zc_snow = snow->zc_snow;
    double *dz_snow = snow->dz_snow;
    double *dz_soil = soil_con->dz_soil;
    double *soil_T = cell->soil_T;
    double *matric = cell->matric;
    double *last_matric = cell->last_matric;
    double *compress = cell->compress;
    double *dcompress_dpsi = cell->dcompress_dpsi;
    double *moist = energy->moist;
    double *last_moist = energy->last_moist;
    double *dCs_dT = energy->dCs_dT;
    double *dCm_dT = energy->dCm_dT;
    double *dfliq_dT = energy->dfliq_dT;
    double *dtheta_dT = energy->dtheta_dT;
    double *dEdW_left = energy->dEdW_left;
    double *dEdW_right = energy->dEdW_right;
    double *dEdT_left = energy->dEdT_left;
    double *dEdT_right = energy->dEdT_right;
    double *dCs_dpsi = energy->dCs_dpsi;
    double *dCm_dpsi = energy->dCm_dpsi;
    double *dCs_dtheta = energy->dCs_dtheta;
    double *dtheta_dpsi = cell->dtheta_dpsi;
    double *dkappa_dT = energy->dkappa_dT;
    double *dkappa_dtheta = energy->dkappa_dtheta;
    double *excess_ice = cell->excess_ice;
    double *theta_liq = snow->theta_liq;
    double *theta_ice = snow->theta_ice;
    double *zc_soil = soil_con->zc_soil;
    double *Zsum_soil = soil_con->Zsum_soil;
    double *Zsum_snow = snow->Zsum_snow;
    double *soil_pore = soil_con->soil_pore;
    double *kappa_node = energy->kappa_node;
    double *kappa_int = energy->kappa_int;
    double *enthalpy = energy->enthalpy;
    // 构造临时节点深度矩阵
    for (i = 0; i < Nsnow; i++) {
        zc_node[i] = zc_snow[i];
        Zsum_node[i] = Zsum_snow[i];
    }
    if (h2osfc > param.TOL_A) {
        zc_node[Nsnow] = 0.5 * h2osfc + Zsum_snow[Nsnow-1];
        Zsum_node[Nsnow] = h2osfc + Zsum_snow[Nsnow-1];
    }
    else {
        zc_node[Nsnow] = zc_soil[0] + Zsum_snow[Nsnow-1];
        Zsum_node[Nsnow] = Zsum_soil[0] + Zsum_snow[Nsnow-1];
    }
    // 计算雪的热导率和热容量
    if (Nsnow > 0) {
        for (i = 0; i < Nsnow; i++) {
            double CP_snow = 92.96 + 7.37 * pack_T[i];
            double air = 1.0 - ice[i] - liq[i];
            Cs_node[i] = max(param.TOL_A, (pack_ice[i] * CP_snow +
                     pack_liq[i] * CONST_CPFWICE) / (dz_snow[i] * coverage));
            double difft = pack_T[i] - CONST_TKFRZ;
            double frac_liq = SnowFreezingCurve(SNOWLIQ_FLAG, pack_T[i], 0.0);
            double integral = 1.0 / param.SNOW_FRZ_SCALE * atan(param.SNOW_FRZ_SCALE * difft);
            Cm_node[i] = (CONST_RHOFW * CONST_CPICE - CONST_RHODAIR * CONST_CPDAIR * 
                         (CONST_RHOFW / CONST_RHOICE)) / (difft - integral) + (CONST_RHOFW * 
                          CONST_CPFWICE - CONST_RHODAIR * CONST_CPDAIR) * integral;
            dCm_dT[i] = (CONST_RHOFW * CONST_CPICE - CONST_RHODAIR * CONST_CPDAIR * 
                        (CONST_RHOFW / CONST_RHOICE)) * (1.0 - frac_liq) + (CONST_RHOFW * 
                         CONST_CPFWICE - CONST_RHODAIR * CONST_CPDAIR) * frac_liq;
            
            dCs_dtheta[i] = CONST_RHOFW * (-CONST_CPICE * (theta_liq[i] - 1.0) + CONST_RHOFW * 
                            theta_liq[i]) + CONST_RHODAIR * ((theta_liq[i] - 1.0) * (CONST_RHOFW / 
                            CONST_RHOICE) - theta_liq[i]) * CONST_CPDAIR;
            if (pack_T[i] < CONST_TKFRZ) {
                dCs_dT[i] = (CONST_RHOFW * (-CONST_CPICE + CONST_CPFWICE) + CONST_CPDAIR * 
                            (CONST_RHOFW / CONST_RHOICE - 1.0) * CONST_CPDAIR) * (moist[i] - last_moist[i]);
                dfliq_dT[i] = SnowFreezingCurve(DERILIQ_FLAG, pack_T[i], 0.0);
                dtheta_dT[i] = dfliq_dT[i] * (ice[i] * CONST_RHOICE / CONST_RHOFW + liq[i]);
                enthliq = CONST_RHOFW * CONST_CPFWICE * moist[i] * integral;
                enthice = CONST_RHOFW * CONST_CPICE * moist[i] * (difft - integral);
                enthair = CONST_RHODAIR * CONST_CPDAIR * (difft - moist[i] * ((CONST_RHOFW / CONST_RHOICE) * 
                          (difft - integral) + integral));
            }
            else {
                dCs_dT[i] = 0.0;
                dfliq_dT[i] = 0.0;
                dtheta_dT[i] = 0.0;
                enthliq = CONST_RHOFW * CONST_CPFWICE * moist[i] * difft;
                enthice = 0.0;
                enthair = CONST_RHODAIR * CONST_CPDAIR * (1.0 - moist[i]) * difft;
            }
            enthalpy[i] = enthliq + enthice + enthair;
            enthalpy[i] /= coverage;
            if (air > 0.0) {
                // 显热贡献
                Cs_node[i] += CONST_RHODAIR * CONST_CPDAIR * air;
                // 潜热贡献
                svp_flags(pack_T[i], pressure, 
                          NULL, &qsaT, 
                          NULL, &qsdT, QSAT | QSDT);
                double air_density = pressure / (CONST_RDAIR * pack_T[i]);
                double dair_dT = -air_density / pack_T[i];
                double drhodT = qsdT * air_density + qsaT * dair_dT;
                if (pack_T[i] < CONST_TKFRZ) {
                    Cs_node[i] += air * CONST_LATSUB * drhodT;
                }
                else {
                    Cs_node[i] += air * CONST_LATVAP * drhodT;
                }
            }
            tmp_density = (pack_ice[i] + pack_liq[i]) / (dz_snow[i] * coverage);
            kappa_node[i] = CONST_KDAIR + (7.75e-5 * tmp_density + 1.105e-6 * 
                            tmp_density * tmp_density) * (CONST_KICE - CONST_KDAIR);
        }
    }

    // 计算表层热导率和热容量
    size_t tmp_Nsnow = Nsnow;
    if (h2osfc > param.TOL_A) {
        if (cell->IS_GLAC) {  
            // 计算冰川热属性
            double h2osfc_ice = cell->h2osfc_ice;
            double h2osfc_liq = cell->h2osfc_liq;
            double h2osfc_T = cell->h2osfc_T;
            Cs_node[Nsnow] = h2osfc_ice * CONST_RHOICE * 
                    CONST_CPICE + h2osfc_liq * CONST_RHOFW * CONST_CPFWICE;
            enthalpy[Nsnow] = (h2osfc_ice * CONST_CPICE + h2osfc_liq * CONST_CPFWICE) * 
                              (h2osfc_T - CONST_TKFRZ) - h2osfc_ice * CONST_LATSUB;
            kappa_node[Nsnow] = 9.828 * exp(-0.0057 * h2osfc_T);
        }
        else if (cell->IS_WET) {
            Cs_node[Nsnow] = CONST_CPFWICE;
            kappa_node[Nsnow] = CONST_KFWICE;
        }
        if (Cs_node[Nsnow] < param.TOL_A || kappa_node[Nsnow] < param.TOL_A) {
            Cs_node[Nsnow] = param.TOL_A;
            kappa_node[Nsnow] = param.TOL_A;
        }
        tmp_Nsnow++;
    }
    double dliq_dtheta, dice_dtheta;
    double dTK_dice, dTK_dliq;

    for (i = 0; i < Nsoil; i++) {
        lidx = tmp_Nsnow + i;
        double dT_dpsi = 0.0;
        double soil_tkfrz = calc_tkfrz(i, &dT_dpsi, energy, cell, soil_con);

        // 计算土壤水分含量对土壤基质势的偏导数
        dtheta_dpsi[i] = SoilWaterRetentionCurve(DERIVW_FLAG, i, liq[i], 
                                                 matric[i], soil_T[i], soil_con);
        if (soil_T[i] < soil_tkfrz) {
            dliq_dtheta = 0.0;
            dice_dtheta = dtheta_dpsi[i];
            dfliq_dT[i] = 0.0;
            dtheta_dT[i] = SoilWaterRetentionCurve(DERIVT_FLAG, i, liq[i], 
                                                   matric[i], soil_T[i], soil_con);
        }
        else {
            dliq_dtheta = dtheta_dpsi[i];
            dice_dtheta = 0.0;
            dfliq_dT[i] = 0.0;
            dtheta_dT[i] = 0.0;
        }
             
        // 计算土壤压缩性
        soil_compress(liq[i], ice[i], Wsat_node[i], matric[i],
                      last_matric[i], &compress[i], &dcompress_dpsi[i]);

        // 土壤节点体积热容
        Cs_node[lidx] = volumetric_heat_capacity(soil_pore[i],
                                                 liq[i], ice[i], 
                                                 soil_T[i], matric[i],
                                                 pressure,
                                                 excess_ice[i],
                                                 vol_sand[i], 
                                                 vol_silt[i], vol_clay[i],
                                                 vol_organic[i], 
                                                 vol_gravel[i]);
        // 土壤节点导热率
        kappa_node[lidx] = soil_conductivity(liq[i], ice[i],
                                             clay_node[i], 
                                             sand_node[i], 
                                             soil_pore[i], 
                                             excess_ice[i],
                                             vol_sand[i], 
                                             vol_silt[i], vol_clay[i], 
                                             vol_organic[i], 
                                             vol_gravel[i],
                                             &dTK_dliq, &dTK_dice);
        // 计算土壤热导率的偏导数
        dkappa_dtheta[i] = dTK_dliq * dliq_dtheta + dTK_dice * dice_dtheta;
        dkappa_dT[i] = dTK_dliq * dtheta_dT[i] - dTK_dice * dtheta_dT[i];
        // 计算土壤体积热容量的偏导数
        double difft = soil_T[i] - CONST_TKFRZ;
        double air = soil_pore[i] - liq[i] - ice[i] - excess_ice[i];
        dCs_dtheta[i] = CONST_RHOFW * CONST_CPFWICE;
        if (soil_T[i] < soil_tkfrz) {
            Cm_node[i] = -CONST_RHODAIR * CONST_CPDAIR * difft + CONST_RHOICE * CONST_CPICE * 
                        (soil_T[i] - soil_tkfrz) + CONST_RHOFW * CONST_CPFWICE * difft;
            dCm_dT[lidx] = -CONST_RHODAIR * CONST_CPDAIR + CONST_RHOICE * CONST_CPICE;

            dCm_dpsi[lidx] = (-CONST_RHOICE * CONST_CPICE + CONST_RHOFW * CONST_CPFWICE) * dT_dpsi;
            dCs_dpsi[i] = (CONST_RHOICE * CONST_CPICE - CONST_RHODAIR * CONST_CPDAIR) * dtheta_dpsi[i];
            dCs_dT[i] = (-CONST_RHOICE * CONST_CPICE + CONST_RHOFW * CONST_CPFWICE) * (moist[i] - last_moist[i]);
            enthliq = CONST_RHOFW * CONST_CPFWICE * liq[i] * difft;
            enthice = CONST_RHOICE * CONST_CPICE * ice[i] * difft;
            enthair = CONST_RHODAIR * CONST_CPDAIR * air * difft;
        }
        else {
            Cm_node[lidx] = (-CONST_RHODAIR * CONST_CPDAIR + CONST_RHOFW * CONST_CPFWICE) * difft;
            dCm_dT[lidx] = -CONST_RHODAIR * CONST_CPDAIR + CONST_RHOFW * CONST_CPFWICE;
            dCm_dpsi[i] = 0.0;
            dCs_dpsi[i] = (CONST_RHOFW * CONST_CPFWICE - CONST_RHODAIR * CONST_CPDAIR) * dtheta_dpsi[i];
            dCs_dT[i] = 0.0;
            enthliq = CONST_RHOFW * CONST_CPFWICE * moist[i] * difft;
            enthice = 0.0;
            enthair = CONST_RHODAIR * CONST_CPDAIR * air * difft;
        }
        double enthsoil = (2800.0 * vol_gravel[i] * CONST_CPGRAVEL +
                           2710.0 * (vol_sand[i] + vol_silt[i] + vol_clay[i]) * CONST_CPMINE +
                           1300.0 * vol_organic[i] * CONST_CPORGANIC) * difft;
        enthalpy[i] = enthliq + enthice + enthair + enthsoil;
    }
    // 计算基岩热导率和热容量
    lidx = Nnode - 1;
    Cs_node[lidx] = CONST_CPSOIL;
    kappa_node[lidx] = CONST_KGRAVEL;

    double dTKdE_left, dTKdE_right;
    double dTKdW_left, dTKdW_right;
    // ===================== 计算界面热导率（调和平均） =====================
    for (i = 0; i < Nnode; i++) {
        if (i < Nsnow) {
            /* 雪层内部界面 */
            dzp = zc_node[i+1] - zc_node[i];
            zdn = Zsum_node[i] - zc_node[i];          /* 界面到下节点距离 */
            zdp = zc_node[i+1] - Zsum_node[i];        /* 界面上节点距离 */
            denom = kappa_node[i] * zdp + kappa_node[i+1] * zdn;
            /* 调和平均公式 */
            if (kappa_node[i] + kappa_node[i+1] > param.TOL_A) {
                k_int = kappa_node[i] * kappa_node[i+1] * dzp / denom;
                dTKdW_right = (kappa_node[i] * dzp - k_int * zdn) / denom * dkappa_dtheta[i+1];
                dTKdW_left  = (kappa_node[i+1] * dzp - k_int * zdp) / denom * dkappa_dtheta[i];
                dTKdE_right = (kappa_node[i] * dzp - k_int * zdn) / denom * dkappa_dT[i+1];
                dTKdE_left  = (kappa_node[i+1] * dzp - k_int * zdp) / denom * dkappa_dT[i];
            } 
            else {
                /* 加权平均公式 */
                k_int = (kappa_node[i] * zdn + kappa_node[i+1] * zdp) / dzp;
                dTKdW_right = zdp / dzp * dkappa_dtheta[i+1];
                dTKdW_left  = zdn / dzp * dkappa_dtheta[i];
                dTKdE_right = zdp / dzp * dkappa_dT[i+1];
                dTKdE_left  = zdn / dzp * dkappa_dT[i];
            }
            /* 除以节点间距，得到等效热导[W/m2/K] */
            kappa_int[i] = k_int / dzp;
        } 
        else if (i == Nsnow && cell->frac_h2o > param.TOL_A) {
            /* 雪-积水界面 */
            dzp = 0.5 * h2osfc + zc_soil[0];
            zdn = zc_soil[0];
            zdp = 0.5 * h2osfc;
            denom = kappa_node[i] * zdp + kappa_node[i+1] * zdn;

            if (kappa_node[i] + kappa_node[i+1] > param.TOL_A) {
                k_int = kappa_node[i] * kappa_node[i+1] * dzp / denom;
                dTKdW_right = (kappa_node[i] * dzp - k_int * zdn) / denom * dkappa_dtheta[i+1];
                dTKdW_left  = (kappa_node[i+1] * dzp - k_int * zdp) / denom * dkappa_dtheta[i];
                dTKdE_right = (kappa_node[i] * dzp - k_int * zdn) / denom * dkappa_dT[i+1];
                dTKdE_left  = (kappa_node[i+1] * dzp - k_int * zdp) / denom * dkappa_dT[i];
            } 
            else {
                k_int = (kappa_node[i] * zdn + kappa_node[i+1] * zdp) / dzp;
                dTKdW_right = zdp / dzp * dkappa_dtheta[i+1];
                dTKdW_left  = zdn / dzp * dkappa_dtheta[i];
                dTKdE_right = zdp / dzp * dkappa_dT[i+1];
                dTKdE_left  = zdn / dzp * dkappa_dT[i];
            }
            kappa_int[i] = k_int / dzp;
        } 
        else if (i < Nnode - 1) {
            /* 土壤层内部界面 */
            lidx = i - tmp_Nsnow;
            zdn = Zsum_soil[lidx] - zc_soil[lidx];
            zdp = zc_soil[lidx+1] - Zsum_soil[lidx];
            dzp = zc_soil[lidx+1] - zc_soil[lidx];
            denom = kappa_node[i] * zdp + kappa_node[i+1] * zdn;

            if (kappa_node[i] + kappa_node[i+1] > param.TOL_A) {
                k_int = kappa_node[i] * kappa_node[i+1] * dzp / denom;
                dTKdW_right = (kappa_node[i] * dzp - k_int * zdn) / denom * dkappa_dtheta[i+1];
                dTKdW_left  = (kappa_node[i+1] * dzp - k_int * zdp) / denom * dkappa_dtheta[i];
                dTKdE_right = (kappa_node[i] * dzp - k_int * zdn) / denom * dkappa_dT[i+1];
                dTKdE_left  = (kappa_node[i+1] * dzp - k_int * zdp) / denom * dkappa_dT[i];
            } 
            else {
                k_int = (kappa_node[i] * zdn + kappa_node[i+1] * zdp) / dzp;
                dTKdW_right = zdp / dzp * dkappa_dtheta[i+1];
                dTKdW_left  = zdn / dzp * dkappa_dtheta[i];
                dTKdE_right = zdp / dzp * dkappa_dT[i+1];
                dTKdE_left  = zdn / dzp * dkappa_dT[i];
            }
            /* 除以节点间距，得到等效热导[W/m2/K] */
            kappa_int[i] = k_int / dzp;
        } 
        else {
            /* 下边界：热导等于最下层节点的热导 */
            kappa_int[Nnode-1] = kappa_node[Nnode-1];
            dTKdW_right = dkappa_dtheta[Nnode-1];
            dTKdW_left  = 0.0;   /* 下边界无上方节点，置零或缺失值 */
            dTKdE_right = dkappa_dT[Nnode-1];
            dTKdE_left  = 0.0;
        }
        // 计算在层间界面处相对于上层和下层状态的通量导数
        dEdW_left[i]  = -dTKdW_left * (T[i+1] - T[i]) / dzp;
        dEdW_right[i] = -dTKdW_right * (T[i+1] - T[i]) / dzp;
        dEdT_left[i]  = -dTKdE_left * (T[i+1] - T[i]) / dzp + kappa_int[i]; 
        dEdT_right[i] = -dTKdE_right * (T[i+1] - T[i]) / dzp - kappa_int[i];
    }
}