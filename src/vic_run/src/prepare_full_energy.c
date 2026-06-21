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
    size_t          i, lidx;
    // 初始化
    size_t Nsnow = snow->Nsnow;
    size_t Nnode = cell->Nnode;
    size_t Nsoil = cell->Nsoil;
    double tmp_ice = 0.;
    double dzp, k_int;
    double tmp_density = 0.;
    // 指针赋值
    double *liq = cell->liq;
    double *ice = cell->ice;
    double *clay_node = soil_con->clay_node;
    double *sand_node = soil_con->sand_node;
    double *silt_node = soil_con->silt_node;
    double *gravel_node = soil_con->gravel_node;
    double *soil_dens_min = soil_con->soil_dens_min;
    double *moist = cell->moist;
    double *Cs_node = energy->Cs_node;
    double *pack_ice = snow->pack_ice;
    double *pack_liq = snow->pack_liq;
    double *zc_snow = snow->zc_snow;
    double *dz_snow = snow->dz_snow;
    double *soil_T = cell->soil_T;
    double *matric = cell->matric;
    double *zc_soil = soil_con->zc_soil;
    double *Zsum_soil = soil_con->Zsum_soil;
    double *Zsum_snow = snow->Zsum_snow;
    double *organic_node = soil_con->organic_node;
    double *theta_ice = snow->theta_ice;
    double *theta_liq = snow->theta_liq;
    double *porosity = snow->porosity;
    double *Wsat_node = soil_con->Wsat_node;
    double *kappa_node = energy->kappa_node;
    double *kappa_int = energy->kappa_int;
    double *bulk_dens_node = soil_con->bulk_dens_node;
    double *soil_dens_org = soil_con->soil_dens_org;
    // 计算雪的热导率和热容量
    if (Nsnow > 0) {
        for (i = 0; i < Nsnow; i++) {
            theta_ice[i] = min(1.0, pack_ice[i] / (dz_snow[i] * CONST_RHOICE));
            porosity[i] = 1.0 - theta_ice[i];
            theta_liq[i] = min(porosity[i], pack_liq[i] / (dz_snow[i] * CONST_RHOFW));

            Cs_node[i] = max(param.TOL_A, (pack_ice[i] * CONST_CPICE +
                     pack_liq[i] * CONST_CPFWICE) / dz_snow[i]);
            
            kappa_node[i] = CONST_KDAIR + (7.75e-5 * tmp_density + 1.105e-6 * 
                            tmp_density * tmp_density) * (CONST_KICE - CONST_KDAIR);
        }
    }

    // 计算表层热导率和热容量
    size_t tmp_Nsnow = Nsnow;
    if (cell->h2osfc > param.TOL_A) {
        if (cell->IS_GLAC) {  
            // 计算冰川热属性
            Cs_node[Nsnow] = cell->h2osfc * CONST_RHOICE * 
                    CONST_CPICE + liq[0] * CONST_RHOFW * CONST_CPFWICE;
            kappa_node[Nsnow] = 9.828 * exp(-0.0057 * cell->h2osfc_T);
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
    for (i = 0; i < Nsoil; i++) {
        lidx = tmp_Nsnow + i;
        tmp_ice = max(moist[i] - liq[i], 0.0);
        // 土壤节点体积热容
        Cs_node[lidx] = volumetric_heat_capacity(Wsat_node[i],
                                                    liq[i], tmp_ice, 
                                                    soil_T[i], 
                                                    moist[i], matric[i],
                                                    pressure,
                                                    organic_node[i],
                                                    bulk_dens_node[i]);
        // 土壤节点导热率
        kappa_node[lidx] = soil_conductivity(liq[i], ice[i],
                                                clay_node[i], sand_node[i], 
                                                silt_node[i], gravel_node[i],
                                                organic_node[i], 
                                                bulk_dens_node[i], 
                                                soil_dens_min[i],
                                                soil_dens_org[i]);
    }
    // 计算基岩热导率和热容量
    lidx = Nnode - 1;
    Cs_node[lidx] = CONST_CPSOIL;
    kappa_node[lidx] = CONST_KGRAVEL;

    // ===================== 计算界面热导率（调和平均） =====================
    for (i = 0; i < Nsnow; i++) {
        if (i == Nsnow - 1) {
            if (cell->h2osfc > param.TOL_A) {
                dzp = zc_snow[i] + 0.5 * cell->h2osfc;
                k_int = kappa_node[i] * kappa_node[i+1] * dzp /
                                (kappa_node[i] * 0.5 * cell->h2osfc +
                                kappa_node[i+1] * zc_snow[i]);
                // 除以节点间距，得到等效热导[W/m2/K]
                kappa_int[i] = k_int / dzp;
            }
            else {
                dzp = zc_soil[0] + zc_snow[i];
                k_int = kappa_node[i] * kappa_node[i+1] * dzp /
                                (kappa_node[i] * zc_soil[0] +
                                kappa_node[i+1] * zc_snow[i]);
                // 除以节点间距，得到等效热导[W/m2/K]
                kappa_int[i] = k_int / dzp;
            }
        }
        else {
            dzp = zc_snow[i+1] - zc_snow[i];
            // 调和平均公式
            k_int = kappa_node[i] * kappa_node[i+1] * dzp /
                            (kappa_node[i] * (zc_snow[i+1] - Zsum_snow[i]) +
                            kappa_node[i+1] * (Zsum_snow[i] - zc_snow[i]));
            // 除以节点间距，得到等效热导[W/m2/K]
            kappa_int[i] = k_int / dzp;
        }
    }
    if (cell->frac_h2o > param.TOL_A) {
        dzp = 0.5 * cell->h2osfc + zc_soil[0];
        k_int = kappa_node[Nsnow] * kappa_node[Nsnow+1] * dzp /
                        (kappa_node[Nsnow] * zc_soil[0] +
                        kappa_node[Nsnow+1] * 0.5 * cell->h2osfc);
        kappa_int[Nsnow] = k_int / dzp;
    }
    for (i = 0; i < Nsoil - 1; i++) {
        lidx = tmp_Nsnow + i;
        dzp = zc_soil[i+1] - zc_soil[i];
        // 调和平均公式
        k_int = kappa_node[lidx] * kappa_node[lidx+1] * dzp /
                        (kappa_node[lidx] * (zc_soil[i+1] - Zsum_soil[i]) +
                        kappa_node[lidx+1] * (Zsum_soil[i] - zc_soil[i]));
        // 除以节点间距，得到等效热导[W/m2/K]
        kappa_int[lidx] = k_int / dzp;
    }
    kappa_int[Nnode-1] = 0.0;
}
