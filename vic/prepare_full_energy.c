/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine returns the soil thermal properties, moisture and ice
 * contents for the top two layers for use with the QUICK_FLUX ground heat flux
 * solution.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    This subroutine returns the soil thermal properties, moisture and
 *           ice contents for the layers.
 *****************************************************************************/
void
prepare_full_energy(unsigned short     veg_class,
                    double             step_dt,
                    cell_data_struct  *cell,
                    energy_bal_struct *energy,
                    snow_data_struct  *snow,
                    soil_con_struct   *soil_con)
{
    extern option_struct options;
    extern parameters_struct param;

    size_t          i;
    // 初始化
    size_t Nsnow = snow->Nsnow;
    size_t Nnode = options.Nnode;
    size_t Total_Layer = Nsnow + Nnode;
    double tmp_ice = 0.;
    double tmp_density = 0.;
    double snow_depth = snow->snow_depth;
    double coverage = snow->coverage;
    // 指针赋值
    double *liq = cell->liq;
    double *moist = cell->moist;
    double *dz_snow = snow->dz_snow;
    double *Cs_node = energy->Cs_node;
    double *pack_ice = snow->pack_ice;
    double *pack_liq = snow->pack_liq;
    double *zc_node = cell->zc_node;
    double *dz_node = cell->dz_node;
    double *Zsum_node = cell->Zsum_node;
    double *quartz_node = soil_con->quartz_node;
    double *organic_node = soil_con->organic_node;
    double *Zsum_soil = soil_con->Zsum_soil;
    double *theta_ice = snow->theta_ice;
    double *theta_liq = snow->theta_liq;
    double *zc_snow = snow->zc_snow;
    double *zc_soil = soil_con->zc_soil;
    double *porosity = snow->porosity;
    double *Wsat_node = soil_con->Wsat_node;
    double *dz_soil = soil_con->dz_soil;
    double *kappa_node = energy->kappa_node;
    double *kappa_int = energy->kappa_int;
    double *fusion_fact = energy->fusion_fact;

    // 计算雪的热导率和热容量
    if (Nsnow > 0) {
        for (i = 0; i < Nsnow; i++) {
            theta_ice[i] = min(1.0, pack_ice[i] / (dz_snow[i] * CONST_RHOICE));
            porosity[i] = 1.0 - theta_ice[i];
            theta_liq[i] = min(porosity[i], pack_liq[i] / (dz_snow[i] * CONST_RHOFW));
            tmp_density = (pack_ice[i] + pack_liq[i]) / dz_snow[i];

            Cs_node[i] = max(param.TOL_A, (CONST_CPICE * theta_ice[i] + 
                                    CONST_CPFWICE * theta_liq[i]));
            
            kappa_node[i] = 2.0e-2 + 2.5e-6 * tmp_density * tmp_density;
        }
    }


    if (veg_class != options.GLACIER_ID) {  
        // 计算土壤热属性
        for (i = 0; i < Nnode; i++) {
            tmp_ice = max(moist[i] - liq[i], 0.0);
            // 土壤节点体积热容
            Cs_node[Nsnow + i] = volumetric_heat_capacity(Wsat_node[i], 
                                                  liq[i], tmp_ice, 
                                                  moist[i], organic_node[i]);
            // 土壤节点导热率
            kappa_node[Nsnow + i] = soil_conductivity(moist[i], liq[i],
                                              Wsat_node[i],
                                              quartz_node[i], organic_node[i]);
        }
    }
    else {
        // 计算冰川热属性
        for (i = 0; i < Nnode; i++) {
            Cs_node[Nsnow + i] = 1.0e6 * (0.8194 + 0.1309 * dz_soil[i]);
            kappa_node[Nsnow + i] = 0.32333 + (0.10073 * dz_soil[i]);
        }
    }

    // ===================== 新增部分：计算界面热导率（调和平均） =====================

    // 计算节点深度和累计深度
    if (Nsnow > 0) {
        for (i = 0; i < Nsnow; i++) {
            zc_node[i] = zc_snow[i];
            dz_node[i] = dz_snow[i];
            Zsum_node[i] = Zsum_soil[i];  
        }
    }
    for (i = 0; i < Nnode; i++) {
        zc_node[Nsnow + i] = zc_soil[i];
        dz_node[Nsnow + i] = dz_soil[i];
        Zsum_node[Nsnow + i] = Zsum_soil[i];
    }

    for (i = 0; i < Total_Layer - 1; i++) { 
        // 调和平均公式
        kappa_int[i] = kappa_node[i] * kappa_node[i + 1] * (zc_node[i + 1] - zc_node[i]) /
                        (kappa_node[i] * (zc_node[i + 1] - Zsum_node[i]) +
                         kappa_node[i + 1] * (Zsum_node[i] - zc_node[i]));
    }
    kappa_int[Total_Layer-1] = 0.0;
    
    // 计算用于融化和冻结的临时变量
    if (Nsnow > 0) {
        for (i = 0; i < Nsnow; i++) {
            fusion_fact[i] = step_dt / (Cs_node[i] + dz_snow[i]);
        }         
    }
    for (i = 0; i < Total_Layer; i++) {
        fusion_fact[Nsnow + i] = step_dt / (Cs_node[Nsnow + i] + dz_node[i]);
    }

    // 雪层/土壤界面处理
    if (Nsnow == 0) {
        kappa_node[0] = (kappa_node[0] * dz_node[0] +
                        0.35 * snow_depth) / (snow_depth + dz_node[0]);
    }
    else {
        kappa_node[Nsnow] = (kappa_node[Nsnow] * dz_node[Nsnow] + 
                        kappa_node[Nsnow - 1] * dz_snow[Nsnow - 1]) / (dz_snow[Nsnow - 1] + dz_node[Nsnow]);
    }

}
