/******************************************************************************
* @section DESCRIPTION
*
* Calculate soil thermal conduction.
******************************************************************************/

#include "vic_run.h"

/******************************************************************************
* @brief    Soil thermal conductivity calculated using Johansen's method.
*
* @note     Reference: Farouki, O.T., "Thermal Properties of Soils" 1986
*               Chapter 7: Methods for Calculating the Thermal Conductivity
*               of Soils
******************************************************************************/
double
soil_conductivity(double liq,
                  double ice,
                  double clay_node,
                  double sand_node,
                  double soil_pore,
                  double excess_ice,
                  double vol_sand,
                  double vol_silt,
                  double vol_clay,
                  double vol_organic,
                  double vol_gravel)
{
    // 形状因子
    const double GA_QUARTZ = 0.144;   // 沙粒形状因子
    const double GA_SILT = 0.144;     // 粉粒形状因子
    const double GA_CLAY = 0.125;     // 黏土形状因子
    const double GA_OM = 0.5;         // 有机质形状因子
    const double GA_ROCK = 0.333;     // 岩石形状因子
        
    double air = 0.0;
    if (excess_ice > 0.0) {
        air  = 0.0;
    } else {
        air = soil_pore - liq - ice;
    }
    if (air < 0.0) {
        air = 0.0;
    }
    
    // 根据黏粒和沙粒含量确定适用下限
    double VLMT = 0.10 + 0.2 * clay_node - 0.1 * sand_node;
    if (VLMT < 0.05) {
        VLMT = 0.05;
    }
    if (VLMT > 0.15) {
        VLMT = 0.15;
    }
    
    // 根据含水量选择计算方法
    double TK;
    
    if (liq > VLMT) {
        // 空气的形状因子随含水量变化
        double GAAIR = 0.035 + 0.298 * (liq - VLMT) / (soil_pore - VLMT);
        double W_air = devries_weight(CONST_KFWICE, CONST_KDAIR, GAAIR);
        
        // 湿润条件下的权重因子
        double W_liq = 1.0;
        double W_ice = 1.0;
        double W_quartz = devries_weight(CONST_KFWICE, CONST_KQUARTZ, GA_QUARTZ);
        double W_silt = devries_weight(CONST_KFWICE,  CONST_KSILT, GA_SILT);
        double W_clay = devries_weight(CONST_KFWICE, CONST_KCLAY, GA_CLAY);
        double W_om = devries_weight(CONST_KFWICE, CONST_KORGANIC, GA_OM);
        double W_rock = devries_weight(CONST_KFWICE, CONST_KGRAVEL, GA_ROCK);
        
        // 加权平均
        double numerator = W_quartz * vol_sand * CONST_KQUARTZ
                         + W_silt * vol_silt *  CONST_KSILT
                         + W_clay * vol_clay * CONST_KCLAY
                         + W_om * vol_organic * CONST_KORGANIC
                         + W_rock * vol_gravel * CONST_KGRAVEL
                         + W_liq * liq * CONST_KFWICE
                         + W_ice * ice * CONST_KICE
                         + W_air * air * CONST_KDAIR;
        
        double denominator = W_quartz * vol_sand
                           + W_silt * vol_silt
                           + W_clay * vol_clay
                           + W_om * vol_organic
                           + W_rock * vol_gravel
                           + W_liq * liq
                           + W_ice * ice
                           + W_air * air;
        
        TK = numerator / denominator;
        
    } 
    else {     
        // 计算干燥土壤热导率
        double W_air_dry = 1.0;
        double W_quartz_dry = devries_weight(CONST_KDAIR, CONST_KQUARTZ, GA_QUARTZ);
        double W_silt_dry = devries_weight(CONST_KDAIR, CONST_KSILT, GA_SILT);
        double W_clay_dry = devries_weight(CONST_KDAIR, CONST_KCLAY, GA_CLAY);
        double W_om_dry = devries_weight(CONST_KDAIR, CONST_KORGANIC, GA_OM);
        double W_rock_dry = devries_weight(CONST_KDAIR, CONST_KGRAVEL, GA_ROCK);
        double W_ice_dry = 1.0;
        
        // 干燥时的空气体积 = 总孔隙度 - 冰
        double V_air_dry = soil_pore - ice;
        if (V_air_dry < 0.0) {
            V_air_dry = 0.0;
        }
        
        double numerator_dry = W_quartz_dry * vol_sand * CONST_KQUARTZ
                             + W_silt_dry * vol_silt * CONST_KSILT
                             + W_clay_dry * vol_clay * CONST_KCLAY
                             + W_om_dry * vol_organic * CONST_KORGANIC
                             + W_rock_dry * vol_gravel * CONST_KGRAVEL
                             + W_ice_dry * ice * CONST_KICE
                             + W_air_dry * V_air_dry * CONST_KDAIR;
        
        double denominator_dry = W_quartz_dry * vol_sand
                               + W_silt_dry * vol_silt
                               + W_clay_dry * vol_clay
                               + W_om_dry * vol_organic
                               + W_rock_dry * vol_gravel
                               + W_ice_dry * ice
                               + W_air_dry * V_air_dry;
        
        double TK_dry = numerator_dry / denominator_dry;
        
        // 计算湿润下限热导率 (含水量 = VLMT)
        double V_air_wet = soil_pore - ice - VLMT;
        if (V_air_wet < 0.0) {
            V_air_wet = 0.0;
        }
        double GAAIR_wet = 0.035;  // 湿润下限时的空气形状因子
        double W_air_wet = devries_weight(CONST_KFWICE, CONST_KDAIR, GAAIR_wet);
        
        // 湿润条件下的权重因子
        double W_water_wet = 1.0;
        double W_ice_wet = 1.0;
        double W_quartz_wet = devries_weight(CONST_KFWICE, CONST_KQUARTZ, GA_QUARTZ);
        double W_silt_wet = devries_weight(CONST_KFWICE, CONST_KSILT, GA_SILT);
        double W_clay_wet = devries_weight(CONST_KFWICE, CONST_KCLAY, GA_CLAY);
        double W_om_wet = devries_weight(CONST_KFWICE, CONST_KORGANIC, GA_OM);
        double W_rock_wet = devries_weight(CONST_KFWICE, CONST_KGRAVEL, GA_ROCK);
        
        double numerator_wet = W_quartz_wet * vol_sand * CONST_KQUARTZ
                             + W_silt_wet * vol_silt * CONST_KSILT
                             + W_clay_wet * vol_clay * CONST_KCLAY
                             + W_om_wet * vol_organic * CONST_KORGANIC
                             + W_rock_wet * vol_gravel * CONST_KGRAVEL
                             + W_water_wet * VLMT * CONST_KFWICE
                             + W_ice_wet * ice * CONST_KICE
                             + W_air_wet * V_air_wet * CONST_KDAIR;
        
        double denominator_wet = W_quartz_wet * vol_sand
                               + W_silt_wet * vol_silt
                               + W_clay_wet * vol_clay
                               + W_om_wet * vol_organic
                               + W_rock_wet * vol_gravel
                               + W_water_wet * VLMT
                               + W_ice_wet * ice
                               + W_air_wet * V_air_wet;
        
        double TK_wet = numerator_wet / denominator_wet;
        
        // 线性插值
        if (VLMT > 0.0) {
            TK = TK_dry + (TK_wet - TK_dry) * liq / VLMT;
        } else {
            TK = TK_dry;
        }
    }
    
    if (TK < CONST_KDAIR) {
        TK = CONST_KDAIR;
    }
    if (TK > CONST_KQUARTZ) {
        TK = CONST_KQUARTZ;
    }
    
    return TK;
}

/******************************************************************************
* @brief    This subroutine calculates the weight of the effective thermal 
            conductivity for each component in the porous medium.
******************************************************************************/
double 
devries_weight(double lambda0, 
               double lambda1,
               double ga) 
{
    double ratio = lambda1 / lambda0 - 1.0;
    double term1 = 2.0 / 3.0 / (1.0 + ratio * ga);
    double term2 = 1.0 / 3.0 / (1.0 + ratio * (1.0 - 2.0 * ga));
    return term1 + term2;
}

/******************************************************************************
* @brief    This subroutine calculates the soil volumetric heat capacity
            based on the fractional volume of its component parts.
******************************************************************************/
double
volumetric_heat_capacity(double soil_pore,
                         double liq,
                         double ice,
                         double soil_T,
                         double matric,
                         double pressure,
                         double excess_ice,
                         double vol_sand,
                         double vol_silt,
                         double vol_clay,
                         double vol_organic,
                         double vol_gravel)
{
    const double BD_gravel_kg = 2800;   // kg/m3
    const double BD_mineral_kg = 2710;  // kg/m3
    const double BD_organic_kg = 1300;  // kg/m3
    double Cs = 0.0;
    double air = 0.0;
    double esat_T = 0.0;
    double qsdT = 0.0;
    double qsaT = 0.0;
    if (excess_ice > 0.0) {
        air = 0.0;
    } else {
        air = soil_pore - liq - ice;
    }
    if (air < 0.0) {
        air = 0.0;
    }
    // Constant values are volumetric heat capacities in J/m^3/K
    Cs = BD_gravel_kg * vol_gravel * CONST_CPGRAVEL +       // gravel
         BD_mineral_kg * (vol_sand + vol_silt + vol_clay) * CONST_CPMINE +  // mineral
         BD_organic_kg * vol_organic * CONST_CPORGANIC;     // organic 
    Cs += CONST_CPFWICE * liq * CONST_RHOFW;                // liquid water
    Cs += CONST_CPDAIR * air * CONST_RHODAIR;               // air
    Cs += CONST_CPICE * ice * CONST_RHOICE;                 // ice
    if (excess_ice > 0.0) {
        Cs += CONST_CPICE * excess_ice * CONST_RHOICE;      // excess ice
    }
    if (matric < 0.0 && air > 0.0) {
        double rel_humid = exp(CONST_MWWV * CONST_G / CONST_RGAS / soil_T * matric);
        // 潜热贡献
        svp_flags(soil_T, pressure,
                  &esat_T, &qsaT, 
                  NULL, &qsdT, 
                  ESAT | QSAT | QSDT);
        double e_actual = esat_T * rel_humid;
        double air_density = (pressure - 0.378 * e_actual) / (CONST_RDAIR * soil_T);
        double dair_dT = -air_density / soil_T;  // 理想气体近似
        double sat_vap_dens_dT = qsdT * air_density + qsaT * dair_dT;
        Cs += air * CONST_LATVAP * rel_humid * sat_vap_dens_dT;
    }

    return (Cs);
}

/******************************************************************************
* @brief    This subroutine sets the thermal node soil parameters to constant
*           values based on those defined for the current grid cells soil type.
*           Thermal node propertiers for the energy balance solution are also
*           set (these constants are used to reduce the solution time required
*           within each iteration).
******************************************************************************/
void
set_node_parameters(size_t  Nbedrock,
                    double *depth,
                    double *Zsum_soil,
                    double *array_node,
                    double *array_layer)
{
    extern option_struct options;

    size_t  nidx, lidx;
    char PAST_BOTTOM = false;
    lidx = 0;
    double Lsum = 0.;
    size_t Nlayer = options.Nlayer;

    /* set node parameters */
    for (nidx = 0; nidx < Nbedrock - 1; nidx++) {
        if (Zsum_soil[nidx] == Lsum + depth[lidx] && nidx != 0 && lidx !=
            Nlayer - 1) {
            /* node on layer boundary */
            array_node[nidx] = (array_layer[lidx] + array_layer[lidx+1]) / 2.0;
        }
        else {
            /* node completely in layer */
            array_node[nidx] = array_layer[lidx];
        }
        if (Zsum_soil[nidx] > Lsum + depth[lidx] && !PAST_BOTTOM) {
            Lsum += depth[lidx];
            lidx++;
            if (lidx == Nlayer) {
                PAST_BOTTOM = true;
                lidx = Nlayer - 1;
            }
        }
    }
}

/******************************************************************************
* @brief    This subroutine determines the moisture and ice contents of each
*           soil thermal node based on the current node temperature and layer
*           moisture content.  Thermal conductivity and volumetric heat
*           capacity are then estimated for each node based on the division of
*           moisture contents.
******************************************************************************/
int
distribute_node_moisture_properties(cell_data_struct *cell,
                                    soil_con_struct  *soil_con)
{
    extern option_struct     options;
    size_t nidx;
    size_t Nsoil = cell->Nsoil;
    double equil_liq = 0.0;
    double *moist = cell->moist;
    double *ice = cell->ice;
    double *liq = cell->liq;
    double *soil_T = cell->soil_T;
    double *Wsat_node = soil_con->Wsat_node;
    double *porosity = cell->porosity;
    double *matric = cell->matric;
    double *excess_ice = cell->excess_ice;

    /* node estimates */
    for (nidx = 0; nidx < Nsoil; nidx++) {
        if (cell->IS_GLAC) {
            ice[nidx] = 1.0;
            moist[nidx] = 1.0;
            porosity[nidx] = 0.0;
        }
        else {
            if (soil_T[nidx] < CONST_TKFRZ) {
                /* compute moisture and ice contents */
                equil_liq = frozen_soil(nidx, CONST_TKFRZ,
                                        soil_T[nidx],
                                        liq, ice,
                                        soil_con);
                liq[nidx] = equil_liq;
                ice[nidx] = (moist[nidx] - liq[nidx]) * CONST_RHOFW / CONST_RHOICE;
                if (ice[nidx] < 0) {
                    ice[nidx] = 0;
                }
            }
            else {
                liq[nidx] = moist[nidx];
                ice[nidx] = 0.0;
            }

            // 计算有效孔隙度和过量冰含量
            double total_vol = liq[nidx] + ice[nidx];
            if (total_vol > Wsat_node[nidx]) {
                excess_ice[nidx] = total_vol - Wsat_node[nidx];
                ice[nidx] = Wsat_node[nidx] - liq[nidx];
                porosity[nidx] = 0.0;
            } 
            else {
                excess_ice[nidx] = 0.0;
                porosity[nidx] = Wsat_node[nidx] - ice[nidx];
                if (porosity[nidx] < 0.0) {
                    porosity[nidx] = 0.0;
                } 
            }
            // 计算土壤基质势
            matric[nidx] = SoilWaterRetentionCurve(MATRIC_FLAG, nidx,
                                                   liq[nidx], 0.0, soil_con);
        }
    }
    return (0);
}

/******************************************************************************
* @brief    This subroutine calculates soil porocity and The volumetric 
            fractions of soil solids needed for soil parameter estimations.
******************************************************************************/
int
calc_solids_fractions(soil_con_struct *soil_con)
{
    // 定义常数
    const double V_pores_gravel = 0.24;  
    const double BD_organic = 1.3;   // g/cm3
    const double BD_mineral = 2.71;  // g/cm3
    const double BD_gravel = 2.80;   // g/cm3
    double *vol_clay = soil_con->vol_clay;
    double *vol_sand = soil_con->vol_sand;
    double *vol_silt = soil_con->vol_silt;
    double *vol_gravel = soil_con->vol_gravel;
    double *vol_organic = soil_con->vol_organic;
    double *clay_node = soil_con->clay_node; 
    double *sand_node = soil_con->sand_node;
    double *silt_node = soil_con->silt_node;
    double *soil_pore = soil_con->soil_pore;
    double *gravel_node = soil_con->gravel_node;
    double *organic_node = soil_con->organic_node;
    double *bulk_dens_avg = soil_con->bulk_dens_avg;
    double *bulk_dens_node = soil_con->bulk_dens_node;

    for (size_t i = 0; i < soil_con->Nbedrock-1; i++) {
        double wf_om_fine = 1.724 * organic_node[i];
        // double vf_om_fine = min(wf_om_fine * bulk_dens_node[i] / BD_organic, 1.0);
        vol_gravel[i] = gravel_node[i] * (1.0 - V_pores_gravel);
        bulk_dens_avg[i] = (1.0 - vol_gravel[i] / (1.0 - V_pores_gravel)) * 
                            bulk_dens_node[i] + vol_gravel[i] * BD_gravel;

        double wf_gravel = vol_gravel[i] * BD_gravel / bulk_dens_avg[i];
        double wf_sand = sand_node[i] * (1.0 - wf_om_fine) * (1.0 - wf_gravel);
        double wf_silt = silt_node[i] * (1.0 - wf_om_fine) * (1.0 - wf_gravel);
        double wf_clay = clay_node[i] * (1.0 - wf_om_fine) * (1.0 - wf_gravel);
        double wf_organ = wf_om_fine * (1.0 - wf_gravel);
        // 各组分在全土中的体积分数（固相）
        vol_sand[i] = wf_sand * bulk_dens_avg[i] / BD_mineral;
        vol_clay[i] = wf_clay * bulk_dens_avg[i] / BD_mineral;
        vol_silt[i] = wf_silt * bulk_dens_avg[i] / BD_mineral;
        vol_organic[i] = wf_organ * bulk_dens_avg[i] / BD_organic;
        double BD_particle_inv = wf_gravel / BD_gravel +
                (1.0 - wf_gravel) * ((1.0 - wf_om_fine) / BD_mineral + wf_om_fine / BD_organic);
        double BD_particle = 1.0 / BD_particle_inv;
        // 土壤孔隙度
        soil_pore[i] = 1.0 - bulk_dens_avg[i] * BD_particle_inv;
        if (soil_pore[i] <= 0.0) {
            log_err("Error: negative soil porosity. bulk density = %.4f, particle density = %.4f", 
                    bulk_dens_avg[i], BD_particle);
        }

        // 一致性检查
        double error = vol_gravel[i] + vol_organic[i] + vol_sand[i] + 
                       vol_clay[i] + vol_silt[i] - (1.0 - soil_pore[i]);
        if (fabs(error) > 1.0e-3) {
            log_err("Error in soil volumetric calculation at layer %zu: "
                    "sum of solid volume fractions minus (1 - porosity) = %.6f, ", i, error);
        }
    }
    return (0);
}
