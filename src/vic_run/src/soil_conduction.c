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
                  double silt_node,
                  double gravel_node,
                  double organic_node,
                  double bulk_dens_node,
                  double soil_dens_min,
                  double soil_dens_org)
{
    double mineral_frac = 0.0;
    double sand_frac = 0.0;
    double silt_frac = 0.0;
    double clay_frac = 0.0;
    
    // 形状因子
    const double GA_QUARTZ = 0.144;   // 沙粒形状因子
    const double GA_SILT = 0.144;     // 粉粒形状因子
    const double GA_CLAY = 0.125;     // 黏土形状因子
    const double GA_OM = 0.5;         // 有机质形状因子
    const double GA_ROCK = 0.333;     // 岩石形状因子
    
    // 扣除有机质和岩石后的矿物质量分数
    mineral_frac = 1.0 - organic_node - gravel_node;
    sand_frac = mineral_frac * sand_node;
    silt_frac = mineral_frac * silt_node;
    clay_frac = mineral_frac * clay_node;
    // 体积分数 = 质量分数 × 容重 / 组分密度
    double V_sand = sand_frac * bulk_dens_node / soil_dens_min;
    double V_silt = silt_frac * bulk_dens_node / soil_dens_min;
    double V_clay = clay_frac * bulk_dens_node / soil_dens_min;
    double V_organic = organic_node * bulk_dens_node / soil_dens_org;
    double V_gravel = gravel_node * bulk_dens_node / soil_dens_min;
    double porosity = 1.0 - V_sand - V_silt - V_clay - V_organic - V_gravel;
    
    double V_air = porosity - liq - ice;
    if (V_air < 0.0) {
        V_air = 0.0;
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
        double GAAIR = 0.035 + 0.298 * (liq - VLMT) / (porosity - VLMT);
        double W_air = devries_weight(CONST_KFWICE, CONST_KDAIR, GAAIR);
        
        // 湿润条件下的权重因子
        double W_water = 1.0;
        double W_ice = 1.0;
        double W_quartz = devries_weight(CONST_KFWICE, CONST_KQUARTZ, GA_QUARTZ);
        double W_silt = devries_weight(CONST_KFWICE,  CONST_KSILT, GA_SILT);
        double W_clay = devries_weight(CONST_KFWICE, CONST_KCLAY, GA_CLAY);
        double W_om = devries_weight(CONST_KFWICE, CONST_KORGANIC, GA_OM);
        double W_rock = devries_weight(CONST_KFWICE, CONST_KGRAVEL, GA_ROCK);
        
        // 加权平均：分子 = Σ (权重 × 体积 × 热导率)
        double numerator = W_quartz * V_sand * CONST_KQUARTZ
                         + W_silt * V_silt *  CONST_KSILT
                         + W_clay * V_clay * CONST_KCLAY
                         + W_om * V_organic * CONST_KORGANIC
                         + W_rock * V_gravel * CONST_KGRAVEL
                         + W_water * liq * CONST_KFWICE
                         + W_ice * ice * CONST_KICE
                         + W_air * V_air * CONST_KDAIR;
        
        // 分母 = Σ (权重 × 体积)
        double denominator = W_quartz * V_sand
                           + W_silt * V_silt
                           + W_clay * V_clay
                           + W_om * V_organic
                           + W_rock * V_gravel
                           + W_water * liq
                           + W_ice * ice
                           + W_air * V_air;
        
        TK = numerator / denominator;
        
    } 
    else {     
        // 5.1 计算干燥土壤热导率 (使用空气为连续相)
        double W_air_dry = 1.0;
        double W_quartz_dry = devries_weight(CONST_KDAIR, CONST_KQUARTZ, GA_QUARTZ);
        double W_silt_dry = devries_weight(CONST_KDAIR, CONST_KSILT, GA_SILT);
        double W_clay_dry = devries_weight(CONST_KDAIR, CONST_KCLAY, GA_CLAY);
        double W_om_dry = devries_weight(CONST_KDAIR, CONST_KORGANIC, GA_OM);
        double W_rock_dry = devries_weight(CONST_KDAIR, CONST_KGRAVEL, GA_ROCK);
        double W_ice_dry = 1.0;
        
        // 干燥时的空气体积 = 总孔隙度 - 冰 (假设无水)
        double V_air_dry = porosity - ice;
        if (V_air_dry < 0.0) {
            V_air_dry = 0.0;
        }
        
        double numerator_dry = W_quartz_dry * V_sand * CONST_KQUARTZ
                             + W_silt_dry * V_silt * CONST_KSILT
                             + W_clay_dry * V_clay * CONST_KCLAY
                             + W_om_dry * V_organic * CONST_KORGANIC
                             + W_rock_dry * V_gravel * CONST_KGRAVEL
                             + W_ice_dry * ice * CONST_KICE
                             + W_air_dry * V_air_dry * CONST_KDAIR;
        
        double denominator_dry = W_quartz_dry * V_sand
                               + W_silt_dry * V_silt
                               + W_clay_dry * V_clay
                               + W_om_dry * V_organic
                               + W_rock_dry * V_gravel
                               + W_ice_dry * ice
                               + W_air_dry * V_air_dry;
        
        double TK_dry = numerator_dry / denominator_dry;
        
        // 5.2 计算湿润下限热导率 (含水量 = VLMT)
        double V_air_wet = porosity - ice - VLMT;
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
        
        double numerator_wet = W_quartz_wet * V_sand * CONST_KQUARTZ
                             + W_silt_wet * V_silt * CONST_KSILT
                             + W_clay_wet * V_clay * CONST_KCLAY
                             + W_om_wet * V_organic * CONST_KORGANIC
                             + W_rock_wet * V_gravel * CONST_KGRAVEL
                             + W_water_wet * VLMT * CONST_KFWICE
                             + W_ice_wet * ice * CONST_KICE
                             + W_air_wet * V_air_wet * CONST_KDAIR;
        
        double denominator_wet = W_quartz_wet * V_sand
                               + W_silt_wet * V_silt
                               + W_clay_wet * V_clay
                               + W_om_wet * V_organic
                               + W_rock_wet * V_gravel
                               + W_water_wet * VLMT
                               + W_ice_wet * ice
                               + W_air_wet * V_air_wet;
        
        double TK_wet = numerator_wet / denominator_wet;
        
        // 5.3 线性插值
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
volumetric_heat_capacity(double Wsat_node,
                         double liq,
                         double ice,
                         double soil_T,
                         double moist,
                         double matric,
                         double pressure,
                         double organic_node,
                         double bulk_dens_node)
{
    double Cs = 0.0;
    double esat_T = 0.0;
    double qsdT = 0.0;
    double mineral = 1.0 - organic_node;
    double air = Wsat_node - ice - liq;
    if (air < 0.0) {
        air = 0.0;
    }
    // Constant values are volumetric heat capacities in J/m^3/K
    Cs =  CONST_CPFWICE * liq * CONST_RHOFW; // liquid water
    Cs += CONST_CPMINE * mineral * bulk_dens_node; // soil
    Cs += CONST_CPDAIR * air * CONST_RHODAIR; // air
    Cs += CONST_CPICE * ice * CONST_RHOICE; // ice
    Cs += CONST_CPORGANIC * organic_node * bulk_dens_node; // organic
    if (matric < 0.0 && air > 0.0 && moist < Wsat_node) {
        double rel_humid = exp(CONST_MWWV * CONST_G / CONST_RGAS / soil_T * matric);
        // 潜热贡献
        svp_flags(soil_T, pressure,
                  &esat_T, NULL, 
                  NULL, &qsdT, 
                  ESAT | QSDT);
        double e_actual = esat_T * rel_humid;
        double air_density = (pressure - 0.378 * e_actual) / (CONST_RDAIR * soil_T);
        Cs += air * CONST_LATVAP * rel_humid * qsdT * air_density;
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
    double *moist = cell->moist;
    double *ice = cell->ice;
    double *liq = cell->liq;
    double *soil_T = cell->soil_T;
    double *Wsat_node = soil_con->Wsat_node;
    double *porosity = cell->porosity;
    double *matric = cell->matric;

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
                double equil_liq = frozen_soil(nidx, CONST_TKFRZ,
                                               soil_T,
                                               liq, ice,
                                               soil_con);
                liq[nidx] = equil_liq;
                ice[nidx] = (moist[nidx] - liq[nidx]) * CONST_RHOFW / CONST_RHOICE;
                if (ice[nidx] < 0) {
                    ice[nidx] = 0;
                }
                matric[nidx] = SoilWaterRetentionCurve(MATRIC_FLAG, nidx,
                                                       equil_liq, 0.0, soil_con);
            }
            porosity[nidx] = Wsat_node[nidx] - ice[nidx];
        }
    }
    return (0);
}
