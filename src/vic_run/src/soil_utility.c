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
soil_conductivity(double  liq,
                  double  ice,
                  double  clay_node,
                  double  sand_node,
                  double  soil_pore,
                  double  excess_ice,
                  double  vol_sand,
                  double  vol_silt,
                  double  vol_clay,
                  double  vol_organic,
                  double  vol_gravel,
                  double *dTK_dliq,
                  double *dTK_dice)
{
    // 形状因子
    const double GA_QUARTZ = 0.144;   // 沙粒形状因子
    const double GA_SILT = 0.144;     // 粉粒形状因子
    const double GA_CLAY = 0.125;     // 黏土形状因子
    const double GA_OM = 0.5;         // 有机质形状因子
    const double GA_ROCK = 0.333;     // 岩石形状因子
    
    double air_raw = soil_pore - liq - ice - excess_ice;
    double air = max(air_raw, 0.0);
    double dair_dliq = (air_raw > 0.0) ? -1.0 : 0.0;
    double dair_dice = (air_raw > 0.0) ? -1.0 : 0.0;
    
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
    (*dTK_dliq) = 0.0;
    (*dTK_dice) = 0.0;
    
    if (liq > VLMT) {
        // 空气的形状因子随含水量变化
        double GAAIR = 0.035 + 0.298 * (liq - VLMT) / (soil_pore - VLMT);
        double dGAAIR_dliq = 0.298 / (soil_pore - VLMT);
        double W_air = devries_weight(CONST_KFWICE, CONST_KDAIR, GAAIR);
        
        // 湿润条件下的权重因子
        double W_liq = 1.0;
        double W_ice = 1.0;
        double W_quartz = devries_weight(CONST_KFWICE, CONST_KQUARTZ, GA_QUARTZ);
        double W_silt = devries_weight(CONST_KFWICE,  CONST_KSILT, GA_SILT);
        double W_clay = devries_weight(CONST_KFWICE, CONST_KCLAY, GA_CLAY);
        double W_om = devries_weight(CONST_KFWICE, CONST_KORGANIC, GA_OM);
        double W_rock = devries_weight(CONST_KFWICE, CONST_KGRAVEL, GA_ROCK);
        double dW_air_dliq = devries_weight_dga(CONST_KFWICE, CONST_KDAIR, GAAIR) * dGAAIR_dliq;
        // 加权平均
        double numer = W_quartz * vol_sand * CONST_KQUARTZ +
                       W_silt * vol_silt *  CONST_KSILT +
                       W_clay * vol_clay * CONST_KCLAY +
                       W_om * vol_organic * CONST_KORGANIC +
                       W_rock * vol_gravel * CONST_KGRAVEL +
                       W_liq * liq * CONST_KFWICE +
                       W_ice * ice * CONST_KICE +
                       W_air * air * CONST_KDAIR;
        
        double denom = W_quartz * vol_sand +
                       W_silt * vol_silt +
                       W_clay * vol_clay +
                       W_om * vol_organic +
                       W_rock * vol_gravel +
                       W_liq * liq +
                       W_ice * ice +
                       W_air * air;
        
        TK = numer / denom;
        // 热导率对水含量的导数
        double dnumer_dliq = W_liq * CONST_KFWICE + W_air * dair_dliq * 
                             CONST_KDAIR + dW_air_dliq * air * CONST_KDAIR;
        double ddenom_dliq = W_liq + W_air * dair_dliq + dW_air_dliq * air;
        (*dTK_dliq) = (dnumer_dliq - TK * ddenom_dliq) / denom;
        // 热导率对冰含量的导数
        double dnumer_dice = W_ice * CONST_KICE + W_air * dair_dice * CONST_KDAIR;
        double ddenom_dice = W_ice + W_air * dair_dice;
        (*dTK_dice) = (dnumer_dice - TK * ddenom_dice) / denom;
    } 
    else {     
        // 计算干燥土壤热导率
        double W_air_dry = devries_weight(CONST_KFWICE, CONST_KDAIR, 0.035);
        double W_quartz_dry = devries_weight(CONST_KDAIR, CONST_KQUARTZ, GA_QUARTZ);
        double W_silt_dry = devries_weight(CONST_KDAIR, CONST_KSILT, GA_SILT);
        double W_clay_dry = devries_weight(CONST_KDAIR, CONST_KCLAY, GA_CLAY);
        double W_om_dry = devries_weight(CONST_KDAIR, CONST_KORGANIC, GA_OM);
        double W_rock_dry = devries_weight(CONST_KDAIR, CONST_KGRAVEL, GA_ROCK);
        double W_ice_dry = 1.0;
        
        // 干燥时的空气体积 = 总孔隙度 - 冰
        double V_air_dry_raw = soil_pore - ice - excess_ice;
        double V_air_dry = max(V_air_dry_raw, 0.0);
        double dVairdry_dice = (V_air_dry_raw > 0.0) ? -1.0 : 0.0;
        
        double numer_dry = W_quartz_dry * vol_sand * CONST_KQUARTZ +
                           W_silt_dry * vol_silt * CONST_KSILT +
                           W_clay_dry * vol_clay * CONST_KCLAY +
                           W_om_dry * vol_organic * CONST_KORGANIC +
                           W_rock_dry * vol_gravel * CONST_KGRAVEL +
                           W_ice_dry * ice * CONST_KICE +
                           W_air_dry * V_air_dry * CONST_KDAIR;
        
        double denom_dry = W_quartz_dry * vol_sand +
                           W_silt_dry * vol_silt +
                           W_clay_dry * vol_clay +
                           W_om_dry * vol_organic +
                           W_rock_dry * vol_gravel +
                           W_ice_dry * ice +
                           W_air_dry * V_air_dry;
        
        double TK_dry = numer_dry / denom_dry;
        // 干燥状态下对冰含量的导数
        double dnumer_dice_dry = W_ice_dry * CONST_KICE +
                                W_air_dry * dVairdry_dice * CONST_KDAIR;

        double ddenom_dice_dry = W_ice_dry +
                                W_air_dry * dVairdry_dice;

        double dTKdry_dice = (dnumer_dice_dry - TK_dry * ddenom_dice_dry) / denom_dry;
        
        // 湿润条件下的权重因子
        double W_water_wet = 1.0;
        double W_ice_wet = 1.0;
        double W_quartz_wet = devries_weight(CONST_KFWICE, CONST_KQUARTZ, GA_QUARTZ);
        double W_silt_wet = devries_weight(CONST_KFWICE, CONST_KSILT, GA_SILT);
        double W_clay_wet = devries_weight(CONST_KFWICE, CONST_KCLAY, GA_CLAY);
        double W_om_wet = devries_weight(CONST_KFWICE, CONST_KORGANIC, GA_OM);
        double W_rock_wet = devries_weight(CONST_KFWICE, CONST_KGRAVEL, GA_ROCK);
        
        // 在当前液态水含量下，计算湿润公式的热导率及其导数
        double GAAIR_current = 0.035 + 0.298 * (liq - VLMT) / (soil_pore - VLMT);
        double dGAAIR_dliq = 0.298 / (soil_pore - VLMT);
        double W_air_current = devries_weight(CONST_KFWICE, CONST_KDAIR, GAAIR_current);
        double dW_air_dliq = devries_weight_dga(CONST_KFWICE, CONST_KDAIR, GAAIR_current) * dGAAIR_dliq;

        double numer_wet_current = W_quartz_wet * vol_sand * CONST_KQUARTZ +
                                   W_silt_wet * vol_silt * CONST_KSILT +
                                   W_clay_wet * vol_clay * CONST_KCLAY +
                                   W_om_wet * vol_organic * CONST_KORGANIC +
                                   W_rock_wet * vol_gravel * CONST_KGRAVEL +
                                   W_water_wet * liq * CONST_KFWICE +
                                   W_ice_wet * ice * CONST_KICE +
                                   W_air_current * air * CONST_KDAIR;

        double denom_wet_current = W_quartz_wet * vol_sand +
                                   W_silt_wet * vol_silt +
                                   W_clay_wet * vol_clay +
                                   W_om_wet * vol_organic +
                                   W_rock_wet * vol_gravel +
                                   W_water_wet * liq +
                                   W_ice_wet * ice +
                                   W_air_current * air;

        double TK_wet_current = numer_wet_current / denom_wet_current;

        double dnumer_dliq = W_water_wet * CONST_KFWICE + W_air_current * dair_dliq * 
                             CONST_KDAIR + dW_air_dliq * air * CONST_KDAIR;
        double ddenom_dliq = W_water_wet + W_air_current * dair_dliq + dW_air_dliq * air;
        double dTKwet_dliq = (dnumer_dliq - TK_wet_current * ddenom_dliq) / denom_wet_current;

        double dnumer_dice = W_ice_wet * CONST_KICE + W_air_current * dair_dice * CONST_KDAIR;
        double ddenom_dice = W_ice_wet + W_air_current * dair_dice;
        double dTKwet_dice = (dnumer_dice - TK_wet_current * ddenom_dice) / denom_wet_current;

        // 平滑过渡
        double delta_liq = 0.1 * VLMT;
        double liq_start = VLMT - delta_liq;

        if (liq <= liq_start) {
            TK = TK_dry;
            (*dTK_dliq) = 0.0;
            (*dTK_dice) = dTKdry_dice;
        }
        else {
            double x = (liq - liq_start) / delta_liq;
            double S = x * x * (3.0 - 2.0 * x);
            double dS_dliq = 6.0 * x * (1.0 - x) / delta_liq;

            TK = (1.0 - S) * TK_dry + S * TK_wet_current;
            (*dTK_dliq) = dS_dliq * (TK_wet_current - TK_dry) + S * dTKwet_dliq;
            (*dTK_dice) = (1.0 - S) * dTKdry_dice + S * dTKwet_dice;
        }
    }

    // 热导率上下限截断
    if (TK < CONST_KDAIR) {
        TK = CONST_KDAIR;
        (*dTK_dliq) = 0.0;
        (*dTK_dice) = 0.0;
    }
    if (TK > CONST_KQUARTZ) {
        TK = CONST_KQUARTZ;
        (*dTK_dliq) = 0.0;
        (*dTK_dice) = 0.0;
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
* @brief    This subroutine calculates the weight of the effective thermal 
            conductivity for each component in the porous medium.
******************************************************************************/
double
devries_weight_dga(double lambda0, 
                   double lambda1,
                   double ga)
{
    double ratio = lambda1 / lambda0 - 1.0;
    double term1 = 1.0 + ratio * ga;
    double term2 = 1.0 + ratio * (1.0 - 2.0 * ga);
    return -(2.0 / 3.0) * ratio / (term1 * term1) + (2.0 / 3.0) * ratio / (term2 * term2);
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
    double esat_T = 0.0;
    double qsdT = 0.0;
    double qsaT = 0.0;
    double air = soil_pore - liq - ice - excess_ice;
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
distribute_node_moisture_properties(cell_data_struct  *cell,
                                    soil_con_struct   *soil_con)
{
    extern option_struct     options;
    size_t i, Nsoil = cell->Nsoil;
    double equil_liq = 0.0;
    double *moist = cell->moist;
    double *ice = cell->ice;
    double *liq = cell->liq;
    double *soil_T = cell->soil_T;
    double *Wsat_node = soil_con->Wsat_node;
    double *soil_porosity = cell->soil_porosity;
    double *matric = cell->matric;
    double *excess_ice = cell->excess_ice;
    double *total_matric = cell->total_matric;

    /* node estimates */
    for (i = 0; i < Nsoil; i++) {
        if (cell->IS_GLAC) {
            ice[i] = 1.0;
            moist[i] = 1.0;
            soil_porosity[i] = 0.0;
        }
        else {
            if (soil_T[i] < CONST_TKFRZ) {
                /* compute moisture and ice contents */
                equil_liq = frozen_soil(i, CONST_TKFRZ,
                                        soil_T[i], 
                                        cell,
                                        soil_con);
                liq[i] = equil_liq;
                ice[i] = (moist[i] - liq[i]) * CONST_RHOFW / CONST_RHOICE;
                if (ice[i] < 0) {
                    ice[i] = 0;
                }
            }
            else {
                liq[i] = moist[i];
                ice[i] = 0.0;
            }

            // 计算有效孔隙度和过量冰含量
            double total_vol = liq[i] + ice[i];
            if (total_vol > Wsat_node[i]) {
                excess_ice[i] = total_vol - Wsat_node[i];
                ice[i] = Wsat_node[i] - liq[i];
                soil_porosity[i] = 0.0;
            } 
            else {
                excess_ice[i] = 0.0;
                soil_porosity[i] = Wsat_node[i] - ice[i];
                if (soil_porosity[i] < 0.0) {
                    soil_porosity[i] = 0.0;
                } 
            }
            // 计算土壤基质势
            matric[i] = SoilWaterRetentionCurve(MATRIC_FLAG, i, liq[i], 
                                                0.0, 0.0, soil_con);
            total_matric[i] = SoilWaterRetentionCurve(MATRIC_FLAG, i, moist[i], 
                                                      0.0, 0.0, soil_con);
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

/******************************************************************************
* @brief    This subroutine compute the critical temperature above which 
*           all water is unfrozen.
******************************************************************************/
double
calc_tkfrz(size_t             nidx,
           double            *dT_dpsi,
           cell_data_struct  *cell,
           soil_con_struct   *soil_con)
{
    double tmp_tkfrz = 0.0;
    double *moist = cell->moist;
    double Wsat_node = soil_con->Wsat_node[nidx];
    if (moist[nidx] > Wsat_node) {
        moist[nidx] = Wsat_node;
    }
    double total_matric = SoilWaterRetentionCurve(MATRIC_FLAG, nidx, 
                                                  moist[nidx], 0.0, 0.0, soil_con);
    if (total_matric < 0.0) {
        double A = CONST_LATICE / CONST_G;
        tmp_tkfrz = CONST_TKTRIP * total_matric / (A - total_matric) + CONST_TKTRIP;
        if (dT_dpsi != NULL) {
            double denom = A - total_matric;
            (*dT_dpsi) = CONST_TKTRIP * A / (denom * denom);
        }
    } 
    else {
        tmp_tkfrz = CONST_TKTRIP;
        if (dT_dpsi != NULL) {
            *dT_dpsi = 0.0;
        }
    }

    return (tmp_tkfrz);
}

/******************************************************************************
* @brief    This subroutine compute the soil compressibility (-) and its 
*           derivative with respect to matric head (m-1).           
******************************************************************************/
int
soil_compress(double  liq,
              double  ice,
              double  Wsat_node,
              double  matric,
              double  last_matric,
              double *compress,
              double *dcompress_dpsi)
{
    extern parameters_struct param;
    
    (*dcompress_dpsi) = param.SOIL_SPECSTORE * (liq + ice) / Wsat_node;
    (*compress) = (matric - last_matric) * (*dcompress_dpsi);

    return (0);
}