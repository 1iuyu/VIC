/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate the aerodynamic resistances.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    Calculate the aerodynamic resistance for each vegetation layer,
 *           based on Monin-Obukhov (M-O) Similarity Theory (MOST).
 *****************************************************************************/
void
ActiveLayer(cell_data_struct *cell,
            soil_con_struct  *soil_con)
{
    size_t Nsoil = cell->Nsoil;
    int nidx;
    size_t Nthaw = 0;
    size_t Nfrost = 0;
    double *Zsum_soil = soil_con->Zsum_soil;
    double *ice = cell->ice;
    double *liq = cell->liq;
    double *moist = cell->moist;

    bool found_frost = false;

    /*************************************************
     * 从底向上扫描：找冻结层和融化层
     *************************************************/
    for (nidx = Nsoil - 1; nidx >= 0; nidx--) {

        if (ice[nidx] > 0.0) {
            if (!found_frost) {
                Nfrost = nidx;
                found_frost = true;
            }
        }
        else {
            if (found_frost && Nthaw == 0) {
                Nthaw = nidx + 1;
                break;
            }
        }
    }

    /*************************************************
     * 计算冻结深度 FDEPTH
     *************************************************/
    double fdepth = 0.0;

    if (Nfrost > 0 || (Nfrost == 0 && ice[0] > 0.0)) {

        size_t nf = Nfrost;

        double fractn = max(0.0, ice[nf] / moist[nf]);

        if (nf == 0) {
            fdepth = fractn * (Zsum_soil[1] - Zsum_soil[0]) / 2.0;
        }
        else if (nf == Nsoil - 1) {
            fdepth = (Zsum_soil[nf] + Zsum_soil[nf - 1]) / 2.0
                + fractn * (Zsum_soil[nf] - Zsum_soil[nf - 1]) / 2.0;
        }
        else {
            fdepth = (Zsum_soil[nf] + Zsum_soil[nf - 1]) / 2.0
                + fractn * (Zsum_soil[nf + 1] - Zsum_soil[nf - 1]) / 2.0;
        }
    }

    /*************************************************
     * 计算融化深度 TDEPTH
     *************************************************/
    double tdepth = 0.0;

    if (Nthaw > 0) {

        size_t nt = Nthaw;

        double fractn = max(0.0, liq[nt] / moist[nt]);;

        if (nt == 0) {
            tdepth = fractn * (Zsum_soil[1] - Zsum_soil[0]) / 2.0;
        }
        else {

            if (nt == Nfrost) {
                double mid = (Zsum_soil[nt] + Zsum_soil[nt - 1]) / 2.0;
                tdepth = mid + fractn * (fdepth - mid);
            }
            else {
                tdepth = (Zsum_soil[nt] + Zsum_soil[nt - 1]) / 2.0
                    + fractn * (Zsum_soil[nt + 1] - Zsum_soil[nt - 1]) / 2.0;
            }
        }
    }
    // 写回结构体
    cell->fdepth = fdepth;
    cell->tdepth = tdepth;
    cell->Nthaw  = Nthaw;
    cell->Nfrost = Nfrost;

}

/******************************************************************************
 * @brief    Calculate the root moisture stress based on the active layer depth and soil moisture conditions.
 *****************************************************************************/
void
calc_root_moist_stress(cell_data_struct *cell,
                       soil_con_struct  *soil_con,
                       veg_lib_struct   *veg_lib)
{
    extern option_struct options;
    size_t i, Nroot;
    double tmp_matric = 0.0;
    double total_unfrozen = 0.0;
    double f_transp = 0.0;
    double smpsc = veg_lib->smpsc; // 叶片气孔关闭时的土壤水势[m]
    double smpso = veg_lib->smpso; // 叶片气孔打开时的土壤水势[m]
    double *soil_T = cell->soil_T;
    double *Ra_root = cell->Ra_root;
    double *porosity = cell->porosity;
    double *liq = cell->liq;
    double *root = cell->root;
    double *matric = cell->matric;
    double *Wsat_node = soil_con->Wsat_node;
    double *Netroot = cell->Netroot;
    double root_unfrozen[MAX_SOILS];
    // 初始化变量
    Nroot = cell->Nroot;
    for (i = 0; i < MAX_SOILS; i++) {
        root_unfrozen[i] = 0.0;
    }
    // 使用冻土活动层计算非冻结的植被根系分数
    if (options.ACTIVE_LAYER) {
        for (i = 0; i < Nroot; i++) {
            if (i < cell->Nthaw) {
                root_unfrozen[i] = root[i];
            }
            else {
                root_unfrozen[i] = 0.0;
            }
            total_unfrozen += root_unfrozen[i];
        }
    }
    else {         
        // 不使用冻土活动层，直接使用土壤温度计算非冻结的植被根系分数                
        for (i = 0; i < Nroot; i++) {
            if (soil_T[i] >= CONST_TKFRZ) {
                root_unfrozen[i] = root[i];
            }
            else {
                root_unfrozen[i] = 0.0;
            }
            total_unfrozen += root_unfrozen[i];
        }
    }
    // 将非冻结的植被根系分数归一化
    if (total_unfrozen > 0.0) {
        for (i = 0; i < Nroot; i++) {
            root_unfrozen[i] = root_unfrozen[i] / total_unfrozen;
        }
    }
    // 计算根部水分胁迫程度
    for (i = 0; i < Nroot; i++) {
        if (liq[i] <= 0.0 || soil_T[i] <= CONST_TKFRZ - 2.0) {
            Netroot[i] = 0.0;
        }
        else {
            tmp_matric = max(matric[i], smpsc);
            Ra_root[i] = min((porosity[i] / Wsat_node[i]) * 
                        (tmp_matric - smpsc) / (smpso - smpsc), 1.0);
            Netroot[i] = root_unfrozen[i] * Ra_root[i];
            f_transp += max(Netroot[i], 0.0);
        }
    }
    
    for (i = 0; i < Nroot; i++) {
        if (f_transp > 0.0) {
            Netroot[i] /= f_transp;
        }
        else {
            Netroot[i] = 0.0;
        }
    }
    
}