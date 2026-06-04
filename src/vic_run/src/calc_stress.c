/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate the transpiration stress using a plant hydraulics approach.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief  Calculate the transpiration stress using a plant hydraulics approach.
 *****************************************************************************/
int
calc_stress(double           *bsun,
            double           *bsha,
            double           *mat_VEG,
            double            thm,
            double            gb_mol,
            double            qsat_T,
            double            Qair_over,
            double            pressure,
            double            air_density,
            double            gs_mol_sun,
            double            gs_mol_sha,
            cell_data_struct *cell,
            soil_con_struct  *soil_con,
            veg_var_struct   *veg_var,
            veg_lib_struct   *veg_lib)
{
    extern parameters_struct param;
    // 初始化输出变量
    size_t i, j;
    size_t iter;
    bool flag = false;
    bool night = false;
    (*bsun) = 1.0;
    (*bsha) = 1.0;
    double qflx_sun = 0.0; // Sunlit leaf transpiration [kg/m2/s]
    double qflx_sha = 0.0; // Shaded leaf transpiration [kg/m2/s]
    double dx[4]; // 水势更新量 [m]
    double mat_A[16]; // 雅可比矩阵
    double mat_RHS[4]; // 右侧项 [m/s]
    double matric50 = veg_lib->matric50; // 50%失水点的水势[m]
    double leaf_sun = veg_var->leaf_sun;
    double leaf_sha = veg_var->leaf_sha;
    for (i = 0; i < 4; i++) {
        dx[i] = 0.0;
        mat_RHS[i] = 0.0;
    }
    for (i = 0; i < 16; i++) {
        mat_A[i] = 0.0;
    }
    // 夜间处理：如果阳叶水势为正，说明植物处于夜间状态，使用阴叶水势作为阳叶水势的近似值
    if (mat_VEG[0] > 0.0) {
        night = true;
        mat_VEG[0] = mat_VEG[1];
    }
    else {
        night = false;
    }
    // 复制以避免重写原始导度值
    double gs0sun = gs_mol_sun;
    double gs0sha = gs_mol_sha;
    
    /* 计算蒸腾需求（未受胁迫的潜在蒸腾） */
    bool CALC_TRANSP = true;
    get_qflx(CALC_TRANSP, gb_mol,
             qsat_T, Qair_over, 
             pressure, air_density, 
             thm, &gs0sun, 
             &gs0sha, &qflx_sun,
             &qflx_sha, veg_var);
    
    // 检查是否存在足够的叶面积和正的蒸腾需求
    if ((leaf_sun > MIN_TOL_LAI || leaf_sha > MIN_TOL_LAI) &&
        (qflx_sun > 0.0 || qflx_sha > 0.0)) {

        iter = 0;
        flag = false;
        
        do {

            iter++;
            // 计算通量残差
            spacF(mat_VEG, mat_A, mat_RHS, 
                  qflx_sun, qflx_sha, cell, 
                  soil_con, veg_var, veg_lib);
            double sum = 0.0;
            for (i = 0; i < 4; i++) {
                sum += mat_RHS[i] * mat_RHS[i];
            }
            double f_norm = sqrt(sum);
            
            // 检查通量是否平衡
            if (f_norm < param.TOL_A * (qflx_sun + qflx_sha)) {
                flag = false;
                break;
            }
            
            // 检查是否超过最大迭代次数
            if (iter > param.MAX_ITER_MOST) {
                flag = false;
                break;
            }
            
            // 根据叶面积情况计算dx
            if (leaf_sun > MIN_TOL_LAI && leaf_sha > MIN_TOL_LAI) {
                int ipiv[4];
                int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, 
                                         4, 1, mat_A, 4, ipiv, mat_RHS, 4);

                if (info != 0) {
                    flag = true;
                    break;
                }
                for (j = 0; j < 4; j++) {
                    dx[j] = mat_RHS[j];
                }
            }
            else {              
                // 根据哪个叶面积为零来调整计算
                int ipiv[3];
                double A_3x3[9], B_3x3[3];
                if (leaf_sha > MIN_TOL_LAI) {
                    // 阳叶LAI=0，求解阴叶-木质部-根系统
                    for (i = 0; i < 3; i++) {
                        for (j = 0; j < 3; j++) {
                            A_3x3[i + j * 3] = mat_A[(1 + i) + (1 + j) * 4];
                        }
                        B_3x3[i] = mat_RHS[1 + i];
                    }
                    
                    // 使用LAPACK求解线性系统 A_3x3 * x = B_3x3
                    int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, 3, 1, A_3x3, 3, ipiv, B_3x3, 3);
                    // 检查求解是否成功
                    if (info != 0) {
                        flag = true;
                        break;
                    }                    
                    // 将解映射回4x1的dx
                    dx[0] = 0.0;
                    dx[1] = B_3x3[0];
                    dx[2] = B_3x3[1];
                    dx[3] = B_3x3[2];
                } 
                else {
                    // 构建简化的3x3系统: [sun, xyl, root]
                    int idx_map[3] = {0, 2, 3};
                    for (j = 0; j < 3; j++) {
                        for (i = 0; i < 3; i++) {
                            A_3x3[i + j * 3] = mat_A[idx_map[i] + idx_map[j] * 4];
                        }
                        B_3x3[j] = mat_RHS[idx_map[j]];
                    }

                    // 使用LAPACK求解线性系统 A_3x3 * x = B_3x3
                    int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, 3, 1, A_3x3, 3, ipiv, B_3x3, 3);
                    // 检查求解是否成功
                    if (info != 0) {
                        flag = true;
                        break;
                    }
                    dx[0] = B_3x3[0];
                    dx[1] = 0.0;
                    dx[2] = B_3x3[1];
                    dx[3] = B_3x3[2];
                }
            }
            
            // 限制步长
            double max_dx = fabs(dx[0]);
            for (j = 1; j < 4; j++) {
                if (fabs(dx[j]) > max_dx) {
                    max_dx = fabs(dx[j]);
                }
            }
            
            if (max_dx > 50.0) {
                double scale = 50.0 / max_dx;
                for (j = 0; j < 4; j++) {
                    dx[j] *= scale;
                }
            }
            
            // 更新水势
            if (leaf_sun > MIN_TOL_LAI && leaf_sha > MIN_TOL_LAI) {
                for (j = 0; j < 4; j++) {
                    mat_VEG[j] += dx[j];
                }
            } 
            else if (leaf_sha > MIN_TOL_LAI) {
                for (j = 0; j < 4; j++) {
                    mat_VEG[j] += dx[j];
                }
                mat_VEG[0] = mat_VEG[2];
            } 
            else {
                /* laisun > 0, laisha == 0 */
                mat_VEG[0] += dx[1];
                mat_VEG[1] = mat_VEG[2];
                mat_VEG[2] += dx[2];
                mat_VEG[3] += dx[3];
            }
            
            /* 检查水势变化是否足够小 */
            sum = 0.0;
            for (j = 0; j < 4; j++) {
                sum += dx[j] * dx[j];
            }
            double dx_norm = sqrt(sum);
            if (dx_norm < 1.0e-12) {
                break;
            }
            
            // 强制保证SPAC梯度朝向大气
            if (mat_VEG[2] > mat_VEG[3]) {
                mat_VEG[2] = mat_VEG[3];
            }
            if (mat_VEG[0] > mat_VEG[2]) {
                mat_VEG[0] = mat_VEG[2];
            }
            if (mat_VEG[1] > mat_VEG[2]) {
                mat_VEG[1] = mat_VEG[2];
            }
        } while (1);
        
    } else {
        flag = true; // 蒸腾通量为零
    }
    
    // 根据求解结果计算胁迫系数
    double qsun, qsha;
    double transp = 0.0;
    if (flag) {
        // 代数求解
        getvegwp(mat_VEG, gb_mol, &gs0sun, &gs0sha, 
                 qsat_T, Qair_over, &transp, 
                 pressure, air_density, thm, 
                 cell, soil_con, veg_var, veg_lib);

        *bsun = plc(mat_VEG[0], matric50);
        *bsha = plc(mat_VEG[1], matric50);
    } 
    else {
        // 计算衰减后的实际通量
        qsun = qflx_sun * plc(mat_VEG[0], matric50);
        qsha = qflx_sha * plc(mat_VEG[1], matric50);
        
        // 反推受胁迫的气孔导度
        CALC_TRANSP = false;
        get_qflx(CALC_TRANSP, gb_mol, 
                 qsat_T, Qair_over, 
                 pressure, air_density, 
                 thm, &gs0sun, 
                 &gs0sha, &qsun,
                 &qsha, veg_var);
        
        if (qflx_sun > 0.0) {
            *bsun = gs0sun / gs_mol_sun;
        } 
        else {
            *bsun = plc(mat_VEG[0], matric50);
        }
        
        if (qflx_sha > 0.0) {
            *bsha = gs0sha / gs_mol_sha;
        } 
        else {
            *bsha = plc(mat_VEG[1], matric50);
        }
    }
    
    // 对极小的胁迫系数进行截断
    if (*bsun < 0.01) {
        *bsun = 0.0;
    }
    if (*bsha < 0.01) {
        *bsha = 0.0;
    }
    
    // 夜间处理
    if (night) {
        gs0sun = (*bsun) * gs_mol_sun;
        gs0sha = (*bsha) * gs_mol_sha;
        getvegwp(mat_VEG, gb_mol, &gs0sun, &gs0sha, 
                 qsat_T, Qair_over, &transp,
                 pressure, air_density, thm, 
                 cell, soil_con, veg_var, veg_lib);

        if (transp < 0.0) {
            transp = 0.0;
        }
        cell->transp = transp;
    }

    return 0;
}

/******************************************************************************
 * @brief    Calculate temperature scaling factor.
 *****************************************************************************/
void get_qflx(bool            CALC_TRANSP,
              double          gb_mol,
              double          qsat_T, 
              double          Qair_over,
              double          pressure, 
              double          air_density, 
              double          thm,
              double         *gs_mol_sun, 
              double         *gs_mol_sha,
              double         *qflx_sun, 
              double         *qflx_sha,
              veg_var_struct *veg_var)
{
    /* 局部变量 */
    double rppdry_sun;    // 阳叶潜在蒸发通过蒸腾的比例 [-]
    double rppdry_sha;    // 阴叶潜在蒸发通过蒸腾的比例 [-]
    double dryFrac = veg_var->dryFrac;
    double NetLAI = veg_var->NetLAI;
    double NetSAI = veg_var->NetSAI;
    double leafsun = veg_var->leaf_sun;
    double leafsha = veg_var->leaf_sha;
    // 计算转换因子
    double CF = pressure / (CONST_RGAS * thm) * 1.0e6;
    // 计算叶片总热导度 [m/s]
    double wtl = (NetLAI + NetSAI) * gb_mol;
    // potential latent energy flux [m/s]
    double efpot = air_density * wtl * (qsat_T - Qair_over) / CONST_RHOFW;
    
    if (CALC_TRANSP) {
        if ((efpot > 0.0) && (NetLAI > 0.0)) {
            // 计算阳叶蒸腾
            if (*gs_mol_sun > 0.0) {
                rppdry_sun = dryFrac / gb_mol * 
                             (leafsun / (1.0 / gb_mol + 1.0 / (*gs_mol_sun))) / 
                             NetLAI;
                *qflx_sun = efpot * rppdry_sun / CF;
            } else {
                *qflx_sun = 0.0;
            }
            // 计算阴叶蒸腾
            if (*gs_mol_sha > 0.0) {
                rppdry_sha = dryFrac / gb_mol * 
                             (leafsha / (1.0 / gb_mol + 1.0 / (*gs_mol_sha))) / 
                             NetLAI;
                *qflx_sha = efpot * rppdry_sha / CF;
            } else {
                *qflx_sha = 0.0;
            }
            
        } else {
            // 潜在通量为负或没有叶片
            *qflx_sun = 0.0;
            *qflx_sha = 0.0;
        }
        
    } 
    else {    
        // 反推阳叶气孔导度
        if (*qflx_sun > 0.0) {
            double denominator = efpot * dryFrac * leafsun - (*qflx_sun) * CF * NetLAI;
            if (fabs(denominator) > 1.0e-15) {
                *gs_mol_sun = gb_mol * (*qflx_sun) * CF * NetLAI / denominator;
                // 确保导度为非负
                if (*gs_mol_sun < 0.0) {
                    (*gs_mol_sun) = 0.0;
                }
            } 
            else {
                (*gs_mol_sun) = 0.0;
            }
        } 
        else {
            (*gs_mol_sun) = 0.0;
        }
        
        // 反推阴叶气孔导度
        if (*qflx_sha > 0.0) {
            double denominator = efpot * dryFrac * leafsha - (*qflx_sha) * CF * NetLAI;
            if (fabs(denominator) > 1.0e-15) {
                *gs_mol_sha = gb_mol * (*qflx_sha) * CF * NetLAI / denominator;

                if (*gs_mol_sha < 0.0) {
                    (*gs_mol_sha) = 0.0;
                }
            } 
            else {
                (*gs_mol_sha) = 0.0;
            }
        } 
        else {
            (*gs_mol_sha) = 0.0;
        }
    }
}
/******************************************************************************
 * @brief    Calculate temperature scaling factor.
 *****************************************************************************/
void spacF(double           *mat_VEG,
           double           *mat_A,
           double           *mat_RHS,
           double            qflx_sun, 
           double            qflx_sha,
           cell_data_struct *cell, 
           soil_con_struct  *soil_con, 
           veg_var_struct   *veg_var,
           veg_lib_struct   *veg_lib)
{
    // 局部变量
    double temp;               /* 临时变量，用于交换 f(sun) 和 f(sha) */
    double k_leaf_sun;         /* 阳叶导度 [m/s] */
    double k_leaf_sha;         /* 阴叶导度 [m/s] */
    double k_xyl;              /* 木质部导度 [m/s] */
    double k_root_soil;        /* 根-土壤导度 */
    double root_soil_pot;      /* 根-土壤水势差累加 */
    double Canopy_Upper = veg_lib->Canopy_Upper;
    double matric50 = veg_lib->matric50;
    double kcano_max = veg_lib->kcano_max;
    double kmax_sun = kcano_max;
    double kmax_sha = kcano_max;
    double kmax_xyl = kcano_max;
    double *zc_soil = soil_con->zc_soil;
    double *matric = cell->matric;
    double *hksr_int = cell->hksr_int;
    double NetSAI = veg_var->NetSAI;
    double leaf_sun = veg_var->leaf_sun;
    double leaf_sha = veg_var->leaf_sha;
    size_t i, j;
    size_t Nroot = cell->Nroot;
    double grav2[MAX_SOILS];
    for (i = 0 ; i < MAX_SOILS; i++) {
        grav2[i] = 0.0;
    }
    // 计算重力势
    double grav1 = Canopy_Upper;
    for (i = 0; i < Nroot; i++) {
        grav2[i] = zc_soil[i];
    }
    // 计算土壤相关项的累加
    k_root_soil = 0.0;
    root_soil_pot = 0.0;
    for (j = 0; j < Nroot; j++) {
        k_root_soil += hksr_int[j];
        root_soil_pot += hksr_int[j] * (mat_VEG[3] + grav2[j] - matric[j]);
    }
    // 计算各段的导度衰减函数
    double fsto1 = plc(mat_VEG[0], matric50);
    double fsto2 = plc(mat_VEG[1], matric50);
    double fx    = plc(mat_VEG[2], matric50);
    double fr    = plc(mat_VEG[3], matric50);
    // 计算各段的导数
    double dfsto1 = d1plc(mat_VEG[0], matric50, 2.0);
    double dfsto2 = d1plc(mat_VEG[1], matric50, 2.0);
    double dfx    = d1plc(mat_VEG[2], matric50, 2.0);
    double dfr    = d1plc(mat_VEG[3], matric50, 2.0);

    // 计算各段的实际导度
    k_leaf_sun = leaf_sun * kmax_sun * fx;
    k_leaf_sha = leaf_sha * kmax_sha * fx;
    k_xyl      = NetSAI * kmax_xyl / Canopy_Upper * fr;
    // 构建雅可比矩阵A
    mat_A[0] = -leaf_sun * kmax_sun * fx - qflx_sun * dfsto1;
    mat_A[1] = 0.0;                      
    mat_A[2] = leaf_sun * kmax_sun * fx;
    mat_A[3] = 0.0;
    mat_A[4] = 0.0;                                            
    mat_A[5] = -leaf_sha * kmax_sha * fx - qflx_sha * dfsto2;
    mat_A[6] = leaf_sha * kmax_sha * fx;
    mat_A[7] = 0.0;
    mat_A[8] = leaf_sun * kmax_sun * dfx * (mat_VEG[2] - mat_VEG[0]) +
               leaf_sun * kmax_sun * fx;
    mat_A[9] = leaf_sha * kmax_sha * dfx * (mat_VEG[2] - mat_VEG[1]) +                   
               leaf_sha * kmax_sha * fx;
    mat_A[10] = -leaf_sun * kmax_sun * dfx * (mat_VEG[2] - mat_VEG[0]) -               
                leaf_sun * kmax_sun * fx -
                leaf_sha * kmax_sha * dfx * (mat_VEG[2] - mat_VEG[1]) -
                leaf_sha * kmax_sha * fx -
                NetSAI * kmax_xyl / Canopy_Upper * fr;
    mat_A[11] = NetSAI * kmax_xyl / Canopy_Upper * fr;
    mat_A[12] = 0.0;                                                              
    mat_A[13] = 0.0;                                                                  
    mat_A[14] = NetSAI * kmax_xyl / Canopy_Upper * dfr * (mat_VEG[3] - mat_VEG[2] - grav1) +
                NetSAI * kmax_xyl / Canopy_Upper * fr;
    mat_A[15] = -NetSAI * kmax_xyl / Canopy_Upper * fr -                                
                 NetSAI * kmax_xyl / Canopy_Upper * dfr * 
                (mat_VEG[3] - mat_VEG[2] - grav1) - k_root_soil;
    // 构建右侧项RHS
    mat_RHS[0] = qflx_sun * fsto1 - k_leaf_sun * (mat_VEG[2] - mat_VEG[0]);
    
    // 阴叶节点
    mat_RHS[1] = qflx_sha * fsto2 - k_leaf_sha * (mat_VEG[2] - mat_VEG[1]);
    mat_RHS[2] = k_leaf_sun * (mat_VEG[2] - mat_VEG[0]) +
                 k_leaf_sha * (mat_VEG[2] - mat_VEG[1]) -
                 k_xyl * (mat_VEG[3] - mat_VEG[2] - grav1);
    
    mat_RHS[3] = k_xyl * (mat_VEG[3] - mat_VEG[2] - grav1) +
                 root_soil_pot;
    
    /* --- 特殊情况处理：阴叶LAI接近零 --- */
    if (leaf_sha < MIN_TOL_LAI) {
        temp = mat_RHS[0];
        mat_RHS[0] = mat_RHS[1];
        mat_RHS[1] = temp;
    }
}
/******************************************************************************
 * @brief    Calculates transpiration and returns corresponding vegwp.
 *****************************************************************************/
void getvegwp(double           *mat_VEG, 
              double            gb_mol, 
              double           *gs_mol_sun, 
              double           *gs_mol_sha,
              double            qsat_T, 
              double            Qair_over, 
              double           *transp,
              double            pressure,
              double            air_density,
              double            thm,
              cell_data_struct *cell, 
              soil_con_struct  *soil_con,
              veg_var_struct   *veg_var,
              veg_lib_struct   *veg_lib)
{
    /* 局部变量 */
    size_t i, j;
    size_t Nroot = cell->Nroot;
    double qflx_sun, qflx_sha;     /* 蒸腾通量 [m/s] */
    double fx;                     /* 木质部到叶的最大导度比例 [-] */
    double fr;                     /* 根到木质部的最大导度比例 [-] */
    double k_soil_root;            /* 土壤-根导度总和 */
    double soil_root_pot;          /* 导度加权土壤水势总和 */
    double *dz_soil = soil_con->dz_soil;
    double *matric = cell->matric;
    double leaf_sun = veg_var->leaf_sun;
    double leaf_sha = veg_var->leaf_sha;
    double Canopy_Upper = veg_lib->Canopy_Upper;
    double matric50 = veg_lib->matric50;
    double kcano_max = veg_lib->kcano_max;
    double *hksr_int = cell->hksr_int;
    double NetSAI = veg_var->NetSAI;
    double kmax_sun = kcano_max;
    double kmax_sha = kcano_max;
    double kmax_root = kcano_max;
    double grav2[MAX_SOILS];
    for (i = 0 ; i < MAX_SOILS; i++) {
        grav2[i] = 0.0;
    }
    // 计算重力势
    double grav1 = Canopy_Upper;
    for (i = 0; i < Nroot; i++) {
        grav2[i] = dz_soil[i];
    }
    
    bool CALC_TRANSP = true;
    get_qflx(CALC_TRANSP, gb_mol, 
             qsat_T, Qair_over, 
             pressure, air_density, 
             thm, gs_mol_sun, 
             gs_mol_sha, &qflx_sun,
             &qflx_sha, veg_var);
    
    // 计算土壤-根导度总和
    k_soil_root = 0.0;
    soil_root_pot = 0.0;
    for (j = 0; j < Nroot; j++) {
        k_soil_root += hksr_int[j];
        soil_root_pot += hksr_int[j] * (matric[j] - grav2[j]);
    }
    
    if (fabs(k_soil_root) == 0.0) {
        /* 无有效土壤导度，使用简单平均 */
        mat_VEG[3] = 0.0;
        for (j = 0; j < Nroot; j++) {
            mat_VEG[3] += (matric[j] - grav2[j]);
        }
        mat_VEG[3] /= (double) Nroot;
    }
    else {
        /* 考虑蒸腾拉力的根水势 */
        mat_VEG[3] = (soil_root_pot - qflx_sun - qflx_sha) / k_soil_root;
    }
    
    /* 计算根部导度衰减 */
    fr = plc(mat_VEG[3], matric50);
    
    if ((NetSAI > 0.0) && (fr > 0.0)) {
        /* 存在茎干且导度不为零，考虑水流阻力 */
        double k_root_eff = fr * kmax_root / Canopy_Upper * NetSAI;
        mat_VEG[2] = mat_VEG[3] - grav1 - (qflx_sun + qflx_sha) / k_root_eff;
    } 
    else {
        // 无茎干或完全栓塞，仅考虑重力
        mat_VEG[2] = mat_VEG[3] - grav1;
    }
    
    /* 计算木质部导度衰减 */
    fx = plc(mat_VEG[2], matric50);
    
    // 计算阴叶水势
    if ((leaf_sha > 0.0) && (fx > 0.0)) {
        double k_sha_eff = fx * kmax_sha * leaf_sha;
        mat_VEG[1] = mat_VEG[2] - qflx_sha / k_sha_eff;
    } 
    else {
        mat_VEG[1] = mat_VEG[2];
    }
    
    // 计算阳叶水势
    if ((leaf_sun > 0.0) && (fx > 0.0)) {
        double k_sun_eff = fx * kmax_sun * leaf_sun;
        mat_VEG[0] = mat_VEG[2] - qflx_sun / k_sun_eff;
    } 
    else {
        mat_VEG[0] = mat_VEG[2];
    }
    
    (*transp) = 0.0;
    for (j = 0; j < Nroot; j++) {
        (*transp) += hksr_int[j] * (matric[j] - mat_VEG[3] - grav2[j]);
    }
}