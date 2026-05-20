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
calc_stress(double            *bsun,
            double            *bsha,
            double             thm,
            double             RS_mol,
            double             qsat_T,
            double             Qair_over,
            double             pressure,
            double             air_density,
            energy_bal_struct *energy,
            cell_data_struct  *cell,
            soil_con_struct   *soil_con,
            veg_var_struct    *veg_var,
            veg_lib_struct    *veg_lib)
{
    extern parameters_struct param;
    // 初始化输出变量
    size_t i, j, k;
    size_t iter;
    bool flag = false;
    (*bsun) = 1.0;
    (*bsha) = 1.0;
    double qflx_sun = 0.0;
    double qflx_sha = 0.0;
    double vegwp[4]; // 水势：0-阳叶，1-阴叶，2-木质部，3-根部
    double dx[4]; // 水势更新量
    double mat_A[16]; // 雅可比矩阵
    double mat_RHS[4]; // 右侧项
    double matric50 = veg_lib->matric50; // 50%失水点的水势
    double leaf_sun = veg_var->leaf_sun;
    double leaf_shade = veg_var->leaf_shade;
    double conduct_max = veg_lib->conduct_max; // 最大导度
    double Canopy_Upper = veg_lib->Canopy_Upper;
    if (energy->aPAR_sunlit > 0.0) {
        vegwp[0] = 1.0;
        vegwp[1] = 0.0;
    }
    // Medlyn intercept of conductance-photosynthesis relationship [umol H2O]
    double gsmin = 200.0;
    // 复制以避免重写原始导度值
    double gs0sun = gsmin;
    double gs0sha = gsmin;
    
    /* 计算蒸腾需求（未受胁迫的潜在蒸腾） */
    bool CALC_TRANSP = true;
    getqflx(CALC_TRANSP, RS_mol, 
            qsat_T, Qair_over, 
            pressure, air_density, 
            thm, &gs0sun, 
            &gs0sha, &qflx_sun,
            &qflx_sha, veg_var);
    
    // 检查是否存在足够的叶面积和正的蒸腾需求
    if ((leaf_sun > MIN_VEG_LAI || leaf_shade > MIN_VEG_LAI) &&
        (qflx_sun > 0.0 || qflx_sha > 0.0)) {

        iter = 0;
        flag = false;
        
        do {

            iter++;
            // 计算通量残差
            spacF(vegwp, mat_RHS, qflx_sun, 
                  qflx_sha, matric50,
                  Canopy_Upper, conduct_max, 
                  cell, soil_con, veg_var);
            double sum = 0.0;
            for (i = 0; i < 4; i++) {
                sum += vegwp[i] * vegwp[i];
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
            
            // 计算雅可比矩阵 A
            spacA(mat_A, qflx_sun, qflx_sha, &flag);
            // 矩阵不可逆，退出迭代采用代数方法求解
            if (flag) {
                break;
            }
            
            // 根据叶面积情况计算dx
            if (leaf_sun > MIN_VEG_LAI && leaf_shade > MIN_VEG_LAI) {
                // 调用LAPACKE函数计算矩阵向量乘积 A * dx = RHS
                int info = LAPACKE_dgemv(LAPACK_COL_MAJOR, 'N', 
                                        4, 4, 1.0, mat_A, 4, mat_RHS, 1, 0.0, dx, 1);

                if (info != 0) {
                    flag = true;
                    break;  // 这里检查是有效的
                }
            } 
            else {
                // 简化为3x3系统
                memset(dx, 0, sizeof(dx));
                
                /* 根据哪个叶面积为零来调整计算 */
                if (leaf_shade > MIN_VEG_LAI) {
                    /* 阳叶LAI=0，求解阴叶-木质部-根系统 */
                    double A_3x3[9], B_3x3[3];
                    bool inv_flag = false;
                    int ipiv[3];
                    for (i = 0; i < 3; i++) {
                        for (j = 0; j < 3; j++) {
                            A_3x3[i * 3 + j] = mat_A[(1 + i) * 4 + (1 + j)];
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
                    /* 阴叶LAI=0，求解阳叶-木质部-根系统 */
                    int ipiv[3];
                    double A_3x3[9], B_3x3[3];
                    bool inv_flag = false;
                    
                    // 构建简化的3x3系统: [sun, xyl, root]
                    for (j = 0; j < 3; j++) {
                        for (i = 0; i < 3; i++) {
                            A_3x3[j * 3 + i] = mat_A[j * 4 + i];
                        }
                    }                   
                    B_3x3[0] = mat_RHS[0];
                    B_3x3[1] = mat_RHS[2];
                    B_3x3[2] = mat_RHS[3];

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
            
            // 限制步长，防止过大的变化
            double max_dx = fabs(dx[0]);
            for (j = 1; j < 4; j++) {
                if (fabs(dx[j]) > max_dx) {
                    max_dx = fabs(dx[j]);
                }
            }
            
            if (max_dx > 50000.0) {
                double scale = 50000.0 / max_dx;
                for (j = 0; j < 4; j++) {
                    dx[j] *= scale;
                }
            }
            
            // 更新水势
            if (leaf_sun > MIN_VEG_LAI && leaf_shade > MIN_VEG_LAI) {
                for (j = 0; j < 4; j++) {
                    vegwp[j] += dx[j];
                }
            } 
            else if (leaf_shade > MIN_VEG_LAI) {
                for (j = 0; j < 4; j++) {
                    vegwp[j] += dx[j];
                }
                vegwp[0] = vegwp[2];
            } 
            else {
                /* laisun > 0, laisha == 0 */
                vegwp[2] += dx[2];
                vegwp[3] += dx[3];
                vegwp[0] += dx[1];
                vegwp[1] = vegwp[2];
            }
            
            /* 检查水势变化是否足够小 */
            double dx_norm = sqrt(sum_of_squares(dx, 4));
            if (dx_norm < 1.0e-9) {
                break;
            }
            
            // 强制保证SPAC梯度朝向大气
            if (vegwp[2] > vegwp[3]) {
                vegwp[2] = vegwp[3];
            }
            if (vegwp[0] > vegwp[2]) {
                vegwp[0] = vegwp[2];
            }
            if (vegwp[1] > vegwp[2]) {
                vegwp[1] = vegwp[2];
            }
        } while (1);
        
    } else {
        flag = true; // 蒸腾通量为零
    }
    
    // 根据求解结果计算胁迫系数
    double qsun, qsha;
    double soilflux = 0.0;
    if (flag) {
        // 代数求解
        getvegwp(vegwp, RS_mol, &gs0sun, &gs0sha, 
                 qsat_T, Qair_over, &soilflux, 
                 Canopy_Upper, matric50, 
                 pressure, conduct_max, air_density, 
                 thm, cell, soil_con, veg_var);
        *bsun = plc(vegwp[0], matric50);
        *bsha = plc(vegwp[1], matric50);
    } 
    else {
        // 计算衰减后的实际通量
        qsun = qflx_sun * plc(vegwp[0], matric50);
        qsha = qflx_sha * plc(vegwp[1], matric50);
        
        // 反推受胁迫的气孔导度
        CALC_TRANSP = false;
        getqflx(CALC_TRANSP, RS_mol, 
                qsat_T, Qair_over, 
                pressure, air_density, 
                thm, &gs0sun, 
                &gs0sha, &qsun,
                &qsha, veg_var);
        
        if (qflx_sun > 0.0) {
            *bsun = gs0sun / gsmin;
        } 
        else {
            *bsun = plc(vegwp[0], matric50);
        }
        
        if (qflx_sha > 0.0) {
            *bsha = gs0sha / gsmin;
        } 
        else {
            *bsha = plc(vegwp[1], matric50);
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
    if (energy->aPAR_sunlit <= 0.0) {
        gs0sun = (*bsun) * gsmin;
        gs0sha = (*bsha) * gsmin;
        getvegwp(vegwp, RS_mol, &gs0sun, &gs0sha, 
                 qsat_T, Qair_over, &soilflux, 
                 Canopy_Upper, matric50, 
                 pressure, conduct_max, air_density, 
                 thm, cell, soil_con, veg_var);
        if (soilflux < 0.0) {
            soilflux = 0.0;
        }
        cell->transp = soilflux;
    }
}

/******************************************************************************
 * @brief    Calculate temperature scaling factor.
 *****************************************************************************/
void get_qflx(bool            CALC_TRANSP,
              double          RS_mol,
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
    double leafshade = veg_var->leaf_shade;
    // 计算转换因子
    double cf = pressure / (CONST_RGAS * thm) * 1.0e6;
    // 计算叶片总热导度 [m/s]
    double wtl = (NetLAI + NetSAI) * RS_mol;
    // potential latent energy flux [kg/m2/s]
    double efpot = air_density * wtl * (qsat_T - Qair_over);
    
    if (CALC_TRANSP) {
        if ((efpot > 0.0) && (NetLAI > 0.0)) {
            // 计算阳叶蒸腾
            if (*gs_mol_sun > 0.0) {
                rppdry_sun = dryFrac / RS_mol * 
                             (leafsun / (1.0 / RS_mol + 1.0 / (*gs_mol_sun))) / 
                             NetLAI;
                *qflx_sun = efpot * rppdry_sun / cf;
            } else {
                *qflx_sun = 0.0;
            }
            // 计算阴叶蒸腾
            if (*gs_mol_sha > 0.0) {
                rppdry_sha = dryFrac / RS_mol * 
                             (leafshade / (1.0 / RS_mol + 1.0 / (*gs_mol_sha))) / 
                             NetLAI;
                *qflx_sha = efpot * rppdry_sha / cf;
            } else {
                *qflx_sha = 0.0;
            }
            
        } else {
            // 潜在通量为负或没有叶片
            *qflx_sun = 0.0;
            *qflx_sha = 0.0;
        }
        
    } else {    
        // 反推阳叶气孔导度
        if (*qflx_sun > 0.0) {
            double denominator = efpot * dryFrac * leafsun - (*qflx_sun) * cf * NetLAI;
            if (fabs(denominator) > 1.0e-15) {
                *gs_mol_sun = RS_mol * (*qflx_sun) * cf * NetLAI / denominator;
                // 确保导度为非负
                if (*gs_mol_sun < 0.0) *gs_mol_sun = 0.0;
            } else {
                *gs_mol_sun = 0.0;
            }
        } else {
            *gs_mol_sun = 0.0;
        }
        
        // 反推阴叶气孔导度
        if (*qflx_sha > 0.0) {
            double denominator = efpot * dryFrac * leafshade - (*qflx_sha) * cf * NetLAI;
            if (fabs(denominator) > 1.0e-15) {
                *gs_mol_sha = RS_mol * (*qflx_sha) * cf * NetLAI / denominator;

                if (*gs_mol_sha < 0.0) {
                    *gs_mol_sha = 0.0;
                }
            } 
            else {
                *gs_mol_sha = 0.0;
            }
        } else {
            *gs_mol_sha = 0.0;
        }
    }
}
/******************************************************************************
 * @brief    Calculate temperature scaling factor.
 *****************************************************************************/
void spacF(double           *vegwp, 
           double           *mat_RHS,
           double            qflx_sun, 
           double            qflx_sha,
           double            matric50,
           double            Canopy_Upper,
           double            conduct_max,
           cell_data_struct *cell, 
           soil_con_struct  *soil_con, 
           veg_var_struct   *veg_var)
{
    /* 局部变量 */
    double grav1;              /* 地表到冠层顶的重力势 [mm H2O] */
    double fsto1;              /* 阳叶蒸腾衰减函数 [-] */
    double fsto2;              /* 阴叶蒸腾衰减函数 [-] */
    double fx;                 /* 木质部到叶的最大导度比例 [-] */
    double fr;                 /* 根到木质部的最大导度比例 [-] */
    double temp;               /* 临时变量，用于交换 f(sun) 和 f(sha) */
    double k_leaf_sun;         /* 阳叶导度 [mm/s] */
    double k_leaf_sha;         /* 阴叶导度 [mm/s] */
    double k_xyl;              /* 木质部导度 [mm/s] */
    double k_root_soil;        /* 根-土壤导度 */
    double root_soil_pot;      /* 根-土壤水势差累加 */
    double kmax_sun = conduct_max;
    double kmax_sha = conduct_max;
    double kmax_xyl = conduct_max;
    double kmax_root = conduct_max;
    double *dz_soil = soil_con->dz_soil;
    double *matric = cell->matric;
    double *hksr_int = cell->hksr_int;
    double NetSAI = veg_var->NetSAI;
    double leaf_sun = veg_var->leaf_sun;
    double leaf_shade = veg_var->leaf_shade;
    size_t i, j;
    size_t Nsoil = cell->Nsoil;
    double grav2[MAX_SOILS];
    for (i = 0 ; i < MAX_SOILS; i++) {
        grav2[i] = 0.0;
    }
    // 计算重力势
    grav1 = Canopy_Upper * MM_PER_M;
    for (i = 0; i < Nsoil; i++) {
        grav2[i] = dz_soil[i] * MM_PER_M;
    }
    // 计算各段的导度衰减函数
    fsto1 = plc(vegwp[0], matric50);     // 阳叶导度衰减
    fsto2 = plc(vegwp[1], matric50);     // 阴叶导度衰减
    fx    = plc(vegwp[2], matric50);     // 木质部导度衰减
    fr    = plc(vegwp[3], matric50);     // 根部导度衰减
    
    // 计算各段的实际导度
    k_leaf_sun = leaf_sun * kmax_sun * fx;           // 阳叶到木质部导度
    k_leaf_sha = leaf_shade * kmax_sha * fx;           // 阴叶到木质部导度
    k_xyl      = NetSAI * kmax_xyl / Canopy_Upper * fr;      // 木质部导度
    
    mat_RHS[0] = qflx_sun * fsto1 - k_leaf_sun * (vegwp[2] - vegwp[0]);
    
    /* --- 阴叶节点 --- */
    mat_RHS[1] = qflx_sha * fsto2 - k_leaf_sha * (vegwp[2] - vegwp[1]);
    mat_RHS[2] = k_leaf_sun * (vegwp[2] - vegwp[0])
           + k_leaf_sha * (vegwp[2] - vegwp[1])
           - k_xyl * (vegwp[3] - vegwp[2] - grav1);
    
    /* 计算土壤相关项的累加 */
    k_root_soil = 0.0;
    root_soil_pot = 0.0;
    for (j = 0; j < Nsoil; j++) {
        k_root_soil += hksr_int[j];
        root_soil_pot += hksr_int[j] * (vegwp[3] + grav2[j] - matric[j]);
    }
    
    mat_RHS[3] = k_xyl * (vegwp[3] - vegwp[2] - grav1)
            + root_soil_pot;
    
    /* --- 特殊情况处理：阴叶LAI接近零 --- */
    if (leaf_shade < MIN_VEG_LAI) {
        temp = mat_RHS[0];
        mat_RHS[0] = mat_RHS[1];
        mat_RHS[1] = temp;
    }
}
/******************************************************************************
 * @brief    Calculate temperature scaling factor.
 *****************************************************************************/
void getvegwp(double           *vegwp, 
              double            RS_mol, 
              double           *gs_mol_sun, 
              double           *gs_mol_sha,
              double            qsat_T, 
              double            Qair_over, 
              double           *soilflux,
              double            Canopy_Upper,
              double            matric50,
              double            pressure,
              double            conduct_max,
              double            air_density, 
              double            thm,
              cell_data_struct *cell, 
              soil_con_struct  *soil_con,
              veg_var_struct   *veg_var)
{
    /* 局部变量 */
    double qflx_sun, qflx_sha;    /* 蒸腾通量 [kg/m2/s] */
    double fx;                     /* 木质部到叶的最大导度比例 [-] */
    double fr;                     /* 根到木质部的最大导度比例 [-] */
    double sum_k_soil_root;        /* 土壤-根导度总和 */
    double sum_k_soil_root_smp;    /* 导度加权土壤水势总和 */
    double *dz_soil = soil_con->dz_soil;
    double *matric = cell->matric;
    double leaf_sun = veg_var->leaf_sun;
    double leaf_shade = veg_var->leaf_shade;
    size_t i, j;
    size_t Nsoil = cell->Nsoil;
    double *hksr_int = cell->hksr_int;
    double NetSAI = veg_var->NetSAI;
    double kmax_sun = conduct_max;
    double kmax_sha = conduct_max;
    double kmax_xyl = conduct_max;
    double kmax_root = conduct_max;
    double grav2[MAX_SOILS];
    for (i = 0 ; i < MAX_SOILS; i++) {
        grav2[i] = 0.0;
    }
    // 计算重力势
    double grav1 = Canopy_Upper * MM_PER_M;
    for (i = 0; i < Nsoil; i++) {
        grav2[i] = dz_soil[i] * MM_PER_M;
    }
    
    bool CALC_TRANSP = true;
    getqflx(CALC_TRANSP, RS_mol, 
            qsat_T, Qair_over, 
            pressure, air_density, 
            thm, &gs_mol_sun, 
            &gs_mol_sha, &qflx_sun,
            &qflx_sha, veg_var);
    
    // 计算土壤-根导度总和
    sum_k_soil_root = 0.0;
    sum_k_soil_root_smp = 0.0;
    for (j = 0; j < Nsoil; j++) {
        sum_k_soil_root += hksr_int[j];
        sum_k_soil_root_smp += hksr_int[j] * (matric[j] - grav2[j]);
    }
    
    if (fabs(sum_k_soil_root) < 1.0e-15) {
        /* 无有效土壤导度，使用简单平均 */
        vegwp[3] = 0.0;
        for (j = 0; j < Nsoil; j++) {
            vegwp[3] += (matric[j] - grav2[j]);
        }
        vegwp[3] /= (double) Nsoil;
    }
    else {
        /* 考虑蒸腾拉力的根水势 */
        vegwp[3] = (sum_k_soil_root_smp - qflx_sun - qflx_sha) / sum_k_soil_root;
    }
    
    /* 计算根部导度衰减 */
    fr = plc(vegwp[3], matric50);
    
    if ((NetSAI > 0.0) && (fr > 0.0)) {
        /* 存在茎干且导度不为零，考虑水流阻力 */
        double k_root_eff = fr * kmax_root / Canopy_Upper * NetSAI;
        vegwp[2] = vegwp[3] - grav1 - (qflx_sun + qflx_sha) / k_root_eff;
    } 
    else {
        // 无茎干或完全栓塞，仅考虑重力
        vegwp[2] = vegwp[3] - grav1;
    }
    
    /* 计算木质部导度衰减 */
    fx = plc(vegwp[2], matric50);
    
    // 计算阴叶水势
    if ((leaf_shade > 0.0) && (fx > 0.0)) {
        double k_sha_eff = fx * kmax_sha * leaf_shade;
        vegwp[1] = vegwp[2] - qflx_sha / k_sha_eff;
    } 
    else {
        vegwp[1] = vegwp[2];
    }
    
    // 计算阳叶水势
    if ((leaf_sun > 0.0) && (fx > 0.0)) {
        double k_sun_eff = fx * kmax_sun * leaf_sun;
        vegwp[0] = vegwp[2] - qflx_sun / k_sun_eff;
    } 
    else {
        vegwp[0] = vegwp[2];
    }
    
    *soilflux = 0.0;
    for (j = 0; j < Nsoil; j++) {
        *soilflux += hksr_int[j] * (matric[j] - vegwp[3] - grav2[j]);
    }
}