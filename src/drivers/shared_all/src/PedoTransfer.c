/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine calculate the hydrological properties of soil moisture 
 * based on different soil compositions
 *****************************************************************************/

#include "vic_driver_shared_all.h"

/******************************************************************************
 * @brief    Calculate soil water infiltration based on different soil composition.
 *****************************************************************************/
void
PedoTransfer(soil_con_struct *soil_con)
{
    double *sand = soil_con->sand_node;
    double *clay = soil_con->clay_node;
    double *orgm = soil_con->organic_node;
    double *Wpwp_node = soil_con->Wpwp_node;
    double *Wsat_node = soil_con->Wsat_node;
    double *expt_node = soil_con->expt_node;
    double *Ksat_node = soil_con->Ksat_node;

    for (size_t i = 0; i < soil_con->Nbedrock; i++) {

        // 输入参数检查与默认值设置
        if (sand[i] <= 0 || clay[i] <= 0) {
            sand[i] = 0.41;
            clay[i] = 0.18;
        }
        
        if (orgm[i] <= 0) {
            orgm[i] = 0.0;
        }
        
        // 计算1500 kPa (永久萎蔫点) 含水量
        double theta_1500t = -0.024 * sand[i] + 0.487 * clay[i] + 0.006 * orgm[i] +
                            0.005 * sand[i] * orgm[i] - 0.013 * clay[i] * orgm[i] + 
                            0.068 * sand[i] * clay[i] + 0.031;
        
        double theta_1500 = theta_1500t + 0.14 * theta_1500t - 0.02;
        theta_1500 = max(theta_1500, pow(10, -5.0));
        
        // 计算33 kPa (田间持水量) 含水量
        double theta_33t = -0.251 * sand[i] + 0.195 * clay[i] + 0.011 * orgm[i] + 
                        0.006 * sand[i] * orgm[i] - 0.027 * clay[i] * orgm[i] + 
                        0.452 * sand[i] * clay[i] + 0.299;
        
        double theta_33 = theta_33t + 1.283 * theta_33t * theta_33t - 
                        0.374 * theta_33t - 0.015;
        theta_33 = max(theta_33, pow(10, -3.0));
        
        // 计算饱和含水量与33 kPa含水量的差值
        double theta_s33t = 0.278 * sand[i] + 0.034 * clay[i] + 0.022 * orgm[i] - 
                            0.018 * sand[i] * orgm[i] - 0.027 * clay[i] * orgm[i] - 
                            0.584 * sand[i] * clay[i] + 0.078;
        
        double theta_s33 = theta_s33t + (0.636 * theta_s33t - 0.107);
        
        // 计算饱和基质势
        double psi_et = -21.67 * sand[i] - 27.93 * clay[i] - 81.97 * theta_s33 + 
                        71.12 * sand[i] * theta_s33 + 8.29 * clay[i] * theta_s33 + 
                        14.05 * sand[i] * clay[i] + 27.16;
        
        double psi_e = psi_et + 0.02 * psi_et * psi_et - 0.113 * psi_et - 0.7;
        
        // 赋值基础变量
        double smcwlt = theta_1500;
        double smcref = theta_33;
        double smcmax = theta_33 + theta_s33 - 0.097 * sand[i] + 0.043;
        
        // 计算孔隙大小分布指数
        double bexp = 3.816712826 / (log(theta_33) - log(theta_1500));
        
        // 饱和基质势
        double psisat = psi_e;
        
        // 饱和导水率
        double dksat = 1930.0 * pow((smcmax - theta_33), (3.0 - 1.0 / bexp));
        
        // 边界值限制
        psisat = max(0.1, psisat);
        psisat = 0.101997 * psisat;
        dksat = dksat / 3600000.0;
        // double dwsat = dksat * psisat * bexp / smcmax;
        
        // 应用物理约束范围
        Wsat_node[i] = max(0.32, min(smcmax, 0.50));
        //Wfc_node[i] = max(0.17, min(smcref, Wsat_node[i]));
        Wpwp_node[i] = max(0.01, min(smcwlt, smcref));
        expt_node[i] = max(2.50, min(bexp, 12.0));
        // psisat_node[i] = max(0.03, min(psisat, 1.00));
        Ksat_node[i] = max(5e-7, min(dksat, 1e-5));
        // Dsat_node[i] = max(1e-6, min(dwsat, 3e-5));
        // bubble_node[i] = 0.32 * expt_node[i] + 4.3;
    }
}