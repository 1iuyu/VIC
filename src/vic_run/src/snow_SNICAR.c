/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine calculate albedo of snow containing impurities and the 
 * evolution of snow effective radius.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief  
 * Compute the albedo of snow containing impurities and the evolution of 
 * snow effective radius.
 *****************************************************************************/
int
snow_SNICAR(size_t             spectral_band,
            double             step_dt,
            double             coszen,
            double             snowfall,
            energy_bal_struct *energy,
            cell_data_struct  *cell,
            snow_data_struct  *snow,
			soil_con_struct   *soil_con)
{
    double zc_wave[5] = {0.5, 0.85, 1.1, 1.35, 3.25};
    // 高斯积分点（弧度）
    static const double difgaus_pt[8] = {
        0.9894009, 0.9445750, 0.8656312, 0.7554044,
        0.6178762, 0.4580168, 0.2816036, 0.0950125
    };
    // 高斯权重
    static const double difgaus_wt[8] = {
        0.0271525, 0.0622535, 0.0951585, 0.1246290,
        0.1495960, 0.1691565, 0.1826034, 0.1894506
    };
    // 非球形雪粒不对称因子参数化系数
    static const double g_wvl[8] = {
        0.25, 0.70, 1.41, 1.90, 2.50, 3.50, 4.00, 5.00
    };
    static const double g_b0[7] = {
        9.76029E-1, 9.67798E-1, 1.00111, 1.00224,
        9.64295E-1, 9.97475E-1, 9.97475E-1
    };
    static const double g_b1[7] = {
        5.21042E-1, 4.96181E-1, 1.83711E-1, 1.37082E-1,
        5.50598E-2, 8.48743E-2, 8.48743E-2
    };
    static const double g_b2[7] = {
        -2.66792E-4, 1.14088E-3, 2.37011E-4, -2.35905E-4,
        8.40449E-4, -4.71484E-4, -4.71484E-4
    };
    size_t i;
    bool virtual_flag;
    double band_wgt[5];
    double DIR_wgt[5] = {1.0, 0.494, 0.181, 0.121, 0.204};
    double DFS_wgt[5] = {1.0, 0.586, 0.202, 0.109, 0.103};
    // 检查是否需要进行辐射传输计算
    if (coszen > 0.0 && snow->swq > MIN_SNOWMASS) {
        size_t tmp_Nsnow = snow->Nsnow;
        if (tmp_Nsnow == 0) {
            virtual_flag = true;
            tmp_Nsnow = 1;

        }
        else {
            virtual_flag = false;
        }
        if (spectral_band == BAND_DIR) {
            for (i =0; i < 5; i++) {
                band_wgt[i] = DIR_wgt[i];
            }
        }
        else if (spectral_band == BAND_DFS) {
            for (i = 0; i < 5; i++) {
                band_wgt[i] = DFS_wgt[i];
            }
        }
        double exp_min = exp(-10.0);
        double mu_not = max(coszen, 0.01);
        for (size_t i = 0; i < 5; i++) {
            size_t err_count = 0;
            // 设置光学性质
            for (i = 0; i < tmp_Nsnow; i++) {

            }
        }

    }
    else if (coszen > 0.0 && snow->swq > 0.0) {


    } 
    else {
        energy->AlbedoSnowDir[BAND_VIS] = 0.0;
        energy->AlbedoSnowDfs[BAND_NIR] = 0.0;
    }

    return (0);
}
