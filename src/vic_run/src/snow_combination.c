/******************************************************************************
* @section DESCRIPTION
*
* Snowpack layer combination process
* Update snow ice, snow water, snow thickness, snow temperature.
******************************************************************************/

#include "vic_run.h"

/******************************************************************************
* @brief    Calculate snow depth via compaction due to destructive metamorphism, 
*           overburden, & melt.
******************************************************************************/
void 
snow_combination(double            dz_soil,
                 cell_data_struct *cell,
                 snow_data_struct *snow)
{
    /* initialize */
    size_t i, j, k;
    size_t lidx;
    double tmp_ice = 0.0;
    double tmp_liq = 0.0;
    double pack_comb = 0.0;
    double pack_transp = 0.0;

    // 定义指针指向结构体中的数组
    double *ice = cell->ice;
    double *liq = cell->liq;
    double *dz_snow = snow->dz_snow;
    double *pack_T = snow->pack_T;
    double *pack_ice = snow->pack_ice;
    double *pack_liq = snow->pack_liq;

    i = 0;
    while (i < snow->Nsnow) {
        if (pack_ice[i] < 0.1) {
            // 如果是顶层（i == 0），将下层合并到顶层
            if (i == 0) {
                if (snow->Nsnow > 1) {
                    // 将第1层(i+1)合并到第0层(i)
                    update_snow_fluxes(&dz_snow[0],
                                       &pack_liq[0],
                                       &pack_ice[0],
                                       &pack_T[0],
                                       dz_snow[1],
                                       pack_liq[1],
                                       pack_ice[1],
                                       pack_T[1]);
                    
                    // 将下面所有层向上移动
                    for (j = i+1; j < snow->Nsnow - 1; j++) {
                        pack_T[j] = pack_T[j + 1];
                        pack_liq[j] = pack_liq[j + 1];
                        pack_ice[j] = pack_ice[j + 1];
                        dz_snow[j] = dz_snow[j + 1];
                    }
                    
                    // 清空最后一个元素
                    j = snow->Nsnow - 1;
                    pack_T[j] = 0.0;
                    pack_liq[j] = 0.0;
                    pack_ice[j] = 0.0;
                    dz_snow[j] = 0.0;
                    
                    snow->Nsnow--;
                }
                else {
                    // 只有一层雪的处理
                    if (pack_ice[i] > 0.0) {
                        pack_comb += pack_liq[i];
                        snow->swq = pack_ice[i];
                        snow->snow_depth = dz_snow[i];
                    } 
                    else {
                        pack_comb += pack_liq[i] + pack_ice[i];
                        if (pack_comb < 0.0) {
                            ice[0] += pack_comb / (dz_soil * MM_PER_M);
                            pack_comb = 0.0;
                        }
                        snow->swq = 0.0;
                        snow->snow_depth = 0.0;
                    }
                    
                    // 清空唯一雪层
                    pack_T[i] = 0.0;
                    pack_liq[i] = 0.0;
                    pack_ice[i] = 0.0;
                    dz_snow[i] = 0.0;
                    snow->Nsnow--;
                }
                
                // 继续检查合并后的顶层（索引0）
                continue;
            }
            // 如果不是顶层（i > 0），将本层合并到上层
            else {
                update_snow_fluxes(&dz_snow[i-1],
                                   &pack_liq[i-1],
                                   &pack_ice[i-1],
                                   &pack_T[i-1],
                                   dz_snow[i],
                                   pack_liq[i],
                                   pack_ice[i],
                                   pack_T[i]);
                
                // 将下面所有层向上移动
                for (j = i; j < snow->Nsnow - 1; j++) {
                    pack_T[j] = pack_T[j + 1];
                    pack_liq[j] = pack_liq[j + 1];
                    pack_ice[j] = pack_ice[j + 1];
                    dz_snow[j] = dz_snow[j + 1];
                }
                
                // 清空最后一个元素
                j = snow->Nsnow - 1;
                pack_T[j] = 0.0;
                pack_liq[j] = 0.0;
                pack_ice[j] = 0.0;
                dz_snow[j] = 0.0;
                
                snow->Nsnow--;
                
                // 回退检查合并后的上层
                i--;
                continue;
            }
        }
        i++;
    }

    snow->pack_comb = pack_comb;
    
    if (ice[0] < 0.0) {
        liq[0] += ice[0] * CONST_RHOICE / CONST_RHOFW;
        ice[0] = 0.0;
    }
    /* if no longer multi-layer */
    if (snow->Nsnow == 0) {
        return;
    }

    snow->swq = 0.0;
    snow->snow_depth = 0.0;

    for (i = 0; i < snow->Nsnow; i++) {
        snow->swq += pack_ice[i] + pack_liq[i];
        snow->snow_depth += dz_snow[i];
        tmp_ice += pack_ice[i];
        tmp_liq += pack_liq[i];
    }

    if (snow->snow_depth < 0.025 && snow->Nsnow > 0) {
        snow->Nsnow = 0;
        snow->swq = tmp_ice;
        pack_transp += tmp_liq;
        if (snow->swq <= 0.0) {
            snow->snow_depth = 0.0;
        }
        for (i = 0; i < MAX_SNOWS; i++) {
            pack_T[i] = 0.0;
            pack_liq[i] = 0.0;
            pack_ice[i] = 0.0;
            dz_snow[i] = 0.0;
        }
    }
    snow->pack_transp = pack_transp;
    
    /* check the snow depth, snow layers combined */
    if (snow->Nsnow > 1) {

        lidx = 0;

        i = 0;
        while (i < snow->Nsnow) {

            if (dz_snow[i] < snow->snow_thresholds[lidx]) {
                size_t source;
                size_t target;         
                // 选择源层和目标层
                if (i == 0) {
                    target = 0;
                    source = 1;
                }
                else {
                    target = i - 1;  // 目标层（上层）
                    source = i;      // 源层（当前薄层）
                }

                update_snow_fluxes(&dz_snow[target],
                                   &pack_liq[target],
                                   &pack_ice[target],
                                   &pack_T[target],
                                   dz_snow[source],
                                   pack_liq[source],
                                   pack_ice[source],
                                   pack_T[source]);

                // 删除源层，将后面的层前移
                for (j = source; j < snow->Nsnow - 1; j++) {
                    pack_T[j] = pack_T[j + 1];
                    pack_ice[j] = pack_ice[j + 1];
                    pack_liq[j] = pack_liq[j + 1];
                    dz_snow[j] = dz_snow[j + 1];
                }

                // 清空最后一层
                j = snow->Nsnow - 1;
                pack_T[j] = 0.0;
                pack_ice[j] = 0.0;
                pack_liq[j] = 0.0;
                dz_snow[j] = 0.0;

                snow->Nsnow--;

                if (snow->Nsnow <= 1) {
                    break;
                }

                /* 回到目标层重新检查 */
                i = target;
                lidx = target;

                continue;
            }

            lidx++;
            i++;
        }
    }
}