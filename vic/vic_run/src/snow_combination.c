/******************************************************************************
* @section DESCRIPTION
*
* Snowpack layer combination process
* Update snow ice, snow water, snow thickness, snow temperature.
******************************************************************************/

#include <vic_run.h>

/******************************************************************************
* @brief    Calculate snow depth via compaction due to destructive metamorphism, 
*           overburden, & melt.
******************************************************************************/
void 
snow_combination(double            dz_soil,
                 cell_data_struct *cell,
                 snow_data_struct *snow)
{

    size_t      i, j, k;
    size_t      lindex;
    size_t      index;
    double      old_Nsnow;
    double      pack_comb;
    double      tmp_ice;
    double      tmp_liq;

    /* initialize */
    old_Nsnow = snow->Nsnow;
    tmp_ice = 0.0;
    tmp_liq = 0.0;

    // 定义指针指向结构体中的数组
    double *ice = cell->ice;
    double *liq = cell->liq;
    double *dz_snow = snow->dz_snow;
    double *pack_T = snow->pack_T;
    double *pack_ice = snow->pack_ice;
    double *pack_liq = snow->pack_liq;

    for (i = 0; i < old_Nsnow; i++) {
        if (pack_ice[i] < 0.1) {
            // 如果不是底层雪层向下合并
            if (i != 0) {
                pack_liq[i - 1] += pack_liq[i];
                pack_ice[i - 1] += pack_ice[i];
                dz_snow[i - 1] += dz_snow[i];
            }
            // 如果是底层雪层向上合并
            else {
                // 如果包含2至3雪层
                if (snow->Nsnow > 1) {
                    pack_liq[i + 1] += pack_liq[i];
                    pack_ice[i + 1] += pack_ice[i];
                    dz_snow[i + 1] += dz_snow[i];
                }
                // 只有1层
                else {
                    if (pack_ice[i] > 0.0) {
                        snow->pack_comb = pack_liq[i];
                        snow->swq = pack_ice[i];
                        snow->snow_depth = dz_snow[i];
                    }
                    else {
                        snow->pack_comb = pack_liq[i] + pack_ice[i];
                        if (snow->pack_comb < 0.0) {
                            ice[0] += snow->pack_comb / (dz_soil * MM_PER_M);
                            snow->pack_comb = 0.0;
                        }
                        snow->swq = 0.0;
                        snow->snow_depth = 0.0;
                    }
                    pack_liq[i] = 0.0;
                    pack_ice[i] = 0.0;
                    dz_snow[i] = 0.0;
                }
            }   // if(i != 0)
            if (i < snow->Nsnow && snow->Nsnow > 1) {
                for (j = i; j < snow->Nsnow - 1; j++) {
                    pack_T[j] = pack_T[j + 1];
                    pack_liq[j] = pack_liq[j + 1];
                    pack_ice[j] = pack_ice[j + 1];
                    dz_snow[j] = dz_snow[j + 1];
                }
            }
            snow->Nsnow -= 1;
        }
    }   // i loop

    if (ice[0] < 0.0) {
        liq[0] += ice[0];
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
        snow->pack_transp = tmp_liq;
        if (snow->swq <= 0.0) {
            snow->snow_depth = 0.0;
        }
    }

    /* check the snow depth, snow layers combined */
    if (snow->Nsnow > 1) {
        old_Nsnow = snow->Nsnow;
        lindex = 0;
        for (i = 0; i < old_Nsnow; i++) {
            if (dz_snow[i] < snow->snow_thresholds[lindex]) {
                if (i == snow->Nsnow - 1) {
                    index = i - 1;
                }
                else if (i == 0) {
                    index = i + 1;
                }
                else {
                    index = i - 1;
                }
                if (dz_snow[i + 1] + dz_snow[i] < 
                            dz_snow[i - 1] + dz_snow[i]) {
                    index = i + 1;
                }
                if (index > i) {
                    j = index;
                    k = i;
                }
                else {
                    j = i;
                    k = index;
                }
                /* update combined snow water & temperature */
                update_surface_fluxes(&dz_snow[j], 
                                      &pack_liq[j],
                                      &pack_ice[j],
                                      &pack_T[j],
                                      dz_snow[k], 
                                      pack_liq[k],
                                      pack_ice[k],
                                      pack_T[k]);
                /* Now shift all elements above this down one. */
                if (j + 1 < snow->Nsnow) {
                    for (k = j + 1; k < snow->Nsnow - 1; k++) {
                        pack_T[k] = pack_T[k + 1];
                        pack_ice[k] = pack_ice[k + 1];
                        pack_liq[k] = pack_liq[k + 1];
                        dz_snow[k] = dz_snow[k + 1];
                    }
                }
                /* Decrease the number of snow layers */
                snow->Nsnow -= 1;
                if (snow->Nsnow <= 1) {
                    return;
                }
            }
            else {
                lindex++;
            }
        }
    }
}