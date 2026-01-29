/******************************************************************************
* @section DESCRIPTION
*
* Snowpack layer division process
* Update snow ice, snow water, snow thickness, snow temperature.
******************************************************************************/

#include <vic_run.h>

/******************************************************************************
* @brief    Calculate snow depth via compaction due to destructive metamorphism, 
*           overburden, & melt.
******************************************************************************/
void 
snow_division(snow_data_struct *snow)

{              

    size_t      i, nidx;
    double      grad_temp;
    double      extra_frac;
    double      tmp_pack_ice[MAX_SNOWS];
    double      tmp_pack_liq[MAX_SNOWS];
    double      tmp_pack_T[MAX_SNOWS];
    double      tmp_depth[MAX_SNOWS];

    double total_depth = 0.;
    double extra_ice = 0.;
    double extra_liq = 0.;
    /* 初始化临时数组 */
    for (i = 0; i < MAX_SNOWS; i++) {
        tmp_pack_ice[i] = 0.;
        tmp_pack_liq[i] = 0.;
        tmp_pack_T[i] = 0.;
        tmp_depth[i] = 0.;
    }

    for (i = 0; i < snow->Nsnow; i++) {
        nidx = snow->Nsnow - 1;
        tmp_depth[i] = snow->dz_snow[nidx];
        tmp_pack_ice[i] = snow->pack_ice[nidx];
        tmp_pack_liq[i] = snow->pack_liq[nidx];
        tmp_pack_T[i] = snow->pack_T[nidx];
    }
    size_t tmp_Nsnow = snow->Nsnow;

    if (tmp_Nsnow == 1) {
        if (tmp_depth[0] > 0.05) {
            tmp_Nsnow = 2;
            tmp_depth[0] /= 2.0;
            tmp_pack_ice[0] /= 2.0;
            tmp_pack_liq[0] /= 2.0;
            tmp_depth[1] = tmp_depth[0];
            tmp_pack_ice[1] = tmp_pack_ice[0];
            tmp_pack_liq[1] = tmp_pack_liq[0];
            tmp_pack_T[1] = tmp_pack_T[0];
        }
    }
    if (tmp_Nsnow > 1) {
        if (tmp_depth[0] > 0.05) {
            total_depth = tmp_depth[0] - 0.05;
            extra_frac = total_depth / tmp_depth[0];
            extra_ice = extra_frac * tmp_pack_ice[0];
            extra_liq = extra_frac * tmp_pack_liq[0];

            extra_frac = 0.05 / tmp_depth[0];
            tmp_pack_ice[0] *= extra_frac;
            tmp_pack_liq[0] *= extra_frac;
            tmp_depth[0] = 0.05;
        }

        /* update combined snow water & temperature */
        update_surface_fluxes(&tmp_depth[1], 
                              &tmp_pack_liq[1],
                              &tmp_pack_ice[1],
                              &tmp_pack_T[1],
                              total_depth, 
                              extra_liq,
                              extra_ice,
                              tmp_pack_T[0]);

        /* subdivide the second snow layer */
        if (tmp_Nsnow <= 2 && tmp_depth[1] > 0.20) {
            tmp_Nsnow = 3;
            grad_temp = tmp_pack_T[0] - tmp_pack_T[1] / 
                             ((tmp_depth[0] + tmp_depth[1]) / 2.0);
            tmp_depth[1] /= 2.0;
            tmp_pack_ice[1] /= 2.0;
            tmp_pack_liq[1] /= 2.0;
            tmp_depth[2] = tmp_depth[1];
            tmp_pack_ice[2] = tmp_pack_ice[1];
            tmp_pack_liq[2] = tmp_pack_liq[1];
            tmp_pack_T[2] = tmp_pack_T[1] - grad_temp * tmp_depth[1] / 2.0;
            if (tmp_pack_T[2] >= CONST_TKFRZ) {
                tmp_pack_T[2] = tmp_pack_T[1];
            }
            else {
                tmp_pack_T[1] = tmp_pack_T[1] + grad_temp * tmp_depth[1] / 2.0;
            }
        }
    }

    if (tmp_Nsnow > 2) {
        if (tmp_depth[1] > 0.20) {
            total_depth = tmp_depth[1] - 0.20;
            extra_frac = total_depth / tmp_depth[1];
            extra_ice = extra_frac * tmp_pack_ice[1];
            extra_liq = extra_frac * tmp_pack_liq[1];
            extra_frac = 0.20 / tmp_depth[1];
            tmp_pack_ice[1] *= extra_frac;
            tmp_pack_liq[1] *= extra_frac;
            tmp_depth[1] = 0.20;

            /* update combined snow water & temperature */
            update_surface_fluxes(&tmp_depth[2], 
                                  &tmp_pack_liq[2],
                                  &tmp_pack_ice[2],
                                  &tmp_pack_T[2],
                                  total_depth, 
                                  extra_liq,
                                  extra_ice,
                                  tmp_pack_T[1]);
        }
    }
    snow->Nsnow = tmp_Nsnow;
    for (i = 0; i < snow->Nsnow; i++) {
        snow->dz_snow[i] = tmp_depth[i];
        snow->pack_ice[i] = tmp_pack_ice[i];
        snow->pack_liq[i] = tmp_pack_liq[i];
        snow->pack_T[i] = tmp_pack_T[i];
    }
}