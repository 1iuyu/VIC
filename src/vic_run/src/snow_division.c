/******************************************************************************
* @section DESCRIPTION
*
* Snowpack layer division process
* Update snow ice, snow water, snow thickness, snow temperature.
******************************************************************************/

#include "vic_run.h"

/******************************************************************************
* @brief    Calculate snow depth via compaction due to destructive metamorphism, 
*           overburden, & melt.
******************************************************************************/
void 
snow_division(snow_data_struct *snow)

{              
    size_t i, tmp_Nsnow;
    double grad_temp;
    double extra_frac;
    double tmp_pack_ice[MAX_SNOWS] = {0};
    double tmp_pack_liq[MAX_SNOWS] = {0};
    double tmp_pack_T[MAX_SNOWS] = {0};
    double tmp_depth[MAX_SNOWS] = {0};
    double *dz_snow = snow->dz_snow;
    double *pack_T = snow->pack_T;
    double *pack_liq = snow->pack_liq;
    double *pack_ice = snow->pack_ice;
    double total_depth = 0.0;
    double extra_ice = 0.0;
    double extra_liq = 0.0;

    for (i = 0; i < snow->Nsnow; i++) {
        tmp_depth[i] = dz_snow[i];
        tmp_pack_ice[i] = pack_ice[i];
        tmp_pack_liq[i] = pack_liq[i];
        tmp_pack_T[i] = pack_T[i];
    }
    tmp_Nsnow = snow->Nsnow;

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
        update_snow_fluxes(&tmp_depth[1], 
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
            update_snow_fluxes(&tmp_depth[2], 
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
        dz_snow[i] = tmp_depth[i];
        pack_ice[i] = tmp_pack_ice[i];
        pack_liq[i] = tmp_pack_liq[i];
        pack_T[i] = tmp_pack_T[i];
    }
}