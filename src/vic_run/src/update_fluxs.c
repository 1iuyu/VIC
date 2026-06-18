/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes the fluxes for each cell in the grid. It sets the
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    This routine initializes roughness lengths and 
 *           forcing heights for each cell in the grid.
 *****************************************************************************/
void
update_fluxes(energy_bal_struct *energy,
              cell_data_struct  *cell,
              snow_data_struct  *snow,
              soil_con_struct   *soil_con)
{
    extern option_struct     options;

    /********************************
       初始化地表温度为各部分的面积加权
    ********************************/
    double coverage = snow->coverage;
    double frac_h2o = cell->frac_h2o;
    double *soil_T = cell->soil_T;
    double *pack_T = snow->pack_T;
    double h2osfc_T = cell->h2osfc_T;
    if (snow->Nsnow > 0) {
        if (frac_h2o > 0.0) {
            energy->Tgrnd = (1.0 - coverage) * soil_T[0] + 
                                    frac_h2o * h2osfc_T + coverage * pack_T[0];
        }
        else {
            energy->Tgrnd = (1.0 - coverage) * soil_T[0] + coverage * pack_T[0];
        }
    }
    else {
        if (frac_h2o > 0.0) {
            energy->Tgrnd = (1.0 - frac_h2o) * soil_T[0] + frac_h2o * h2osfc_T;
        }
        else {
            energy->Tgrnd = soil_T[0];
        }
    }

    /********************************
     计算活动层深度并更新根区分数
    ********************************/
    if (options.ACTIVE_LAYER) {
        ActiveLayer(cell, soil_con);
    }    
}