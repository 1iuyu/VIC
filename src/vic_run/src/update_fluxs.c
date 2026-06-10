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
    if (snow->Nsnow > 0) {
        energy->Tgrnd = (1.0 - coverage) * 
                cell->soil_T[0] + coverage * snow->pack_T[0];
    }
    else {
        energy->Tgrnd = cell->soil_T[0];
    }

    if (cell->IS_VEG) {

        /********************************
         计算活动层深度并更新根区分数
        ********************************/
        if (options.ACTIVE_LAYER) {
            ActiveLayer(cell, soil_con);
        }
    }
    
}