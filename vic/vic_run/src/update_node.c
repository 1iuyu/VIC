/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine update the depth of layer nodes in the vertical direction.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    This routine update the depth of layer nodes in the vertical
 *           direction.
 *****************************************************************************/
void
update_node(cell_data_struct  *cell,
            snow_data_struct  *snow,
            soil_con_struct   *soil_con)
{
    size_t i, k, lidx;
    size_t Nsnow = snow->Nsnow;
    size_t Nsoil = cell->Nsoil;
    size_t Nnode = Nsnow + Nsoil;
    double *dz_node = cell->dz_node;
    double *dz_snow = snow->dz_snow;
    double *zc_node = cell->zc_node;
    double *zc_snow = snow->zc_snow;
    double *Zsum_node = cell->Zsum_node;
    double *Zsum_snow = snow->Zsum_snow;
    double *dz_soil = soil_con->dz_soil;
    /********************************
      更新结构体中的节点位置和界面位置
    ********************************/
    cell->Nnode = Nnode;
    // 温度节点[0 到 Nsnow-1]
    for (i = 0; i < Nsnow; i++) {
        dz_node[i] = dz_snow[i];
    }
    // 土壤节点[Nsnow 到 Nsnow+Nsoil-1]
    for (i = 0; i < Nsoil; i++) {
        lidx = Nsnow + i;
        dz_node[lidx] = dz_soil[i];
    }

    /* 计算节点总深度和节点总高度 */
    if (Nsnow > 0) {
        for (k = 0; k < Nsnow; k++) {
            lidx = Nsnow - 1 - k;
            Zsum_node[k] = Zsum_snow[lidx];
        }
    }
    Zsum_node[Nsnow] = 0.0;
    double sum_dz = 0.;
    for (k = Nsnow; k < Nnode; k++) {
        sum_dz += dz_node[k];
        Zsum_node[k+1] = sum_dz;
    }
    // 节点中心坐标
    if (Nsnow > 0) {
        for (k = 0; k < Nsnow; k++) {
            lidx = Nsnow - 1 - k;
            zc_node[k] = zc_snow[lidx];
        }
    }
    for (k = Nsnow; k < Nnode; k++) {
        zc_node[k] = Zsum_node[k+1] - dz_node[k] / 2.0;
    }
}