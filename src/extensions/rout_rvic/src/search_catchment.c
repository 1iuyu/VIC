/******************************************************************************
 * @section DESCRIPTION
 *
 * SearchCatchment - Find number of cells upstream current gage location, and their
 * row and col number
 *****************************************************************************/

#include "rout.h"

/*****************************************************************************
 * @brief  Search Catchment number
*****************************************************************************/
void SearchCatchment(int *direction)
{
    extern domain_struct global_domain;
    extern rout_struct  *rout;
    size_t i, j;
    size_t n_nx;
    n_nx = global_domain.n_nx;

    // 计算 downstream（torow, tocol）
    for (i = 0, j = 0; i < global_domain.ncells_total; i++) {
        size_t x_index = i % n_nx;
        size_t y_index = i / n_nx;
        if (direction[i] != 0 && global_domain.locations[i].run) {
            switch (direction[i]) {
            case 1:
                rout[j].rout_param.source_torow = y_index - 1;
                rout[j].rout_param.source_tocol = x_index;
            break;
            case 2:
                rout[j].rout_param.source_torow = y_index + 1;
                rout[j].rout_param.source_tocol = x_index + 1;
            break;
            case 3:
                rout[j].rout_param.source_torow = y_index;
                rout[j].rout_param.source_tocol = x_index + 1;
            break;
            case 4:
                rout[j].rout_param.source_torow = y_index + 1;
                rout[j].rout_param.source_tocol = x_index + 1;
            break;
            case 5:
                rout[j].rout_param.source_torow = y_index + 1;
                rout[j].rout_param.source_tocol = x_index;
            break;
            case 6:
                rout[j].rout_param.source_torow = y_index + 1;
                rout[j].rout_param.source_tocol = x_index - 1;
            break;
            case 7:
                rout[j].rout_param.source_torow = y_index;
                rout[j].rout_param.source_tocol = x_index - 1;
            break;
            case 8:
                rout[j].rout_param.source_torow = y_index - 1;
                rout[j].rout_param.source_tocol = x_index - 1;
            break;
            default:
                rout[j].rout_param.source_torow = 0;
                rout[j].rout_param.source_tocol = 0;
            }
            j++;
        }
    }
}