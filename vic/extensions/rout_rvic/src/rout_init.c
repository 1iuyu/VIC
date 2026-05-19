/******************************************************************************
 * @section DESCRIPTION
 *
 * Initialize routing model parameters
 *****************************************************************************/

#include "rout.h"

/******************************************************************************
 * @brief    Initialize routing model parameters
 *****************************************************************************/
void
rout_init(void)
{
    extern int              mpi_rank;
    extern rout_struct      rout;
    extern domain_struct    global_domain;
    extern domain_struct    local_domain;
    extern filenames_struct filenames;
    int                     status;

    if (mpi_rank == VIC_MPI_ROOT) {
        int    *direction = NULL;
        double *dvar = NULL;
        int    *ivar = NULL;
        size_t  i, j, d, idx;
        int     ii, jj;
        int     n_nx, n_ny;
        size_t  d2count[2];
        size_t  d2start[2];
        d2start[0] = 0;
        d2start[1] = 0;
        d2count[0] = global_domain.n_ny;
        d2count[1] = global_domain.n_nx;

        n_nx = global_domain.n_nx;
        n_ny = global_domain.n_ny;

        // allocate memory for variables to be read
        direction = malloc(global_domain.ncells_total * sizeof(*direction));
        check_alloc_status(direction, "Memory allocation error.");
        dvar = malloc(global_domain.ncells_total * sizeof(*dvar));
        check_alloc_status(dvar, "Memory allocation error.");
        ivar = malloc(global_domain.ncells_total * sizeof(*ivar));
        check_alloc_status(ivar, "Memory allocation error.");
        
        // open parameter file
        status = nc_open(filenames.rout_params.nc_filename, NC_NOWRITE,
                         &(filenames.rout_params.nc_id));
        check_nc_status(status, "Error opening %s",
                        filenames.rout_params.nc_filename);

        get_nc_field_int(&(filenames.rout_params), "direction", d2start,
                                   d2count, ivar);
        for (i = 0; i < global_domain.ncells_total; i++) {
            direction[i] = (int) ivar[i];
        }
        // 验证flow direction的有效性
        for (i = 0; i < global_domain.ncells_total; i++) {
            if (direction[i] < 0 || direction[i] > 8) {
                log_err("Invalid flow direction value %d at grid cell %zu. "
                           "Flow direction should be encoded as D8 (0-8).", 
                                        direction[i], i);
            }
        }
        // 根据flow direction构建source_torow和source_tocol数组
        // D8编码：0-无流向，1-北，2-东北，3-东，4-东南，5-南，6-西南，7-西，8-西北
        SearchCatchment(direction);
        // 下游索引
        for (i = 0; i < local_domain.ncells_active; i++) {
            // 格点的全局索引
            ii = rout.rout_param.source_torow[i];
            jj = rout.rout_param.source_tocol[i];
            idx = ii * n_nx + jj;
            if (ii>=0 && ii<n_ny && jj>=0 && jj<n_nx) {
                rout.rout_param.downstream[i] = global_domain.locations[idx].global_idx;  // 下游索引
            }
            else {
                rout.rout_param.downstream[i] = 99999;   // 出口或边界
            }
        }
        // 计算入度
        for (i = 0; i < local_domain.ncells_active; i++) {
            rout.rout_param.indegree[i] = 0;  // 初始化入度为0
        }
        for (j = 0; j < local_domain.ncells_active; j++) {
            d = rout.rout_param.downstream[j];  // 当前单元格的下游索引
            if (d < local_domain.ncells_active) {
                rout.rout_param.indegree[d]++;  // 下游单元格的入度增加
            }
        }

        // 初始化队列（入度为0的节点）
        int queue_front = 0, queue_rear = 0;
        for (j = 0; j < local_domain.ncells_active; j++) {
            if (rout.rout_param.indegree[j] == 0) {
                rout.rout_param.queue[queue_rear++] = j;  // 入度为0的单元格加入队列
            }
        }

        // 拓扑排序
        size_t count = 0;
        while (queue_front < queue_rear) {
            j = rout.rout_param.queue[queue_front++];  // 出队
            rout.rout_param.routing_order[count++] = j;  // 记录排序结果
            
            d = rout.rout_param.downstream[j];
            if (d < local_domain.ncells_active) {
                rout.rout_param.indegree[d]--;  // 减少下游的入度
                if (rout.rout_param.indegree[d] == 0) {
                    rout.rout_param.queue[queue_rear++] = d;  // 如果入度变为0，加入队列
                }
            }
        }
        if (count != local_domain.ncells_active) {
            log_err("Error: The flow direction network contains a cycle.");
        }

        // 初始化：每个格点的汇水面积 = 自身面积
        for (i = 0; i < local_domain.ncells_active; i++) {
            rout.acc_area[i] = local_domain.locations[i].area;
        }

        // 按拓扑序计算（从上游到下游）
        for (j = 0; j < local_domain.ncells_active; j++) {
            i = rout.rout_param.routing_order[j];   // 当前格点
            d = rout.rout_param.downstream[i];      // 下游格点

            if (d < local_domain.ncells_active) {
                rout.acc_area[d] += rout.acc_area[i];  // 把自己的汇水面积加到下游
            }
        }
        for (i = 0; i < local_domain.ncells_active; i++) {
            rout.main_channel.rwidth[i] = 0.001 * pow(rout.acc_area[i], 0.5);
            rout.main_channel.rwidth0[i] = 5.0 * rout.main_channel.rwidth[i];
            rout.main_channel.rdepth[i] = pow(rout.main_channel.rwidth[i], 0.3333);
        }

        // channel length: river length (m)
        get_nc_field_double(&(filenames.rout_params),
                            "rlen", d2start, d2count, dvar);
        for (i = 0, j = 0; i < global_domain.ncells_total; i++) {
            if (global_domain.locations[i].run) {
                rout.main_channel.rlen[j] = (double) dvar[i];
                j++;
            }
        }
        // sub_channel length: (m)
        get_nc_field_double(&(filenames.rout_params),
                            "tlen", d2start, d2count, dvar);
        for (i = 0, j = 0; i < global_domain.ncells_total; i++) {
            if (global_domain.locations[i].run) {
                rout.sub_channel.tlen[j] = (double) dvar[i];
                j++;
            }
        }
        // drainage density: rainage density within the cell, [1/m]
        for (i = 0; i < local_domain.ncells_active; i++) {
            rout.total_length[i] = rout.main_channel.rlen[i] + rout.sub_channel.tlen[i];
            rout.drainage_density[i] = rout.total_length[i] / local_domain.locations[i].area;
        }
        // rslpsqrt: sqrt channel slope
        get_nc_field_double(&(filenames.rout_params),
                            "rslpsqrt", d2start, d2count, dvar);
        for (i = 0, j = 0; i < global_domain.ncells_total; i++) {
            if (global_domain.locations[i].run) {
                rout.main_channel.rslpsqrt[j] = (double) dvar[i];
                j++;
            }
        }
        // tslpsqrt: sqrt sub-channel slope
        get_nc_field_double(&(filenames.rout_params),
                            "tslpsqrt", d2start, d2count, dvar);
        for (i = 0, j = 0; i < global_domain.ncells_total; i++) {
            if (global_domain.locations[i].run) {
                rout.sub_channel.tslpsqrt[j] = (double) dvar[i];
                j++;
            }
        }
        // hslpsqrt: sqrt hillslope slope
        get_nc_field_double(&(filenames.rout_params),
                            "hslpsqrt", d2start, d2count, dvar);
        for (i = 0, j = 0; i < global_domain.ncells_total; i++) {
            if (global_domain.locations[i].run) {
                rout.hillslope.hslpsqrt[j] = (double) dvar[i];
                j++;
            }
        }

        // 主河道和子河道长度验证
        double hlen_max = 0.0;
        double rlen_min = 0.0;
        double river_depth = 0.0;
        double hydrR = 0.0;
        double velocity = 0.0;
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (rout.main_channel.rlen[i] > 0.0) {
                rout.hillslope.hlen[i] = local_domain.locations[i].area /
                                                    rout.total_length[i] / 2.0;
                hlen_max = max(1000.0, sqrt(local_domain.locations[i].area));
                if (rout.hillslope.hlen[i] > hlen_max) {
                    rout.hillslope.hlen[i] = hlen_max;
                }
                rlen_min = sqrt(local_domain.locations[i].area);
                if (rout.main_channel.rlen[i] < rlen_min) {
                    rout.main_channel.rlen[i] = rlen_min;
                }
                if (rout.sub_channel.twidth[i] < 0.0) {
                    rout.sub_channel.twidth[i] = 0.0;
                }
                if (rout.sub_channel.tlen[i] > 0.0) {
                    rout.sub_channel.twidth[i] = 0.001 * pow(local_domain.locations[i].area, 0.5) * 0.6;
                    if (rout.sub_channel.twidth[i] < 0.0) {
                        rout.sub_channel.twidth[i] = 0.0;
                    }
                }
            }
            else {
                rout.hillslope.hlen[i] = 0.0;
                rout.sub_channel.tlen[i] = 0.0;
                rout.sub_channel.twidth[i] = 0.0;
            }
            rout.hillslope.nh[i] = 0.4;
            rout.hillslope.yh[i] = 0.0; // 坡面水深初始化为0.001m
            rout.hillslope.wh[i] = 0.0; // 坡面水量初始化为0.001m
            rout.sub_channel.nt[i] = 0.05;
            rout.sub_channel.wt[i] = rout.sub_channel.tlen[i] * rout.sub_channel.twidth[i] * 0.5;  // 初始化子河道水量
            if (rout.sub_channel.tlen[i] > 0.0 && rout.sub_channel.wt[i] > 0.0) {
                rout.sub_channel.mt[i] = GRMR(rout.sub_channel.wt[i], rout.sub_channel.tlen[i]);   // 过水面积（单位长度）
                rout.sub_channel.yt[i] = GRHT(rout.sub_channel.mt[i], rout.sub_channel.twidth[i]); // 水深
                rout.sub_channel.pt[i] = GRPT(rout.sub_channel.yt[i], rout.sub_channel.twidth[i]); // 湿周
                rout.sub_channel.rt[i] = GRRR(rout.sub_channel.mt[i], rout.sub_channel.pt[i]);     // 水力半径
            }
            else {
                rout.sub_channel.mt[i] = 0.0;
                rout.sub_channel.yt[i] = 0.0;
                rout.sub_channel.pt[i] = 0.0;
                rout.sub_channel.rt[i] = 0.0;
            }
            rout.main_channel.nr[i] = 0.05;
            river_depth = rout.main_channel.rdepth[i] * 0.5;
            rout.main_channel.wr[i] = rout.main_channel.rlen[i] * rout.main_channel.rwidth[i] * river_depth; // 初始化主河道水量
            hydrR = rout.main_channel.rwidth[i] * river_depth / (rout.main_channel.rwidth[i] + 2.0 * river_depth);
            velocity = CRVRMAN_nosqrt(rout.main_channel.rslpsqrt[i], rout.main_channel.nr[i], hydrR);
            rout.main_channel.erout[i] = -velocity * rout.main_channel.rwidth[i] * river_depth; // 计算主河道出流量
            rout.total_storage_prev += rout.hillslope.wh[i] * local_domain.locations[i].area * 
                                        local_domain.locations[i].frac + rout.sub_channel.wt[i] + rout.main_channel.wr[i];
            // 更新主河道水力参数
            if(rout.main_channel.rlen[i] > 0.0 && rout.main_channel.wr[i] > 0.0) {
                rout.main_channel.mr[i] = GRMR(rout.main_channel.wr[i], rout.main_channel.rlen[i]);
                rout.main_channel.yr[i] = GRHR(rout.main_channel.mr[i], rout.main_channel.rwidth[i], 
                                        rout.main_channel.rwidth0[i], rout.main_channel.rdepth[i]);
                rout.main_channel.pr[i] = GRPR(rout.main_channel.yr[i], rout.main_channel.rwidth[i],
                                        rout.main_channel.rwidth0[i], rout.main_channel.rdepth[i]);
                rout.main_channel.rr[i] = GRRR(rout.main_channel.mr[i], rout.main_channel.pr[i]);
            }
            else {
                rout.main_channel.mr[i] = 0.0;
                rout.main_channel.yr[i] = 0.0;
                rout.main_channel.pr[i] = 0.0;
                rout.main_channel.rr[i] = 0.0;
            }
        }

        // 计算sub-channel和main-channel的汇流时间
        double river_T = 0.0;
        double sub_T = 0.0;
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (rout.main_channel.rlen[i] > 0.0) {
                river_T = rout.acc_area[i] * rout.main_channel.rslpsqrt[i] / 
                                    (rout.main_channel.rlen[i] * rout.main_channel.rwidth[i]);
                if (river_T >= 10.0) {
                    rout.river_steps[i] = log10(river_T) * DLevelH2R + 1;
                }
                else {
                    rout.river_steps[i] = DLevelH2R + 1;
                }
            }
            if (rout.river_steps[i] < 1) {
                rout.river_steps[i] = 1;
            }
            if (rout.sub_channel.tlen[i] > 0.0) {
                sub_T = local_domain.locations[i].area * rout.sub_channel.tslpsqrt[i] / 
                                    (rout.sub_channel.tlen[i] * rout.sub_channel.twidth[i]);
                if (sub_T >= 10.0) {
                    rout.sub_steps[i] = log10(sub_T) * DLevelH2R + 1;
                }
                else {
                    rout.sub_steps[i] = DLevelH2R + 1;
                }
            }
            if (rout.sub_steps[i] < 1) {
                rout.sub_steps[i] = 1;
            }
        }

        // close parameter file
        status = nc_close(filenames.rout_params.nc_id);
        check_nc_status(status, "Error closing %s",
                        filenames.rout_params.nc_filename);

        // cleanup
        free(dvar);
        free(ivar);
        free(direction);
    }
}
