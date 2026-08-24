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
    extern rout_struct     *rout;
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
            ii = rout[i].rout_param.source_torow;
            jj = rout[i].rout_param.source_tocol;
            idx = ii * n_nx + jj;
            if (ii>=0 && ii<n_ny && jj>=0 && jj<n_nx) {
                rout[i].rout_param.downstream = global_domain.locations[idx].global_idx;  // 下游索引
            }
            else {
                rout[i].rout_param.downstream = 99999;   // 出口或边界
            }
        }
        // 计算入度
        for (i = 0; i < local_domain.ncells_active; i++) {
            rout[i].rout_param.indegree = 0;  // 初始化入度为0
        }
        for (j = 0; j < local_domain.ncells_active; j++) {
            d = rout[j].rout_param.downstream;  // 当前单元格的下游索引
            if (d < local_domain.ncells_active) {
                rout[d].rout_param.indegree++;  // 下游单元格的入度增加
            }
        }

        // 初始化队列（入度为0的节点）
        int queue_front = 0, queue_rear = 0;
        for (j = 0; j < local_domain.ncells_active; j++) {
            if (rout[j].rout_param.indegree == 0) {
                rout[queue_rear++].rout_param.queue = j;  // 入度为0的单元格加入队列
            }
        }

        // 拓扑排序
        size_t count = 0;
        while (queue_front < queue_rear) {
            j = rout[queue_front++].rout_param.queue;  // 出队
            rout[count++].rout_param.routing_order = j;  // 记录排序结果
            
            d = rout[j].rout_param.downstream;
            if (d < local_domain.ncells_active) {
                rout[d].rout_param.indegree--;  // 减少下游的入度
                if (rout[d].rout_param.indegree == 0) {
                    rout[queue_rear++].rout_param.queue = d;  // 如果入度变为0，加入队列
                }
            }
        }
        if (count != local_domain.ncells_active) {
            log_err("Error: The flow direction network contains a cycle.");
        }

        // 初始化：每个格点的汇水面积 = 自身面积
        for (i = 0; i < local_domain.ncells_active; i++) {
            rout[i].acc_area = local_domain.locations[i].area;
        }

        // 按拓扑序计算（从上游到下游）
        for (j = 0; j < local_domain.ncells_active; j++) {
            i = rout[j].rout_param.routing_order;   // 当前格点
            d = rout[i].rout_param.downstream;      // 下游格点

            if (d < local_domain.ncells_active) {
                rout[d].acc_area += rout[i].acc_area;  // 把自己的汇水面积加到下游
            }
        }
        for (i = 0; i < local_domain.ncells_active; i++) {
            rout[i].main_channel.rwidth = 0.001 * pow(rout[i].acc_area, 0.5);
            rout[i].main_channel.rwidth0 = 5.0 * rout[i].main_channel.rwidth;
            rout[i].main_channel.rdepth = pow(rout[i].main_channel.rwidth, 0.3333);
        }

        // channel length: river length (m)
        get_nc_field_double(&(filenames.rout_params),
                            "rlen", d2start, d2count, dvar);
        for (i = 0, j = 0; i < global_domain.ncells_total; i++) {
            if (global_domain.locations[i].run) {
                rout[j].main_channel.rlen = (double) dvar[i];
                j++;
            }
        }
        // sub_channel length: (m)
        get_nc_field_double(&(filenames.rout_params),
                            "tlen", d2start, d2count, dvar);
        for (i = 0, j = 0; i < global_domain.ncells_total; i++) {
            if (global_domain.locations[i].run) {
                rout[j].sub_channel.tlen = (double) dvar[i];
                j++;
            }
        }
        // drainage density: rainage density within the cell, [1/m]
        for (i = 0; i < local_domain.ncells_active; i++) {
            rout[i].total_length = rout[i].main_channel.rlen + rout[i].sub_channel.tlen;
            rout[i].drainage_density = rout[i].total_length / local_domain.locations[i].area;
        }
        // rslpsqrt: sqrt channel slope
        get_nc_field_double(&(filenames.rout_params),
                            "rslpsqrt", d2start, d2count, dvar);
        for (i = 0, j = 0; i < global_domain.ncells_total; i++) {
            if (global_domain.locations[i].run) {
                rout[j].main_channel.rslpsqrt = (double) dvar[i];
                j++;
            }
        }
        // tslpsqrt: sqrt sub-channel slope
        get_nc_field_double(&(filenames.rout_params),
                            "tslpsqrt", d2start, d2count, dvar);
        for (i = 0, j = 0; i < global_domain.ncells_total; i++) {
            if (global_domain.locations[i].run) {
                rout[j].sub_channel.tslpsqrt = (double) dvar[i];
                j++;
            }
        }
        // hslpsqrt: sqrt hillslope slope
        get_nc_field_double(&(filenames.rout_params),
                            "hslpsqrt", d2start, d2count, dvar);
        for (i = 0, j = 0; i < global_domain.ncells_total; i++) {
            if (global_domain.locations[i].run) {
                rout[j].hillslope.hslpsqrt = (double) dvar[i];
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
            if (rout[i].main_channel.rlen > 0.0) {
                rout[i].hillslope.hlen = local_domain.locations[i].area /
                                                    rout[i].total_length / 2.0;
                hlen_max = max(1000.0, sqrt(local_domain.locations[i].area));
                if (rout[i].hillslope.hlen > hlen_max) {
                    rout[i].hillslope.hlen = hlen_max;
                }
                rlen_min = sqrt(local_domain.locations[i].area);
                if (rout[i].main_channel.rlen < rlen_min) {
                    rout[i].main_channel.rlen = rlen_min;
                }
                if (rout[i].sub_channel.twidth < 0.0) {
                    rout[i].sub_channel.twidth = 0.0;
                }
                if (rout[i].sub_channel.tlen > 0.0) {
                    rout[i].sub_channel.twidth = 0.001 * pow(local_domain.locations[i].area, 0.5) * 0.6;
                    if (rout[i].sub_channel.twidth < 0.0) {
                        rout[i].sub_channel.twidth = 0.0;
                    }
                }
            }
            else {
                rout[i].hillslope.hlen = 0.0;
                rout[i].sub_channel.tlen = 0.0;
                rout[i].sub_channel.twidth = 0.0;
            }
            rout[i].hillslope.nh = 0.4;
            rout[i].hillslope.yh = 0.0; // 坡面水深初始化为0.001m
            rout[i].hillslope.wh = 0.0; // 坡面水量初始化为0.001m
            rout[i].sub_channel.nt = 0.05;
            rout[i].sub_channel.wt = rout[i].sub_channel.tlen * rout[i].sub_channel.twidth * 0.5;  // 初始化子河道水量
            if (rout[i].sub_channel.tlen > 0.0 && rout[i].sub_channel.wt > 0.0) {
                rout[i].sub_channel.mt = GRMR(rout[i].sub_channel.wt, rout[i].sub_channel.tlen);   // 过水面积（单位长度）
                rout[i].sub_channel.yt = GRHT(rout[i].sub_channel.mt, rout[i].sub_channel.twidth); // 水深
                rout[i].sub_channel.pt = GRPT(rout[i].sub_channel.yt, rout[i].sub_channel.twidth); // 湿周
                rout[i].sub_channel.rt = GRRR(rout[i].sub_channel.mt, rout[i].sub_channel.pt);     // 水力半径
            }
            else {
                rout[i].sub_channel.mt = 0.0;
                rout[i].sub_channel.yt = 0.0;
                rout[i].sub_channel.pt = 0.0;
                rout[i].sub_channel.rt = 0.0;
            }
            rout[i].main_channel.nr = 0.05;
            river_depth = rout[i].main_channel.rdepth * 0.5;
            rout[i].main_channel.wr = rout[i].main_channel.rlen * rout[i].main_channel.rwidth * river_depth; // 初始化主河道水量
            hydrR = rout[i].main_channel.rwidth * river_depth / (rout[i].main_channel.rwidth + 2.0 * river_depth);
            velocity = CRVRMAN_nosqrt(rout[i].main_channel.rslpsqrt, rout[i].main_channel.nr, hydrR);
            rout[i].main_channel.erout = -velocity * rout[i].main_channel.rwidth * river_depth; // 计算主河道出流量
            rout[i].total_storage_prev += rout[i].hillslope.wh * local_domain.locations[i].area * 
                                        local_domain.locations[i].frac + rout[i].sub_channel.wt + rout[i].main_channel.wr;
            // 更新主河道水力参数
            if(rout[i].main_channel.rlen > 0.0 && rout[i].main_channel.wr > 0.0) {
                rout[i].main_channel.mr = GRMR(rout[i].main_channel.wr, rout[i].main_channel.rlen);
                rout[i].main_channel.yr = GRHR(rout[i].main_channel.mr, rout[i].main_channel.rwidth, 
                                        rout[i].main_channel.rwidth0, rout[i].main_channel.rdepth);
                rout[i].main_channel.pr = GRPR(rout[i].main_channel.yr, rout[i].main_channel.rwidth,
                                        rout[i].main_channel.rwidth0, rout[i].main_channel.rdepth);
                rout[i].main_channel.rr = GRRR(rout[i].main_channel.mr, rout[i].main_channel.pr);
            }
            else {
                rout[i].main_channel.mr = 0.0;
                rout[i].main_channel.yr = 0.0;
                rout[i].main_channel.pr = 0.0;
                rout[i].main_channel.rr = 0.0;
            }
        }

        // 计算sub-channel和main-channel的汇流时间
        double river_T = 0.0;
        double sub_T = 0.0;
        for (i = 0; i < local_domain.ncells_active; i++) {
            if (rout[i].main_channel.rlen > 0.0) {
                river_T = rout[i].acc_area * rout[i].main_channel.rslpsqrt / 
                                    (rout[i].main_channel.rlen * rout[i].main_channel.rwidth);
                if (river_T >= 10.0) {
                    rout[i].river_steps = log10(river_T) * DLevelH2R + 1;
                }
                else {
                    rout[i].river_steps = DLevelH2R + 1;
                }
            }
            if (rout[i].river_steps < 1) {
                rout[i].river_steps = 1;
            }
            if (rout[i].sub_channel.tlen > 0.0) {
                sub_T = local_domain.locations[i].area * rout[i].sub_channel.tslpsqrt / 
                                    (rout[i].sub_channel.tlen * rout[i].sub_channel.twidth);
                if (sub_T >= 10.0) {
                    rout[i].sub_steps = log10(sub_T) * DLevelH2R + 1;
                }
                else {
                    rout[i].sub_steps = DLevelH2R + 1;
                }
            }
            if (rout[i].sub_steps < 1) {
                rout[i].sub_steps = 1;
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
