/******************************************************************************
 * @section DESCRIPTION
 *
 * Run routing over the domain.
 *****************************************************************************/

#include "rout.h"

/******************************************************************************
* @brief        This subroutine controls the MOSART routing over the domain.
******************************************************************************/
void
rout_run(void)
{
    extern double     ***out_data;
    extern domain_struct local_domain;
    extern rout_struct   rout;
    extern force_data_struct  *force;
    extern parameters_struct   param;
    extern global_param_struct global_param;
    int    ErrorFlag;
    size_t i, j, d;
    size_t step_dt = global_param.step_dt;
    
    double total_storage = 0.0;     // 当前总水量
    double total_input = 0.0;       // 输入（水文产生）
    double total_output = 0.0;      // 流出（出口）
    // Read from runoff and baseflow from out_data and sum to runoff
    for (i = 0; i < local_domain.ncells_active; i++) {
        rout.runoff[i] = out_data[i][OUT_RUNOFF][0];
        rout.baseflow[i] = out_data[i][OUT_BASEFLOW][0];
        if (rout.baseflow[i] < 0.0) {
            rout.baseflow[i] = 0.0;
        }
    }
    // 初始化upstream
    for (i = 0; i < local_domain.ncells_active; i++) {
        rout.upstream[i] = 0.0;
        rout.discharge[i] = 0.0;
    }
    // 存在模拟区域外的上游来水，需要将其加入到upstream中
    for (i = 0; i < local_domain.ncells_active; i++) {
        rout.upstream[i] = force[i].channel_in[NR];
    }
    // Run the convolution on the master node
    for (i = 0; i < local_domain.ncells_active; i++) {
        // 格点汇流顺序
        j = rout.rout_param.routing_order[i];
        // 下游格点索引
        d = rout.rout_param.downstream[j];

        ErrorFlag = Euler_Routing(j, step_dt);

        total_input += rout.sub_channel.etin[j];

        if (ErrorFlag) {
            fprintf(stderr, "Error in Euler_Routing\n");
            exit(EXIT_FAILURE);
        }
        // 将当前格点的流量加到下游格点的流量上
        if (d < local_domain.ncells_active) {
            rout.upstream[d] += rout.discharge[j];
        }
        else {
            // ====== 流域出口 ======
            total_output += rout.discharge[j];
        }
    }
    // 统计总存储（所有水体）
    for (i = 0; i < local_domain.ncells_active; i++) {
        // 计算总存储[m/s]
        total_storage +=
            (rout.hillslope.wh[i] * 
            local_domain.locations[i].area * 
            local_domain.locations[i].frac) +      // 坡面水
            rout.sub_channel.wt[i] +    // 子河道
            rout.main_channel.wr[i];    // 主河道
    }
    // 计算存储变化
    double storage_change = total_storage - rout.total_storage_prev;

    double balance_error =
        (total_input - total_output) * step_dt - storage_change;
    
    // 更新总存储
    rout.total_storage_prev = total_storage;

    if (fabs(balance_error) > param.TOL_A) {
        log_warn("Water balance error too large: %.3e", balance_error);
    }
    
    // Write to output struct
    for (i = 0; i < local_domain.ncells_active; i++) {
        out_data[i][OUT_DISCHARGE][0] = rout.discharge[i];
    }
}
