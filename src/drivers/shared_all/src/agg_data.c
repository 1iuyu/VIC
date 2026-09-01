/******************************************************************************
 * @section DESCRIPTION
 *
 * Initialize output structures.
 *****************************************************************************/

#include "vic_driver_shared_all.h"

/******************************************************************************
 * @brief    Perform temporal aggregation on stream data
 *****************************************************************************/
void
agg_stream_data(stream_struct     *stream,
                dmy_struct        *dmy_current,
                veg_con_struct   **veg_con,
                double         ****out_data)
{
    alarm_struct          *alarm;
    size_t                 i, j, h, k;
    size_t                 nelem;
    size_t                 nveg;
    size_t                 max_nelem;
    unsigned int           varid;
    bool                   alarm_now;
    bool                   is_first = true;
    double cell_value[MAX_NODES] = {0};
    alarm = &(stream->agg_alarm);
    alarm->count++;
    alarm_now = raise_alarm(alarm, dmy_current);

    if (alarm->count == 1) {
        stream->time_bounds[0] = (*dmy_current);
    }

    if (alarm_now) {
        stream->time_bounds[1] = (*dmy_current);
    }

    for (i = 0; i < stream->ngridcells; i++) {
        nveg = stream->nveg[i];
        for (j = 0; j < stream->nvars; j++) {
            varid = stream->varid[j];
            // 提前确定聚合模式
            bool is_end = (stream->aggtype[j] == AGG_TYPE_END) && alarm_now;
            bool is_beg = (stream->aggtype[j] == AGG_TYPE_BEG) && (alarm->count == 1);
            bool is_sum = (stream->aggtype[j] == AGG_TYPE_SUM);
            bool is_avg = (stream->aggtype[j] == AGG_TYPE_AVG);
            bool is_max = (stream->aggtype[j] == AGG_TYPE_MAX);
            bool is_min = (stream->aggtype[j] == AGG_TYPE_MIN);
            bool apply_avg_div = is_avg && alarm_now;
            /* * CELL variables have only one output value per cell. 
               * HRU variables have one output value for each HRU. */ 
            if (stream->domain[j] == OUT_DOMAIN_HRU) {

                for (h = 0; h < nveg; h++) {
                    nelem = get_output_nelem(varid, i, h);
                    // Instantaneous at the end of the period
                    if (is_end) {
                        for (k = 0; k < nelem; k++) {
                            stream->aggdata[i][j][h][k][0] = out_data[i][varid][h][k];
                        }
                    }
                    // Instantaneous at the beginning of the period
                    else if (is_beg) {
                        for (k = 0; k < nelem; k++) {
                            stream->aggdata[i][j][h][k][0] = out_data[i][varid][h][k];
                        }
                    }
                    // Sum over the period
                    else if (is_sum || is_avg) {
                        for (k = 0; k < nelem; k++) {
                            stream->aggdata[i][j][h][k][0] += out_data[i][varid][h][k];
                        }
                    }
                    // Maximum over the period
                    else if (is_max) {
                        for (k = 0; k < nelem; k++) {
                            stream->aggdata[i][j][h][k][0] =
                                max(stream->aggdata[i][j][h][k][0], out_data[i][varid][h][k]);
                        }
                    }
                    // Minimum over the period
                    else if (is_min) {
                        for (k = 0; k < nelem; k++) {
                            stream->aggdata[i][j][h][k][0] =
                                min(stream->aggdata[i][j][h][k][0], out_data[i][varid][h][k]);
                        }
                    }
                    // Average over the period if counter is full
                    if (apply_avg_div) {
                        for (k = 0; k < nelem; k++) {
                            stream->aggdata[i][j][h][k][0] /= (double) alarm->count;
                        }
                    }
                }
            }
            // CELL output: first aggregate all HRUs using Cv.
            else {
                // 先初始化 cell_value 数组
                if (!is_first) {
                    for (k = 0; k < max_nelem; k++) {
                        cell_value[k] = 0;
                    }
                }
                max_nelem = 0;
                for (h = 0; h < nveg; h++) {
                    nelem = get_output_nelem(varid, i, h);
                    if (nelem > max_nelem) {
                        max_nelem = nelem;
                    }
                    double Cv = veg_con[i][h].Cv;
                    for (k = 0; k < nelem; k++) {
                        cell_value[k] += out_data[i][varid][h][k] * Cv;
                    }
                    is_first = false;
                }
                for (k = 0; k < max_nelem; k++) {
                    // The CELL result is stored at h = 0.
                    if (is_end) {
                        stream->aggdata[i][j][0][k][0] = cell_value[k];
                    }
                    else if (is_beg) {
                        stream->aggdata[i][j][0][k][0] = cell_value[k];
                    }
                    else if (is_sum || is_avg) {
                        stream->aggdata[i][j][0][k][0] += cell_value[k];
                    }
                    else if (is_max) {
                        stream->aggdata[i][j][0][k][0] = max(stream->aggdata[i][j][0][k][0], cell_value[k]);
                    }
                    else if (is_min) {
                        stream->aggdata[i][j][0][k][0] = min(stream->aggdata[i][j][0][k][0], cell_value[k]);
                    }
                    if (apply_avg_div) {
                        stream->aggdata[i][j][0][k][0] /= (double) alarm->count;
                    }
                }      
            }
        }
    }
}
