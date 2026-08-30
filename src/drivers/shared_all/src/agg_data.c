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
    extern metadata_struct out_metadata[N_OUTVAR_TYPES];

    alarm_struct          *alarm;
    size_t                 i, j, h, k;
    size_t                 nelem;
    size_t                 nveg;
    unsigned int           varid;
    bool                   alarm_now;
    double                 Cv;
    double cell_value[MAX_NODES] = {0};
    alarm = &(stream->agg_alarm);
    alarm->count++;
    alarm_now = raise_alarm(alarm, dmy_current);

    if (alarm->count == 1) {
        stream->time_bounds[0] = *dmy_current;
    }

    if (alarm_now) {
        stream->time_bounds[1] = *dmy_current;
    }

    for (i = 0; i < stream->ngridcells; i++) {
        nveg = stream->nveg[i];
        for (j = 0; j < stream->nvars; j++) {
            varid = stream->varid[j];
            nelem = out_metadata[varid].nelem;

            /* * CELL variables have only one output value per cell. 
               * HRU variables have one output value for each HRU. */ 
            if (stream->domain[j] == OUT_DOMAIN_HRU) {

                for (h = 0; h < nveg; h++) {
                    // Instantaneous at the end of the period
                    if ((stream->aggtype[j] == AGG_TYPE_END) && (alarm_now)) {
                        for (k = 0; k < nelem; k++) {
                            stream->aggdata[i][j][h][k][0] = out_data[i][varid][h][k];
                        }
                    }
                    // Instantaneous at the beginning of the period
                    else if ((stream->aggtype[j] == AGG_TYPE_BEG) &&
                            (alarm->count == 1)) {
                        for (k = 0; k < nelem; k++) {
                            stream->aggdata[i][j][h][k][0] = out_data[i][varid][h][k];
                        }
                    }
                    // Sum over the period
                    else if ((stream->aggtype[j] == AGG_TYPE_SUM) ||
                            (stream->aggtype[j] == AGG_TYPE_AVG)) {
                        for (k = 0; k < nelem; k++) {
                            stream->aggdata[i][j][h][k][0] += out_data[i][varid][h][k];
                        }
                    }
                    // Maximum over the period
                    else if (stream->aggtype[j] == AGG_TYPE_MAX) {
                        for (k = 0; k < nelem; k++) {
                            stream->aggdata[i][j][h][k][0] =
                                max(stream->aggdata[i][j][h][k][0], out_data[i][varid][h][k]);
                        }
                    }
                    // Minimum over the period
                    else if (stream->aggtype[j] == AGG_TYPE_MIN) {
                        for (k = 0; k < nelem; k++) {
                            stream->aggdata[i][j][h][k][0] =
                                min(stream->aggdata[i][j][h][k][0], out_data[i][varid][h][k]);
                        }
                    }
                    // Average over the period if counter is full
                    if ((stream->aggtype[j] == AGG_TYPE_AVG) && (alarm_now)) {
                        for (k = 0; k < nelem; k++) {
                            stream->aggdata[i][j][h][k][0] /= (double) alarm->count;
                        }
                    }
                }
            }
            // CELL output: first aggregate all HRUs using Cv.
            else {
                for (k = 0; k < nelem; k++) {
                    for (h = 0; h < nveg; h++) {
                        cell_value[k] += out_data[i][varid][h][k] * veg_con[i][h].Cv;
                    }
                    // The CELL result is stored at h = 0.
                    if ((stream->aggtype[j] == AGG_TYPE_END) && alarm_now) {
                        stream->aggdata[i][j][0][k][0] = cell_value[k];
                    }
                    else if ((stream->aggtype[j] == AGG_TYPE_BEG) && (alarm->count == 1)) {
                        stream->aggdata[i][j][0][k][0] = cell_value[k];
                    }
                    else if ((stream->aggtype[j] == AGG_TYPE_SUM) ||
                                (stream->aggtype[j] == AGG_TYPE_AVG)) {
                        stream->aggdata[i][j][0][k][0] += cell_value[k];
                    }
                    else if (stream->aggtype[j] == AGG_TYPE_MAX) {
                        stream->aggdata[i][j][0][k][0] = max(stream->aggdata[i][j][0][k][0], cell_value[k]);
                    }
                    else if (stream->aggtype[j] == AGG_TYPE_MIN) {
                        stream->aggdata[i][j][0][k][0] = min(stream->aggdata[i][j][0][k][0], cell_value[k]);
                    }
                    if ((stream->aggtype[j] == AGG_TYPE_AVG) && alarm_now) {
                        stream->aggdata[i][j][0][k][0] /= (double) alarm->count;
                    }
                }      
            }
        }
    }
}
