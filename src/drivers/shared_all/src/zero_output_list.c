/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine resets the values of all output variables to 0.
 *****************************************************************************/

#include "vic_driver_shared_all.h"

/******************************************************************************
 * @brief    This routine resets the values of all output variables to 0.
 *****************************************************************************/
void
zero_output_list(size_t         nveg,
                 double      ***out_data,
                 stream_struct *streams)
{
    extern metadata_struct out_metadata[N_OUTVAR_TYPES];

    size_t                 varid, i, j, k;

    for (i = 0; i < streams->nvars; i++) {
        varid = streams->varid[i];
        for (j = 0; j < nveg; j++) {
            for (k = 0; k < out_metadata[varid].nelem; k++) {
                out_data[varid][j][k] = 0.0;
            }
        }
    }
}
