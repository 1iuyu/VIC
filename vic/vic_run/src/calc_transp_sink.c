/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine computes the water sink due to plant transpiration. 
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief  
 * Compute the water sink due to plant transpiration. 
 *****************************************************************************/
int
calc_transp_sink(cell_data_struct *cell,
                 soil_con_struct  *soil_con,
                 veg_var_struct   *veg_var)
{
    double grav2 = 0.0;
    double *transp_sink = cell->transp_sink;
    double *hksr_int = cell->hksr_int;
    double *vegwp = veg_var->vegwp;
    double *zc_soil = soil_con->zc_soil;
    double *matric = cell->matric;
    // hydraulic redistribution duo to plant transpiration
    for (size_t i = 0; i < cell->Nroot; i++) {
        grav2 = zc_soil[i];
        transp_sink[i] = hksr_int[i] * (matric[i] - vegwp[3] - grav2);
    }
    
    return (0);
}