/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine calculate albedo of snow containing impurities and the 
 * evolution of snow effective radius.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief  
 * Compute the albedo of snow containing impurities and the evolution of 
 * snow effective radius.
 *****************************************************************************/
int
snow_SNICAR(double             step_dt,
            double             coszen,
            double             snowfall,
            energy_bal_struct *energy,
            cell_data_struct  *cell,
            snow_data_struct  *snow,
			soil_con_struct   *soil_con)
{
    //
}