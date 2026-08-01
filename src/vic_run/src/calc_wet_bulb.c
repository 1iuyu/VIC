/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine calculates the wet-bulb temperature based on the air 
 * temperature and relative humidity.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief  This subroutine calculates the wet-bulb temperature.
 *****************************************************************************/
int
calc_wet_bulb(double  air_temp,
              double  pressure,
              double  rel_humid,
              double  Qair,
              double *wet_bulb)
{

    *wet_bulb = air_temp + pressure + rel_humid + Qair;
    
    return (0);
}