/******************************************************************************
* @section DESCRIPTION
*
* Calculates the heat flux advected from precipitation to glacier ground
******************************************************************************/

#include <vic_run.h>

/******************************************************************************
* @brief    Calculates the heat flux advected from precipitation
******************************************************************************/
void
AdvectedEnergyGlac(double   Tgrnd,            // ground temperature [K]
                   double   air_temp, 
                   double   RainFall,
                   double   SnowFall,
                   double  *AdvectGrnd)

{

    /* initialization */
    *AdvectGrnd = 0.;
    
    /* Heat advection for liquid rainfall */
    *AdvectGrnd = RainFall / MM_PER_M *
                            CONST_CPFWICE * (air_temp - Tgrnd);

    /* Heat advection for snowfall */
    *AdvectGrnd += SnowFall / MM_PER_M *
                            CONST_CPICE * (air_temp - Tgrnd);

    /* Put some artificial limits here for stability */
    *AdvectGrnd  = max(-20.0, min(20.0, *AdvectGrnd));
}