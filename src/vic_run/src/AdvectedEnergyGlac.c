/******************************************************************************
* @section DESCRIPTION
*
* Calculates the heat flux advected from precipitation to glacier ground
******************************************************************************/

#include "vic_run.h"

/******************************************************************************
* @brief    Calculates the heat flux advected from precipitation
******************************************************************************/
void
AdvectedEnergyGlac(double             air_temp,
                   double             snowfall,
                   double             rainfall,
                   energy_bal_struct *energy)

{
    double Tgrnd = energy->Tgrnd;
    /* initialization */
    double AdvectGrnd = 0.;
    
    /* Heat advection for liquid rainfall */
    AdvectGrnd = rainfall / MM_PER_M *
                            CONST_CPFWICE * (air_temp - Tgrnd);

    /* Heat advection for snowfall */
    AdvectGrnd += snowfall / MM_PER_M *
                            CONST_CPICE * (air_temp - Tgrnd);

    /* Put some artificial limits here for stability */
    energy->AdvectGrnd = max(-20.0, min(20.0, AdvectGrnd));
}