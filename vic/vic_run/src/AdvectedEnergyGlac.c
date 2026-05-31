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
AdvectedEnergyGlac(size_t             hidx,
                   force_data_struct *force,
                   energy_bal_struct *energy)

{
    double Tgrnd = energy->Tgrnd;
    double rainfall = force->rainf[hidx];
    double snowfall = force->snowf[hidx];
    double air_temp = force->air_temp[hidx];

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