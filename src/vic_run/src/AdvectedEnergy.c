/******************************************************************************
* @section DESCRIPTION
*
* Calculates the heat flux advected from precipitation to vegetation and ground
******************************************************************************/

#include "vic_run.h"

/******************************************************************************
* @brief    Calculates the heat flux advected from precipitation
******************************************************************************/
void
AdvectedEnergy(double             air_temp,
               double             snowfall,
               double             rainfall,
               energy_bal_struct *energy,
               veg_var_struct    *veg_var)

{
    double fcanopy = veg_var->fcanopy;
    double Tgrnd = energy->Tgrnd;
    double Tfoliage = energy->Tfoliage;
    double RainThroughFall = veg_var->RainThroughFall;
    double SnowThroughFall = veg_var->SnowThroughFall;
    double RainDrip = veg_var->RainDrip;
    double SnowDrip = veg_var->SnowDrip;

    /* initialization */
    double AdvectOver = 0.;
    double AdvectGrnd = 0.;
    double AdvectSub = 0.;
    double Advectedair2can = 0.;
    double Advectedair2grnd = 0.;
    double Advectedcan2grnd = 0.;

    /* Heat advection for liquid rainfall */
    Advectedair2can = rainfall * fcanopy * 
                                    (CONST_CPFWICE / 1000) * (air_temp - Tfoliage);
    Advectedair2grnd = RainThroughFall * 
                                    (CONST_CPFWICE / 1000) * (air_temp - Tgrnd);
    Advectedcan2grnd = RainDrip * 
                                    (CONST_CPFWICE / 1000) * (Tfoliage - Tgrnd);

    /* Heat advection for snowfall */
    Advectedair2can += snowfall * fcanopy * 
                                        (CONST_CPICE / 1000) * (air_temp - Tfoliage);
    Advectedair2grnd += SnowThroughFall * 
                                        (CONST_CPICE / 1000) * (air_temp - Tgrnd);
    Advectedcan2grnd += SnowDrip * (CONST_CPICE / 1000) * (Tfoliage - Tgrnd); //???

    /* net heat advection */
    AdvectOver = Advectedair2can - Advectedcan2grnd;
    AdvectGrnd = Advectedair2grnd;
    AdvectSub = Advectedcan2grnd;

    /* adjust for fcanopy */
    if (fcanopy > 0.0 && fcanopy < 1.0) {
        AdvectGrnd /= fcanopy;
        AdvectSub /= (1 - fcanopy); 
    }
    else if (fcanopy <= 0) {
        AdvectOver = 0.0;
        AdvectSub = 0.0;
        AdvectGrnd += AdvectSub;
    }
    else {
        AdvectGrnd = 0.0;
    }

    /* Put some artificial limits here for stability */
    energy->AdvectOver = max(-20.0, min(20.0, AdvectOver));
    energy->AdvectGrnd = max(-20.0, min(20.0, AdvectGrnd));
    energy->AdvectSub = max(-20.0, min(20.0, AdvectSub));
}