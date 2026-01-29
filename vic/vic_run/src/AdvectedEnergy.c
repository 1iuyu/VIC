/******************************************************************************
* @section DESCRIPTION
*
* Calculates the heat flux advected from precipitation to vegetation and ground
******************************************************************************/

#include <vic_run.h>

/******************************************************************************
* @brief    Calculates the heat flux advected from precipitation
******************************************************************************/
void
AdvectedEnergy(double    fcanopy,
               double    Tcanopy,          // canopy air temperature
               double    Tgrnd,            // soil surface temperature
               double    air_temp, 
               double    RainFall,
               double    SnowFall,
               double    RainThroughFall,
               double    SnowThroughFall,
               double    RainDrip,
               double    SnowDrip,
               double   *AdvectOver,
               double   *AdvectGrnd,
               double   *AdvectSub)

{
    double      Advectedair2can;
    double      Advectedair2grnd;
    double      Advectedcan2grnd;

    /* initialization */
    (*AdvectOver) = 0.;
    (*AdvectGrnd) = 0.;
    (*AdvectSub) = 0.;
    Advectedair2can = 0.;
    Advectedair2grnd = 0.;
    Advectedcan2grnd = 0.;

    /* Heat advection for liquid rainfall */
    Advectedair2can = RainFall * fcanopy * 
                                    (CONST_CPFWICE / 1000) * (air_temp - Tcanopy);
    Advectedair2grnd = RainThroughFall * 
                                    (CONST_CPFWICE / 1000) * (air_temp - Tgrnd);
    Advectedcan2grnd = RainDrip * 
                                    (CONST_CPFWICE / 1000) * (Tcanopy - Tgrnd);

    /* Heat advection for snowfall */
    Advectedair2can += SnowFall * fcanopy * 
                                        (CONST_CPICE / 1000) * (air_temp - Tcanopy);
    Advectedair2grnd += SnowThroughFall * 
                                        (CONST_CPICE / 1000) * (air_temp - Tgrnd);
    Advectedcan2grnd += SnowDrip * (CONST_CPICE / 1000) * (Tcanopy - Tgrnd); //???

    /* net heat advection */
    (*AdvectOver) = Advectedair2can - Advectedcan2grnd;
    (*AdvectGrnd) = Advectedair2grnd;
    (*AdvectSub) = Advectedcan2grnd;

    /* adjust for fcanopy */
    if (fcanopy > 0. && fcanopy < 1.) {
        (*AdvectGrnd) /= fcanopy;
        (*AdvectSub) /= (1 - fcanopy); 
    }
    else if (fcanopy <= 0) {
        (*AdvectOver) = 0.;
        (*AdvectSub) = 0.;
        (*AdvectGrnd) += (*AdvectSub);
    }
    else {
        (*AdvectGrnd) = 0.;
    }

    /* Put some artificial limits here for stability */
    (*AdvectOver) = max(-20.0, min(20.0, *AdvectOver));
    (*AdvectGrnd) = max(-20.0, min(20.0, *AdvectGrnd));
    (*AdvectSub) = max(-20.0, min(20.0, *AdvectSub));
}