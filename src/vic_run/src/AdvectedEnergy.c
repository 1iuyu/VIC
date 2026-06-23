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
               cell_data_struct  *cell,
               veg_var_struct    *veg_var)
{
    /* initialization */
    double AdvectOver = 0.0;
    double AdvectGrnd = 0.0;
    double AdvectSub = 0.0;
    double fcanopy = veg_var->fcanopy;
    double Tgrnd = energy->Tgrnd;
    
    if (cell->IS_VEG) {
        double fcanopy = veg_var->fcanopy;
        double Tfoliage = energy->Tfoliage;
        double RainThroughFall = veg_var->RainThroughFall;
        double SnowThroughFall = veg_var->SnowThroughFall;
        double RainDrip = veg_var->RainDrip;
        double SnowDrip = veg_var->SnowDrip;
        double Advectedair2can = 0.0;
        double Advectedair2grnd = 0.0;
        double Advectedcan2grnd = 0.0;
        /* Heat advection for liquid rainfall */
        Advectedair2can = rainfall * fcanopy * 
                                        (CONST_CPFWICE / MM_PER_M) * (air_temp - Tfoliage);
        Advectedair2grnd = RainThroughFall * 
                                        (CONST_CPFWICE / MM_PER_M) * (air_temp - Tgrnd);
        Advectedcan2grnd = RainDrip * 
                                        (CONST_CPFWICE / MM_PER_M) * (Tfoliage - Tgrnd);

        /* Heat advection for snowfall */
        Advectedair2can += snowfall * fcanopy * 
                                            (CONST_CPICE / MM_PER_M) * (air_temp - Tfoliage);
        Advectedair2grnd += SnowThroughFall * 
                                            (CONST_CPICE / MM_PER_M) * (air_temp - Tgrnd);
        Advectedcan2grnd += SnowDrip * (CONST_CPICE / MM_PER_M) * (Tfoliage - Tgrnd); //???

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
    }
    else if (cell->IS_GLAC || cell->IS_WET) {
        /* Heat advection for liquid rainfall */
        AdvectGrnd = rainfall / MM_PER_M * CONST_CPFWICE * (air_temp - Tgrnd);
        /* Heat advection for snowfall */
        AdvectGrnd += snowfall / MM_PER_M * CONST_CPICE * (air_temp - Tgrnd);       
    }

    /* Put some artificial limits here for stability */
    energy->AdvectOver = max(-20.0, min(20.0, AdvectOver));
    energy->AdvectGrnd = max(-20.0, min(20.0, AdvectGrnd));
    energy->AdvectSub = max(-20.0, min(20.0, AdvectSub));
    if (cell->IS_VEG) {
        energy->advection = fcanopy * energy->AdvectSub + (1.0 - fcanopy) * 
                            energy->AdvectGrnd + energy->AdvectOver;
    }
    else {
        energy->advection = energy->AdvectGrnd;
    }

}