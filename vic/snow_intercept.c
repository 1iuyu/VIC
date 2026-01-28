/******************************************************************************
* @section DESCRIPTION
*
* Calculates the interception and subsequent release of by the forest canopy
* using an energy balance approach.
******************************************************************************/

#include <vic_run.h>

/******************************************************************************
* @brief    Calculate snow interception and release by the canopy
******************************************************************************/
int
snow_intercept(double            step_dt,
               double            Tcanopy,  // canopy air temperature
               double           *SnowFall,
               double           *RainFall,
               double            Wind,
               snow_data_struct *snow,
               veg_var_struct   *veg_var)
{

    extern parameters_struct param;

    double      wetFrac;
    double      DeltaSnowInt; /* Change in the physical swe of snow
                                              interceped on the branches. (m) */
    double      SnowDriptemp;
    double      SnowDripwind;
    double      SnowDrip;     /* Amount of drip from intercepted snow as a
                                                               result of snowmelt (m) */
    double      tmp_SnowDrip;
    double      RainDrip;
    double      MaxRainInt; /* rain interception capacity (mm) */
    double      MaxSnowInt; /* snow interception capacity (mm) */
    double      DeltaRainInt;
    double      RainThroughFall; /* Amount of rain reaching to the ground (m)
                                               */
    double      SnowThroughFall; /* Amount of snow reaching to the ground (m)
                                               */
    /* Initialization */
    double Wdew = veg_var->Wdew;
    double NetLAI = veg_var->NetLAI;
    double NetSAI = veg_var->NetSAI;
    double fcanopy = veg_var->fcanopy;
    double int_snow = veg_var->int_snow;    // intercepted rain
    double int_rain = veg_var->int_rain;    // intercepted snow
    double new_snow_density = snow->new_snow_density;

    /* Determine the maximum canopy capacity for snow interception [mm] */
    MaxSnowInt = fcanopy * 6.6 * (0.27 + 46.0 / new_snow_density) *
                                                    (NetLAI + NetSAI);

    /* Calculate snow interception. */
    if (NetLAI + NetSAI > 0.) {
        DeltaSnowInt = fcanopy * (*SnowFall);
        DeltaSnowInt = min(DeltaSnowInt, (MaxSnowInt - int_snow) / step_dt * 
                               (1.0 - exp(-*SnowFall * step_dt / MaxSnowInt)));
        if (DeltaSnowInt < 0.) {
            DeltaSnowInt = 0.;
        }
        /* Reduce the amount of intercepted snow if windy and cold. */
        SnowDriptemp = (Tcanopy - 270.15) / 1.87e5;
        SnowDripwind = Wind / 1.56e5;
        tmp_SnowDrip = max(0, int_snow) * (SnowDriptemp + SnowDripwind);
        if (tmp_SnowDrip > int_snow / step_dt + DeltaSnowInt) {
            tmp_SnowDrip = int_snow / step_dt + DeltaSnowInt;
        }
        SnowDrip = (*SnowFall * fcanopy - DeltaSnowInt) + tmp_SnowDrip;
        SnowThroughFall = (*SnowFall) * (1 - fcanopy);

        /* now update snowfall and total accumulated intercepted snow amounts */
        int_snow += (DeltaSnowInt - tmp_SnowDrip) * step_dt;
        if (int_snow < 0.) {
           int_snow = 0.;
        }
    }
    else {
        DeltaSnowInt = 0.;
        SnowDrip = 0.;
        SnowThroughFall = (*SnowFall);
        /* canopy gets buried by snow */
        if (int_snow > 0.) {
            SnowDrip += int_snow / step_dt;
            int_snow = 0.;
        }
    }

    /* Calculate amount of rain intercepted on branches and stored in
       intercepted snow. */
    MaxRainInt = fcanopy * Wdew * (NetLAI + NetSAI);
    if ((NetLAI + NetSAI) > 0.) {
        DeltaRainInt = fcanopy * (*RainFall);
        DeltaRainInt = min(DeltaRainInt, (MaxRainInt - int_rain) / step_dt *
                                        (1.0 - exp(-*RainFall * step_dt / MaxRainInt)));
        if (DeltaRainInt < 0.) {
            DeltaRainInt = 0.;
        }

        RainDrip = (*RainFall) * fcanopy - DeltaRainInt;
        RainThroughFall = (*RainFall) * (1 - fcanopy);
        int_rain += DeltaRainInt * step_dt;
        if (int_rain < 0.) {
            int_rain = 0.;
        }
    }
    else {
        DeltaRainInt = 0.;
        RainDrip = 0.;
        RainThroughFall = (*RainFall);
        /* canopy gets buried by rain */
        if (int_rain > 0.) {
            RainDrip += int_rain / step_dt;
            int_rain = 0.; 
        }
    }

    /* Return structure */
    snow->snow_canopy = int_snow + int_rain;
    veg_var->int_snow = int_snow;
    veg_var->int_rain = int_rain;
    veg_var->RainThroughFall = RainThroughFall;
    veg_var->SnowThroughFall = SnowThroughFall;
    veg_var->RainDrip = RainDrip;
    veg_var->SnowDrip = SnowDrip;
    veg_var->MaxSnowInt = MaxSnowInt;
    veg_var->MaxRainInt = MaxRainInt;
    
    /* rain or snow on the ground [mm/s] */
    (*RainFall) = RainDrip + RainThroughFall; 

    (*SnowFall) = SnowDrip + SnowThroughFall;

    /* wetted fraction of canopy */
    if (int_snow > 0.0) {
        wetFrac = max(0.0, int_snow) / max(MaxSnowInt, param.TOL_A);
    }
    else {
        wetFrac = max(0.0, int_rain) / max(MaxRainInt, param.TOL_A);
    }
    veg_var->wetFrac = pow(min(wetFrac, 1.0), 0.667);

    /* Update snow water equivalent and snow depth */
    snow->delta_depth = (*SnowFall) / new_snow_density;

    return (0);
}
