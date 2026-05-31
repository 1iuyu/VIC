/******************************************************************************
* @section DESCRIPTION
*
* Calculates the interception and subsequent release of by the forest canopy
* using an energy balance approach.
******************************************************************************/

#include "vic_run.h"

/******************************************************************************
* @brief    Calculate snow interception and release by the canopy
******************************************************************************/
int
snow_intercept(size_t             hidx,
               double             step_dt,
               double             Tcanopy,  // canopy air temperature
               force_data_struct *force,
               snow_data_struct  *snow,
               veg_var_struct    *veg_var)
{
    extern option_struct     options;   
    extern parameters_struct param;

    double      wetFrac;
    double      DeltaSnowInt; /* Change in the physical swe of snow
                                              interceped on the branches. (m) */
    double      SnowDriptemp;
    double      SnowDripwind;
    double      SnowDrip;     /* Amount of drip from intercepted snow as a
                                                               result of snowmelt (m) */
    double      SnowUnload = 0.0;   /* Amount of snow unloaded by snowfall (m) */
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
    double wind = force->wind[hidx];
    double snowfall = force->snowf[hidx];
    double rainfall = force->rainf[hidx];
    /* Determine the maximum canopy capacity for snow interception [mm] */
    MaxSnowInt = fcanopy * 6.6 * (0.27 + 46.0 / new_snow_density) *
                                                    (NetLAI + NetSAI);

    /* Calculate snow interception. */
    if (NetLAI + NetSAI > 0.) {
        DeltaSnowInt = fcanopy * snowfall;
        DeltaSnowInt = min(DeltaSnowInt, (MaxSnowInt - int_snow) / step_dt * 
                            (1.0 - exp(-snowfall * step_dt / MaxSnowInt)));
        DeltaSnowInt = max(DeltaSnowInt, 0.0);
        /* Reduce the amount of intercepted snow if windy and cold. */
        SnowDriptemp = max(0.0, (Tcanopy - 270.15) / 1.87e5);
        SnowDripwind = max(0.0, wind / 1.56e5);
        SnowUnload = max(0, int_snow) * (SnowDriptemp + SnowDripwind);
        SnowUnload = min(SnowUnload, int_snow / step_dt + DeltaSnowInt);
        SnowDrip = (snowfall * fcanopy - DeltaSnowInt);
        SnowThroughFall = snowfall * (1 - fcanopy);

        /* now update snowfall and total accumulated intercepted snow amounts */
        int_snow += (DeltaSnowInt - SnowUnload) * step_dt;
        if (int_snow < 0.) {
            int_snow = 0.;
        }
    }
    else {
        DeltaSnowInt = 0.;
        SnowDrip = 0.;
        SnowThroughFall = snowfall;
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
        DeltaRainInt = fcanopy * rainfall;
        DeltaRainInt = min(DeltaRainInt, (MaxRainInt - int_rain) / step_dt *
                                        (1.0 - exp(-rainfall * step_dt / MaxRainInt)));
        DeltaRainInt = max(DeltaRainInt, 0.0);
        RainDrip = rainfall * fcanopy - DeltaRainInt;
        RainThroughFall = rainfall * (1 - fcanopy);
        int_rain += DeltaRainInt * step_dt;
        if (int_rain < 0.) {
            int_rain = 0.;
        }
    }
    else {
        DeltaRainInt = 0.;
        RainDrip = 0.;
        RainThroughFall = rainfall;
        /* canopy gets buried by rain */
        if (int_rain > 0.) {
            RainDrip += int_rain / step_dt;
            int_rain = 0.; 
        }
    }

    /* wetted fraction of canopy */
    if (options.CANOPY_INTERCEP == NOAH) {
        if (int_snow > 0.0) {
            wetFrac = max(0.0, int_snow) / max(MaxSnowInt, param.TOL_A);
        }
        else {
            wetFrac = max(0.0, int_rain) / max(MaxRainInt, param.TOL_A);
        }
        veg_var->wetFrac = pow(min(wetFrac, 1.0), 0.667);
        veg_var->dryFrac = (1.0 - wetFrac) * NetLAI / (NetLAI + NetSAI);
    }
    else if (options.CANOPY_INTERCEP == CLM) {
        double f_snow = 0.0;
        double dryFrac = 0.0;
        if (NetLAI + NetSAI > 0.0) {
            if (int_rain + int_snow > 0.0) {
                wetFrac = pow(((int_rain + int_snow) / (NetLAI + NetSAI)), 0.66666);
                if (int_snow > 0.0) {
                    f_snow = pow((int_snow / (NetLAI + NetSAI)), 0.15);
                    if (f_snow > 1.0) {
                        f_snow = 1.0;
                    }
                }
                else {
                    f_snow = 0.0;
                }
            }
            else {
                wetFrac = 0.0;
                f_snow = 0.0;
            }
            dryFrac = (1.0 - wetFrac) * NetLAI / (NetLAI + NetSAI);
        }
        else {
            dryFrac = 0.0;
            wetFrac = 0.0;
        }
        veg_var->wetFrac = wetFrac;
        veg_var->dryFrac = dryFrac;
    }
    else {
        log_err("Unknown CANOPY_INTERCEP option");
    }
    /* Return structure */
    veg_var->canopy_swq = int_snow + int_rain;
    veg_var->int_snow = int_snow;
    veg_var->int_rain = int_rain;
    veg_var->RainThroughFall = RainThroughFall;
    veg_var->SnowThroughFall = SnowThroughFall;
    veg_var->RainDrip = RainDrip;
    veg_var->SnowDrip = SnowDrip;
    veg_var->MaxSnowInt = MaxSnowInt;
    veg_var->MaxRainInt = MaxRainInt;
    veg_var->SnowUnload = SnowUnload;
    /* rain or snow on the ground [mm/s] */
    force->rainf[hidx] = RainDrip + RainThroughFall;
    force->snowf[hidx] = SnowDrip + SnowThroughFall + SnowUnload;
    /* Update snow water equivalent and snow depth */
    if (force->snowf[hidx] > 0.0) {
        snow->delta_depth = force->snowf[hidx] / new_snow_density;
    }
    else {
        snow->delta_depth = 0.0;
    }
    
    return (0);
}
