/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate snow pack energy balance
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief  
 * Solves the glacier energy balance. See SnowPackEnergyBalance for comparison.
 *****************************************************************************/
double
GlacierEnergyBalance(double  TSurf,
                      va_list ap)
{
    extern option_struct     options;
    extern parameters_struct param;
    /* Define Variable Argument List */
    double Dt;                      /* Model time step (sec) */
    double Ra;                      /* Aerodynamic resistance (s/m) */
    double *Ra_used;                /* Aerodynamic resistance (s/m) after stability correction */

    /* Vegetation Parameters */
    double Z;                       /* Reference height (m) */
    double *Z0;                     /* surface roughness height (m) */

    /* Atmospheric Forcing Variables */
    double AirDens;                 /* Density of air (kg/m3) */
    double EactAir;                 /* Actual vapor pressure of air (Pa) */
    double LongSnowIn;              /* Incoming longwave radiation (W/m2) */
    double Lv;                      /* Latent heat of vaporization (J/kg3) */
    double Press;                   /* Air pressure (Pa) */
    double Rain;                    /* Rain fall (m/timestep) */
    double NetShortUnder;           /* Net incident shortwave radiation (W/m2) */
    double Vpd;                     /* Vapor pressure deficit (Pa) */
    double Wind;                    /* Wind speed (m/s) */

    /* Snowpack Variables */
    double OldTSurf;                /* Surface temperature during previous time step */
    double IceDepth;                /* Depth of glacier surface layer (m) */
    double IceWE;                   /* Liquid water in the glacier surface layer (m) */

    /* Energy Balance Components */
    double Tair;                    /* Canopy air / Air temperature (C) */
    double TGrnd;                   /* Ground surface temperature (C) */

    double *AdvectedEnergy;         /* Energy advected by precipitation (W/m2) */
    double *DeltaColdContent;       /* Change in cold content of surface layer (W/m2) */
    double *GroundFlux;             /* Ground Heat Flux (W/m2) */
    double *LatentHeat;             /* Latent heat exchange at surface (W/m2) */
    double *LatentHeatSub;          /* Latent heat of sublimation exchange at surface (W/m2) */
    double *NetLongUnder;           /* Net longwave radiation at snowpack surface (W/m^2) */
    double *SensibleHeat;           /* Sensible heat exchange at surface (W/m2) */
    double *vapor_flux;              /* Mass flux of water vapor to or from 
                                    the intercepted snow (m/timestep) */
    /* Internal Routine Variables */
    double Density;                 /* Density of water/ice at TSurf (kg/m3) */
    double NetRad;                  /* Net radiation exchange at surface (W/m2) */
    double RestTerm;                /* Rest term in surface energy balance (W/m2) */
    double TMean;                   /* Average ice surface layer temperature for time step (C) */
    double OldTMean;
    double Tmp;
    double VaporMassFlux;           /* Mass flux of water vapor to or from the intercepted snow (kg/m2s) */
    double Fbal;                    /* Energy balance at glacier surface */
    double temp_IceDepth;           /* Local parameter to hold scaled value for IceDepth */

    /* General Model Parameters */
    Dt = (double) va_arg(ap, double);
    Ra = (double) va_arg(ap, double);
    Ra_used = (double *) va_arg(ap, double *);

    /* Vegetation Parameters */
    Z = (double) va_arg(ap, double);
    Z0 = (double *) va_arg(ap, double *);

    /* Atmospheric Forcing Variables */
    AirDens = (double) va_arg(ap, double);
    EactAir = (double) va_arg(ap, double);
    LongSnowIn = (double) va_arg(ap, double);
    Lv = (double) va_arg(ap, double);
    Press = (double) va_arg(ap, double);
    Rain = (double) va_arg(ap, double);
    NetShortUnder = (double) va_arg(ap, double);
    Vpd = (double) va_arg(ap, double);
    Wind = (double) va_arg(ap, double);

    /* Snowpack Variables */
    OldTSurf = (double) va_arg(ap, double);
    IceDepth = (double) va_arg(ap, double);
    IceWE = (double) va_arg(ap, double);

    /* Energy Balance Components */
    Tair = (double) va_arg(ap, double);
    TGrnd = (double) va_arg(ap, double);

    AdvectedEnergy = (double *) va_arg(ap, double *);
    DeltaColdContent = (double *) va_arg(ap, double *);
    GroundFlux = (double *) va_arg(ap, double *);
    LatentHeat = (double *) va_arg(ap, double *);
    LatentHeatSub = (double *) va_arg(ap, double *);
    NetLongUnder = (double *) va_arg(ap, double *);
    SensibleHeat = (double *) va_arg(ap, double *);
    vapor_flux = (double *) va_arg(ap, double *);

    /* Calculate active temp for energy balance as average of old and new  */

    TMean = (TSurf + TGrnd) / 2;
    OldTMean = (OldTSurf + TGrnd) / 2;
    Density = CONST_RHOFW;
    temp_IceDepth = IceDepth / MM_PER_M;

    /* Correct aerodynamic conductance for stable conditions
       Note: If air temp >> glacier temp then aero_cond -> 0 (i.e. very stable)
       velocity (vel_2m) is expected to be in m/sec */

    /* Apply the stability correction to the aerodynamic resistance
       NOTE: In the old code 2m was passed instead of Z-Displacement.  I (bart)
       think that it is more correct to calculate ALL fluxes at the same
       reference level */

    if (Wind > 0.0) {
        Ra_used[0] = Ra / StabilityCorrection(Z, 0.f, TSurf, Tair, Wind, Z0[2]);
    }
    else {
        Ra_used[0] = param.HUGE_RESIST;
    }

    /* Calculate longwave exchange and net radiation */

    Tmp = TSurf + CONST_TKFRZ;
    (*NetLongUnder) = LongSnowIn - calc_outgoing_longwave(Tmp,
                                                          param.EMISS_SNOW);
    NetRad = NetShortUnder + (*NetLongUnder);

    /* Calculate the sensible heat flux */
    *SensibleHeat = calc_sensible_heat(AirDens, Tair, TSurf, Ra_used[0]);

    /* Convert sublimation terms from m/timestep to kg/m2s */
    VaporMassFlux = *vapor_flux * Density / Dt;

    /* Calculate the mass flux of ice to or from the surface layer */

    /* Calculate the saturated vapor pressure,
       (Equation 3.32, Bras 1990) */
     latent_heat_from_glacier(AirDens, EactAir, Lv, Press, Ra_used[0], TSurf, 
                              Vpd, LatentHeat, LatentHeatSub, &VaporMassFlux);

    /* Convert sublimation terms from kg/m2s to m/timestep */
    *vapor_flux = VaporMassFlux * Dt / Density;

    /* Calculate advected heat flux from rain
       Equation 7.3.12 from H.B.H. for rain falling on melting snowpack */
    if (TSurf == 0) {
        *AdvectedEnergy = (CONST_CPFW * CONST_RHOFW * (Tair) * Rain) / (Dt);
    }
    else {
        *AdvectedEnergy = 0.;
    }

    /* Calculate change in cold content */
    *DeltaColdContent = CONST_VCPICE_WQ * temp_IceDepth * 
                        (TMean - OldTMean) / (Dt);

    /* Calculate Ground Heat Flux */
    /* Estimate of ice thermal conductivity (at atmospheric pressure) adapted from Slack (1980), Table 1; assumes
         linear relationship between TSurf and K above -75C */
    *GroundFlux = (GLAC_K_ICE + TSurf * (-0.0142)) * (TGrnd - TSurf) / temp_IceDepth;

    /* Calculate energy balance error at the glacier surface */
    Fbal = NetRad + (*SensibleHeat) + (*LatentHeat) + (*LatentHeatSub) + (*AdvectedEnergy);
    RestTerm = Fbal - *DeltaColdContent + *GroundFlux;

    /* Melting occurs when surface at melting point and surface energy flux is positive */
    if (TSurf == 0.0 && RestTerm >= 0.) {
      RestTerm = 0.;
    }

    return RestTerm;
}