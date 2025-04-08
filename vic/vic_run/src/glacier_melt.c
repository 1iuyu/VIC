/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate glacier accumulation and melt
 *****************************************************************************/

#include <vic_run.h>

/*****************************************************************************
  * @brief   Calculate glacier accumulation and melt using an energy balance
  *          approach for a two layer model
 *****************************************************************************/
int
glacier_melt(double            Le,
             double            NetShortGlac,  // net SW at absorbed by glacier
             double            Tgrnd,
             double           *Z0,  // roughness
             double            aero_resist,  // aerodynamic resistance
             double           *aero_resist_used,  // stability-corrected aerodynamic resistance
             double            air_temp,  // air temperature
             double            delta_t,  // time step in secs
             double            density,  // atmospheric density
             double            displacement,  // surface displacement
             double            LongUnderIn,  // incoming longwave radiation
             double            pressure, 
             double            rainfall,
             double            vp,
             double            vpd,
             double            wind,
             double            z2,
             double           *NetLongGlac,
             double           *OldTSurf,
             double           *melt,
             double           *save_Qnet,
             double           *save_advection,
             double           *save_deltaCC_glac,
             double           *save_glacier_melt_energy,
             double           *save_grnd_flux,
             double           *save_latent,
             double           *save_latent_sub,
             double           *save_sensible,
             glac_data_struct *glacier,
             soil_con_struct  *soil_con,
             snow_data_struct *snow)
{
    extern option_struct     options;
    extern parameters_struct param;

    double                     error;
 //   double                     MassBalanceError;       /* Mass balance error (m) */
    double                     Qnet;                   /* Net energy exchange at the surface (W/m2) */
    double                     GlacMelt = 0.;               /* Amount of ice melt during time interval (m water equivalent) */
    double                     GlacCC = 0.;                 /* Cold content of glacier surface layer (J) */
    double                     RainFall;
    double                     advection;
    double                     deltaCC_glac;
    double                     latent_heat;
    double                     latent_heat_sub;
    double                     sensible_heat;
    double                     melt_energy = 0.;
    double                     grnd_flux;

    RainFall = rainfall / MM_PER_M; /* convert to m */

    (*OldTSurf) = glacier->surf_temp;

    /* Calculate the surface energy balance for surf_temp = 0.0 */
    Qnet = CalcGlacierEnergyBalance((double) 0.0, delta_t, aero_resist,
                                     aero_resist_used, z2, Z0,
                                     density, vp, LongUnderIn, Le, pressure,
                                     RainFall, NetShortGlac, vpd,
                                     wind, (*OldTSurf),
                                     soil_con->GLAC_SURF_THICK,
                                     soil_con->GLAC_SURF_WE,
                                     air_temp, Tgrnd,
                                     &advection,
                                     &deltaCC_glac,
                                     &grnd_flux, &latent_heat,
                                     &latent_heat_sub, NetLongGlac,
                                     &sensible_heat,
                                     &glacier->vapor_flux);

    /* If Qnet == 0.0, then set the surface temperature to 0.0 */
    if (Qnet == 0.0) {
      glacier->surf_temp = 0.;
      melt_energy = NetShortGlac + (*NetLongGlac) + sensible_heat
          + latent_heat + latent_heat_sub + advection - deltaCC_glac;
      GlacMelt = melt_energy / (CONST_LATICE * CONST_RHOFW) * delta_t;
      GlacCC = 0.;
    }

    /* Else, GlacierEnergyBalance(T=0.0) RestTerm <= 0.0 */
    else {
        /* Calculate surface layer temperature using "Brent method" */
        glacier->surf_temp = root_brent(
            (double) (glacier->surf_temp - param.SNOW_DT),
            (double) (glacier->surf_temp + param.SNOW_DT),
            GlacierEnergyBalance,
            delta_t, aero_resist, aero_resist_used, z2, Z0,
            density, vp, LongUnderIn, Le, pressure,
            RainFall, NetShortGlac, vpd,
            wind, (*OldTSurf),
            soil_con->GLAC_SURF_THICK,
            soil_con->GLAC_SURF_WE,
            air_temp, Tgrnd,
            &advection,
            &deltaCC_glac,
            &grnd_flux, &latent_heat,
            &latent_heat_sub, NetLongGlac,
            &sensible_heat,
            &glacier->vapor_flux);
        if (glacier->surf_temp <= -998) {
            if (options.TFALLBACK) {
                glacier->surf_temp = *OldTSurf;
                glacier->surf_temp_fbflag = 1;
                glacier->surf_temp_fbcount++;
            }

            else {           
                error = ErrorGlacierEnergyBalance(glacier->surf_temp,
                                                  delta_t,
                                                  aero_resist, 
                                                  aero_resist_used, 
                                                  displacement, z2, Z0,
                                                  density, vp, LongUnderIn, 
                                                  Le, pressure, RainFall, 
                                                  NetShortGlac, vpd, wind,
                                                  (*OldTSurf), soil_con->GLAC_SURF_THICK, 
                                                  soil_con->GLAC_SURF_WE, air_temp, 
                                                  Tgrnd, &advection, &deltaCC_glac, 
                                                  &grnd_flux, &latent_heat, 
                                                  &latent_heat_sub, NetLongGlac, 
                                                  &sensible_heat, 
                                                  &glacier->vapor_flux);
                return (error);
            }
        }
        if (glacier->surf_temp > -998 && glacier->surf_temp < 999) {
            Qnet = CalcGlacierEnergyBalance(glacier->surf_temp,
                                            delta_t, aero_resist,
                                            aero_resist_used, z2, Z0,
                                            density, vp, LongUnderIn, 
                                            Le, pressure,
                                            RainFall, NetShortGlac, vpd,
                                            wind, (*OldTSurf),
                                            soil_con->GLAC_SURF_THICK,
                                            soil_con->GLAC_SURF_WE,
                                            air_temp, Tgrnd,
                                            &advection,
                                            &deltaCC_glac,
                                            &grnd_flux, &latent_heat,
                                            &latent_heat_sub, NetLongGlac,
                                            &sensible_heat,
                                            &glacier->vapor_flux);


            /* since we iterated, the surface layer is below freezing and no snowmelt */
            GlacMelt = 0.0;
            GlacCC = CONST_VCPICE_WQ * glacier->surf_temp * soil_con->GLAC_SURF_THICK / 1000.;
        }
    }

    melt[0] = GlacMelt;

    glacier->cold_content = GlacCC;
    glacier->vapor_flux *= -1.;
    *save_advection = advection;
    *save_deltaCC_glac = deltaCC_glac;
    *save_glacier_melt_energy = melt_energy;
    *save_grnd_flux = grnd_flux;
    *save_latent = latent_heat;
    *save_latent_sub = latent_heat_sub;
    *save_sensible = sensible_heat;
    *save_Qnet = Qnet;

    return (0);
}

/******************************************************************************
 * @brief    Dummy function to make a direct call to SnowEnergyBalance()
 *           possible.
 *****************************************************************************/
double
CalcGlacierEnergyBalance(double Tsurf,
                          ...)
{
    va_list ap;                 /* Used in traversing variable argument list
                                 */
    double  Qnet;                /* Net energy exchange at the SnowPack snow
                                    surface (W/m^2) */

    va_start(ap, Tsurf);
    Qnet = GlacierEnergyBalance(Tsurf, ap);
    va_end(ap);

    return Qnet;
}

/******************************************************************************
 * @brief    Pass snow pack energy balance terms to
 *           ErrorPrintSnowPackEnergyBalance
 *****************************************************************************/
int
ErrorGlacierEnergyBalance(double Tsurf,
                           ...)
{
    va_list ap;                 /* Used in traversing variable argument list
                                 */
    int     error;             /* error from ErrorPrintSnowPackEnergyBalance */

    va_start(ap, Tsurf);
    error = ErrorPrintGlacierEnergyBalance(Tsurf, ap);
    va_end(ap);

    return (error);
}

/******************************************************************************
 * @brief    Print snow pack energy balance terms
 *****************************************************************************/
int
ErrorPrintGlacierEnergyBalance(double  TSurf,
                               va_list ap)
{

    /* General Model Parameters */
    double  Dt;                      /* Model time step (sec) */
    double  Ra;                      /* Aerodynamic resistance (s/m) */
    double *Ra_used;                 /* Aerodynamic resistance (s/m) after stability correction */

    /* Vegetation Parameters */
 //   double Displacement,            /* Displacement height (m) */
    double  Z;                       /* Reference height (m) */
    double *Z0;                      /* surface roughness height (m) */

    /* Atmospheric Forcing Variables */
    double  AirDens;                 /* Density of air (kg/m3) */
    double  EactAir;                 /* Actual vapor pressure of air (Pa) */
    double  LongSnowIn;              /* Incoming longwave radiation (W/m2) */
    double  Lv;                      /* Latent heat of vaporization (J/kg3) */
    double  Press;                   /* Air pressure (Pa) */
    double  Rain;                    /* Rain fall (m/timestep) */
    double  NetShortUnder;           /* Net incident shortwave radiation (W/m2) */
    double  Vpd;                     /* Vapor pressure deficit (Pa) */
    double  Wind;                    /* Wind speed (m/s) */

    /* Snowpack Variables */
    double  OldTSurf;                /* Surface temperature during previous time step */
    double  IceDepth;                /* Depth of glacier surface layer (m) */
    double  IceWE;                   /* Liquid water in the glacier surface layer (m) */

    /* Energy Balance Components */
    double  Tair;                    /* Canopy air / Air temperature (C) */
    double  TGrnd;                   /* Ground surface temperature (C) */

    double *AdvectedEnergy;         /* Energy advected by precipitation (W/m2) */
    double *DeltaColdContent;       /* Change in cold content of surface layer (W/m2) */
    double *GroundFlux;             /* Ground Heat Flux (W/m2) */
    double *LatentHeat;             /* Latent heat exchange at surface (W/m2) */
    double *LatentHeatSub;          /* Latent heat of sublimation exchange at surface (W/m2) */
    double *NetLongUnder;           /* Net longwave radiation at snowpack surface (W/m^2) */
    double *SensibleHeat;           /* Sensible heat exchange at surface (W/m2) */
    double *VaporMassFlux;             /* Mass flux of water vapor to or from the intercepted snow (m/timestep) */

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
    VaporMassFlux = (double *) va_arg(ap, double *);

    /* print variables */
    log_warn("fglacier_melt failed to converge to a solution in "
             "root_brent.  Variable values will be dumped to the "
             "screen, check for invalid values.");

    /* general model terms */
    fprintf(LOG_DEST, "Dt = %f\n", Dt);

    /* land surface parameters */
    fprintf(LOG_DEST, "Ra = %f\n", Ra);
    fprintf(LOG_DEST, "Ra_used = %f\n", *Ra_used);
    fprintf(LOG_DEST, "Z = %f\n", Z);
    fprintf(LOG_DEST, "Z0 = %f\n", Z0);

    /* meteorological terms */
    fprintf(LOG_DEST, "AirDens = %f\n", AirDens);
    fprintf(LOG_DEST, "EactAir = %f\n", EactAir);
    fprintf(LOG_DEST, "LongSnowIn = %f\n", LongSnowIn);
    fprintf(LOG_DEST, "Lv = %f\n", Lv);
    fprintf(LOG_DEST, "Press = %f\n", Press);
    fprintf(LOG_DEST, "Rain = %f\n", Rain);
    fprintf(LOG_DEST, "NetShortUnder = %f\n", NetShortUnder);
    fprintf(LOG_DEST, "Vpd = %f\n", Vpd);
    fprintf(LOG_DEST, "Wind = %f\n", Wind);

    /* glacier terms */
    fprintf(LOG_DEST, "OldTSurf = %f\n", OldTSurf);
    fprintf(LOG_DEST, "IceDepth = %f\n", IceDepth);
    fprintf(LOG_DEST, "IceWE = %f\n", IceWE);
    fprintf(LOG_DEST, "Tair = %f\n", Tair);
    fprintf(LOG_DEST, "TGrnd = %f\n", TGrnd);
    fprintf(LOG_DEST, "AdvectedEnergy = %f\n", AdvectedEnergy[0]);
    fprintf(LOG_DEST, "DeltaColdContent = %f\n", DeltaColdContent[0]);
    fprintf(LOG_DEST, "GroundFlux = %f\n", GroundFlux[0]);
    fprintf(LOG_DEST, "LatentHeat = %f\n", LatentHeat[0]);
    fprintf(LOG_DEST, "LatentHeatSub = %f\n", LatentHeatSub[0]);
    fprintf(LOG_DEST, "NetLong = %f\n",NetLongUnder[0]);
    fprintf(LOG_DEST, "SensibleHeat = %f\n", SensibleHeat[0]);
    fprintf(LOG_DEST, "VaporMassFlux = %f\n", VaporMassFlux[0]);

    log_warn("Finished dumping glacier_melt variables.\nTry increasing "
             "SNOW_DT to get model to complete cell.\nThencheck output "
             "for instabilities.");

    return(ERROR);

}
