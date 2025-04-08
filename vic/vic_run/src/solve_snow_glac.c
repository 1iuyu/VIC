/******************************************************************************
* @section DESCRIPTION
*
* This routine was written to handle the various calls and data
* handling needed to solve the various components of the new VIC
* snow code for both the full_energy and water_balance models.
******************************************************************************/

#include <vic_run.h>

/******************************************************************************
* @brief        This routine was written to handle the various calls and data
*               handling needed to solve the various components of the new VIC
*               snow code for both the full_energy and water_balance models.
******************************************************************************/
double
solve_snow_glac(double             BareAlbedo,
                double             Tgrnd,          // soil surface temperature
                double             air_temp,          // air temperature
                double            *AlbedoUnder,
                double            *Le,
                double            *LongUnderIn,          // surface incomgin LW
                double            *NetLongGlac,          // net LW at glacier surface
                double            *NetShortGlac,         // net SW at glacier surface
                double            *ShortUnderIn,          // surfave incoming SW
                double            *Torg_snow,
                double            *aero_resist,
                double            *aero_resist_used,
                double            *coverage,          // best guess snow coverage
                double            *delta_coverage,          // cover fract change
                double            *displacement,
                double            *melt_energy,
                double            *ppt,
                double            *rainfall,
                double            *ref_height,
                double            *roughness,
                double            *snow_inflow,
                double            *snowfall,
                double            *wind,
                unsigned short     iveg,
                unsigned short     band,
                double             dt,
                size_t             hidx,
                int               *UnderStory,
                dmy_struct        *dmy,
                force_data_struct *force,
                energy_bal_struct *energy,
                snow_data_struct  *snow,
                soil_con_struct   *soil_con,
                glac_data_struct  *glacier)
{
    extern option_struct     options;
    extern parameters_struct param;

    int                      ErrorFlag;
    double                   ShortOverIn;
    double                   melt;
    double                   old_coverage;
    double                   old_depth;
    double                   old_swq;
    double                   tmp_grnd_flux;
//    double                   store_snowfall;
    int                      month;
    int                      day_in_year;
    double                   density;
    double                   longwave;
    double                   pressure;
    double                   shortwave;
    double                   vp;
    double                   vpd;

    month = dmy->month;
    day_in_year = dmy->day_in_year;

    density = force->density[hidx];
    longwave = force->longwave[hidx];
    pressure = force->pressure[hidx];
    shortwave = force->shortwave[hidx];
    vp = force->vp[hidx];
    vpd = force->vpd[hidx];
    /* initialize moisture variables */
    melt = 0.;
    *ppt = 0.;

    /* initialize storage for energy consumed in changing snowpack
       cover fraction */
    (*melt_energy) = 0.;

    /** Compute latent heats **/
    (*Le) = calc_latent_heat_of_vaporization(air_temp);

    /* initialize understory radiation inputs */
    (*ShortUnderIn) = shortwave;
    (*LongUnderIn) = longwave;

    snow->snow = true; // snow is present during time step

    old_coverage = snow->coverage; // store previous coverage fraction
    energy->NetLongOver = 0;
    energy->LongOverIn = 0;
    (*snow_inflow) = (*rainfall) + (*snowfall);
    old_swq = snow->swq;

    (*UnderStory) = 2;   /* ground snow is present of accumulating
                                    during time step */

    if (options.SPATIAL_SNOW) {
        /* make snowpack uniform at mean depth */
        if (*snowfall > 0) {
            snow->coverage = 1;
        }
        if (snow->coverage > 0 && *snowfall == 0) {
            if (snow->coverage < 1) {
                /* rain falls evenly over grid cell */
                *ppt = *rainfall * (1.0 - snow->coverage);
                *rainfall *= snow->coverage;
            }
        }
    }

    /** compute understory albedo and net shortwave radiation **/
    if (snow->swq > 0 && snowfall == 0) {
        // age snow albedo if no new snowfall
        // ignore effects of snow dropping from canopy; only consider fresh snow from sky
        snow->last_snow++;
        snow->albedo = snow_albedo(*snowfall, 
                                    soil_con->NEW_SNOW_ALB,
                                    snow->swq,
                                    snow->albedo,
                                    snow->coldcontent, dt,
                                    snow->last_snow, snow->MELTING);
        (*AlbedoUnder) =
                (*coverage * snow->albedo + (1. - *coverage) * BareAlbedo);
        }
        else {
            // set snow albedo to new snow albedo
            snow->last_snow = 0;
            snow->albedo = soil_con->NEW_SNOW_ALB;
            (*AlbedoUnder) = snow->albedo;
        }
        (*NetShortGlac) = (1.0 - *AlbedoUnder) * (*ShortUnderIn);

        /** Call snow pack accumulation and ablation algorithm **/
        ErrorFlag = snow_melt_glac((*Le), (*NetShortGlac), Tgrnd,
                                  roughness, aero_resist[*UnderStory],
                                  aero_resist_used,
                                  air_temp, *coverage, dt,
                                  density, displacement,
                                  *LongUnderIn, pressure, *rainfall, *snowfall,
                                  vp, vpd, wind[*UnderStory],
                                  ref_height[*UnderStory],
                                  NetLongGlac, Torg_snow, &melt, &energy->error,
                                  &energy->advected_sensible,
                                  &energy->advection,
                                  &energy->deltaCC, &energy->grnd_flux,
                                  &energy->latent,
                                  &energy->latent_sub, &energy->refreeze_energy,
                                  &energy->sensible, iveg, band,
                                  snow, glacier);
        if (ErrorFlag == ERROR) {
            return (ERROR);
        }

        // store melt water
        *ppt += melt;

        // store snow albedo
        energy->AlbedoUnder = *AlbedoUnder;

        /** Compute Snow Parameters **/
        if (snow->swq > 0.) {
            /** Calculate Snow Density **/
            if (snow->surf_temp <= 0) {
                // snowpack present, compress and age density
                snow->density = snow_density(snow, *snowfall, old_swq,
                                                 air_temp, dt);
            }
            else
            // no snowpack present, start with new snow density
            if (snow->last_snow == 0) {
                snow->density = new_snow_density(air_temp);
            }

            /** Calculate Snow Depth (H.B.H. 7.2.1) **/
            old_depth = snow->depth;
            snow->depth = CONST_RHOFW * snow->swq / snow->density;

            /** Record if snowpack is melting this time step **/
            if (snow->coldcontent >= 0 && (
                    (soil_con->lat >= 0 && (day_in_year > 60 && // ~ March 1
                                                day_in_year < 273)) || // ~ October 1
                    (soil_con->lat < 0 && (day_in_year < 60 || // ~ March 1
                                               day_in_year > 273)) // ~ October 1
                    )) {
                snow->MELTING = true;
            }
            else if (snow->MELTING && *snowfall > param.SNOW_TRACESNOW) {
                snow->MELTING = false;
            }


            /** Check for Thin Snowpack which only Partially Covers Grid Cell
                exists only if not snowing and snowpack has started to melt **/
            if (options.SPATIAL_SNOW) {
                snow->coverage = calc_snow_coverage(&snow->store_snow,
                                                    soil_con->max_snow_distrib_slope,
                                                    old_coverage, snow->swq,
                                                    old_swq, snow->depth, old_depth,
                                                    melt / MM_PER_M + snow->vapor_flux,
                                                    &snow->max_snow_depth, *snowfall,
                                                    &snow->store_swq,
                                                    &snow->snow_distrib_slope,
                                                    &snow->store_coverage);
            }
            else {
                if (snow->swq > 0) {
                    snow->coverage = 1.;
                }
                else {
                    snow->coverage = 0.;
                }
            }
        }
        else {
            snow->coverage = 0.;
        }

        *delta_coverage = old_coverage - snow->coverage;

        if (*delta_coverage != 0) {
            /* returns mixed surface albedo if snow cover fraction has
                decreased (old_coverage is cover fraction for previous
                time step, snow->coverage is cover fraction for current
                time step. */
            if (old_coverage > snow->coverage) {
                /* melt has occured */
                *coverage = (old_coverage);
                (*AlbedoUnder) = (*coverage - snow->coverage) /
                                    (1. - snow->coverage) * snow->albedo;
                (*AlbedoUnder) += (1. - *coverage) /
                                    (1. - snow->coverage) * BareAlbedo;

                /* compute snowpack energy used in reducing coverage area */
                (*melt_energy) = (*delta_coverage) *
                                    (energy->advection - energy->deltaCC +
                                    energy->latent + energy->latent_sub +
                                    energy->sensible +
                                    energy->refreeze_energy +
                                    energy->advected_sensible);
            }
            else if (old_coverage < snow->coverage) {
                *coverage = snow->coverage;
                *delta_coverage = 0;
            }
            else {
                *coverage = snow->coverage;
                *delta_coverage = 0.;
            }
        }
        else if (old_coverage == 0 && snow->coverage == 0) {
            // snow falls and melts all in one time step
            *delta_coverage = 1.;
            *coverage = 0.;
            (*melt_energy) = (energy->advection - energy->deltaCC +
                                energy->latent + energy->latent_sub +
                                energy->sensible + energy->refreeze_energy +
                                energy->advected_sensible);
        }

        /** Compute energy balance components for snowpack */

        (*NetLongGlac) *= (snow->coverage);
        (*NetShortGlac) *= (snow->coverage);
        energy->latent *= (snow->coverage + *delta_coverage);
        energy->latent_sub *= (snow->coverage + *delta_coverage);
        energy->sensible *= (snow->coverage + *delta_coverage);

        if (snow->swq == 0) {
            /** Reset Snow Pack Variables after Complete Melt **/

            /*** NOTE *coverage should not be zero the time step the
                     snowpack melts - FIX THIS ***/

            snow->density = 0.;
            snow->depth = 0.;
            snow->surf_water = 0;
            snow->pack_water = 0;
            snow->surf_temp = 0;
            snow->pack_temp = 0;
            snow->coverage = 0;
            snow->snow_distrib_slope = 0;
            snow->store_snow = true;
            snow->MELTING = false;
        }

        *snowfall = 0; /* all falling snow has been added to the pack */
        *rainfall = 0; /* all rain has been added to the pack */

        energy->melt_energy *= -1.;

        return(melt);
    }
