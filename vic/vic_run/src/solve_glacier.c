/******************************************************************************
 * @section DESCRIPTION
 *  Calculate glacier accumulation and melt. Returns snow, veg_var, 
 *  and energy variables for each elevation band.  
 *  Variable ppt[] is defined for elevation bands with snow.
 *****************************************************************************/

#include <vic_run.h>

/*****************************************************************************
  * @brief   solve the various components of the new VIC glacier code for 
  *          both the full_energy and water_balance models.
 *****************************************************************************/
double 
solve_glacier(
              double               BareAlbedo,
              double               Tgrnd,               // glacier slab temperature
              double               air_temp,            // air temperature
              double              *AlbedoUnder,
              double              *Le,
              double              *LongUnderIn,         // surface incoming LW
              double              *NetLongGlac,         // net LW at glacier surface
              double              *NetShortGlac,        // net SW at glacier surface
              double              *ShortUnderIn,        // surface incoming SW
              double              *Torg_snow,
              double              *aero_resist,
              double              *aero_resist_used,
              double              *displacement,
              double              *melt_energy,
              double              *ppt,
              double              *rainfall,
              double              *ref_height,
              double              *Z0,
              double              *wind,
              int                  dt,
              int                  hidx,
              int                 *UnderStory,
              force_data_struct   *force,
              energy_bal_struct   *energy,
              glac_data_struct    *glacier,
              soil_con_struct     *soil_con,
              snow_data_struct    *snow) 
{

  int                 ErrorFlag;
//  double              ShortOverIn;
  double              melt;
//  double              tmp_grnd_flux;
//  double              store_snowfall;
//  double              tmp_ref_height;
  double              density;
  double              longwave;
  double              pressure;
  double              shortwave;
  double              vp;
  double              vpd;

  density   = force->density[hidx];
  longwave  = force->longwave[hidx];
  pressure  = force->pressure[hidx];
  shortwave = force->shortwave[hidx];
  vp        = force->vp[hidx];
  vpd       = force->vpd[hidx];

  /* initialize moisture variables */
  melt     = 0.;
  *ppt     = 0.;

  /* initialize storage for energy consumed in changing snowpack cover fraction */
  (*melt_energy)     = 0.;

  /** Compute latent heats **/
  (*Le) = calc_latent_heat_of_vaporization(air_temp);

  /* initialize glacier surface radiation inputs */
  (*ShortUnderIn) = shortwave;
  (*LongUnderIn)  = longwave;

  /** compute net shortwave radiation **/
  (*AlbedoUnder) = BareAlbedo;
  (*NetShortGlac) = (1.0 - *AlbedoUnder) * (*ShortUnderIn);

  /* glacier is present */
  *UnderStory = 0;
   /** Call glacier ablation algorithm **/
   ErrorFlag = glacier_melt((*Le), (*NetShortGlac), Tgrnd,
                            Z0, aero_resist[*UnderStory],
                            aero_resist_used,
                            air_temp, dt, density,
                            displacement[*UnderStory], *LongUnderIn,
                            pressure, *rainfall, vp, vpd, 
                            wind[*UnderStory], ref_height[*UnderStory],
                            NetLongGlac, Torg_snow, &melt,
                        	  &energy->error, &energy->advection,
                            &energy->deltaCC_glac,
                        	  &energy->glacier_melt_energy, &energy->grnd_flux,
                            &energy->latent, &energy->latent_sub,
                            &energy->sensible, glacier, soil_con, snow);

   if (ErrorFlag == ERROR) {
      return (ERROR);
   }
    

   // store melt water and rainfall
   *ppt = (melt + *rainfall / MM_PER_M); /* convert rainfall to m */

   // store glacier albedo
   energy->AlbedoUnder = *AlbedoUnder;

   *rainfall = 0; /* all rain has been added to the glacier */

  return(melt);

}
