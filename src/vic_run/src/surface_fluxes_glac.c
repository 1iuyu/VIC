/******************************************************************************
* @section DESCRIPTION
*
* This routine computes all surface fluxes, and solves the snow accumulation
* and ablation algorithm. Solutions are for the current snow band and
* vegetation type.
******************************************************************************/

#include "vic_run.h"

/******************************************************************************
* @brief  This is a modified version of surface_fluxes.c specific for glaciers.
******************************************************************************/
int
surface_fluxes_glac(size_t             hidx,
                    double             step_dt,
                    double             air_temp,
                    double             snowfall,
                    double             rainfall,
                    force_data_struct *force,
                    energy_bal_struct *energy,
                    cell_data_struct  *cell,
                    snow_data_struct  *snow,
                    soil_con_struct   *soil_con)
{
    extern parameters_struct param;
    extern option_struct     options;
    int    ErrorFlag;
    size_t i;
    double shortwave_dir[MAX_SWBANDS];
    double shortwave_dfs[MAX_SWBANDS];

    double coverage = snow->coverage;
    double *AlbedoGrndDir = energy->AlbedoGrndDir;
    double *AlbedoGrndDfs = energy->AlbedoGrndDfs;

    /***************************
      Surface radiative fluxes
    ***************************/
    double NetShortGrnd = 0.0;
    double ReflShortSurf = 0.0;
    for (i = 0; i < options.Nswband; i++) {

        ReflShortSurf +=
            shortwave_dir[i] * AlbedoGrndDir[i] +
                shortwave_dfs[i] * AlbedoGrndDfs[i];
        energy->NetShortGrnd = NetShortGrnd;
        energy->ReflShortSurf = ReflShortSurf;
        energy->NetShortSurf = NetShortGrnd;
    }

    /******************************
      Compute longwave emissivity
    ******************************/

    double Ra_evap = 1.0;
    double rh_grnd = 1.0;
    
    /******************************
      Set psychrometric variables
    ******************************/
    energy->LatentVapGrnd = CONST_LATSUB;

    return(0);
}
