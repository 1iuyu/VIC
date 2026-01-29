/******************************************************************************
* @section DESCRIPTION
*
* Collection of snow utilities.
******************************************************************************/

#include <vic_run.h>

/******************************************************************************
* @brief     Estimate the density of new snow.
*
* @note     DENS_SNTHRM = Algorithm is taken from SNTHERM89, adjusted for an
*               essentially single-layer model.
*           DENS_BRAS   = Original algorithm, originally based on a plot of
*               seasonal variation of typical snow densities found in Bras
*               (figure 6.10, p 258).  Because this equation was developed by
*               regression against data from Southern Manitoba, it therefore is
*               limited in applicability to other regions.
******************************************************************************/
double
new_snow_density(double air_temp)

{
    extern parameters_struct param;
    extern option_struct     options;

    double                   density_new;

    density_new = 0.0;
    air_temp = K_TO_C(air_temp);
    
    if (options.SNOW_DENSITY == DENS_SNTHRM) {
        // new snow density based on Hedstrom and Pomeroy (1998)
        density_new = param.SNOW_NEW_SNT_C1 + param.SNOW_NEW_SNT_C2 * exp(
            air_temp / param.SNOW_NEW_SNT_C3);
    }
    else if (options.SNOW_DENSITY == DENS_BRAS) {
        // equation 6.2 in Bras 1990
        air_temp = C_TO_F(air_temp);
        if (air_temp > 0) {
            density_new = param.SNOW_NEW_SNOW_DENSITY + GRAMS_PER_KG *
                          (air_temp / param.SNOW_NEW_BRAS_DENOM) *
                          (air_temp / param.SNOW_NEW_BRAS_DENOM);
        }
        else {
            density_new = param.SNOW_NEW_SNOW_DENSITY;
        }
    }
    else {
        log_err("Unknown SNOW_DENSITY option");
    }

    // cap new snow density to prevent the calculation from
    // becoming unphysical
    if (density_new > param.SNOW_NEW_SNOW_DENS_MAX) {
        density_new = param.SNOW_NEW_SNOW_DENS_MAX;
    }

    return (density_new);
}

/******************************************************************************
* @brief    This subroutine computes the snow pack surface albedo.'
*
* @note     Computes albedo as a function of snow age and season, based on the
*           algorithm of the US Army Corps of Engineers.
******************************************************************************/
void
snow_albedo(double      coszen,
            double      SnowAge_fact,
            double     *AlbedoSnowDir,
            double     *AlbedoSnowDfs)
{
    extern parameters_struct param;

    double      solar_fact;
    double      coszen_fact;
    
    /** New Snow **/
    solar_fact = param.SNOW_COSZEN_B;
    coszen_fact = (1. + 1. / solar_fact) / (1. + 2 * solar_fact * coszen) - 
                                                             (1. / solar_fact);
    if (coszen_fact < 0.) {
        coszen_fact = 0.;
    }
    AlbedoSnowDfs[0] = param.SNOW_NEW_SNOW_VIS * 
                                 (1. - param.SNOW_AGE_DFS_VIS * SnowAge_fact);
    AlbedoSnowDfs[1] = param.SNOW_NEW_SNOW_NIR * 
                                 (1. - param.SNOW_AGE_DFS_NIR * SnowAge_fact);
    AlbedoSnowDir[0] = AlbedoSnowDfs[0] + param.SNOW_AGE_DIR_VIS * 
                                        coszen_fact * (1. - AlbedoSnowDfs[0]);
    AlbedoSnowDir[1] = AlbedoSnowDfs[1] + param.SNOW_AGE_DIR_NIR * 
                                        coszen_fact * (1. - AlbedoSnowDfs[1]);

}

/******************************************************************************
* @brief    This subroutine update snow water equivalent and snow depth
*
* @note     Computes albedo as a function of snow age and season, based on the
*           algorithm of the US Army Corps of Engineers.
******************************************************************************/
void
update_snow(double            air_temp,
            double            step_dt,
            double            snowfall,
            snow_data_struct *snow)
{

    /* 定义局部变量指向结构体成员 */
    double *dz_snow = snow->dz_snow;
    double *pack_T = snow->pack_T;
    double *pack_ice = snow->pack_ice;
    double *pack_liq = snow->pack_liq;

    /* initialize */
    size_t i = 0;
    size_t new_layer = 0;
    double delta_depth = snow->delta_depth;
    /* shallow snow or no layer */
    if (snow->Nsnow == 0 && snowfall > 0.) {
        snow->snow_depth += delta_depth * step_dt;
        snow->swq += snowfall * step_dt;
    }

    if (snow->Nsnow == 0 && snow->snow_depth > 0.025) {
        snow->Nsnow = 1;
        new_layer = 1;
        dz_snow[0] = snow->snow_depth;
        snow->snow_depth = 0.;
        pack_T[0] = min(CONST_TKTRIP, air_temp);
        pack_ice[0] = snow->swq;
        pack_liq[0] = 0.;
    }

    if (snow->Nsnow > 0 && new_layer == 0 && snowfall > 0.) {
        i = snow->Nsnow - 1;
        pack_ice[i] += snowfall * step_dt;
        dz_snow[i] += delta_depth * step_dt;
    }
}

/******************************************************************************
* @brief    This subroutine update snow water and temperature for 
*           combined snowpack layer
******************************************************************************/
void
update_surface_fluxes(double    *depth1,
                      double    *liq1,
                      double    *ice1,
                      double    *temp1,
                      double     depth2,
                      double     liq2,
                      double     ice2,
                      double     temp2)
{

    double  enthalpy1;
    double  enthalpy2;
    double  comb_depth;
    double  comb_ice;
    double  comb_liq;
    double  comb_enthalpy;
    double  comb_temp; 

    comb_depth = (*depth1) + depth2;
    comb_ice = (*ice1) + ice2;
    comb_liq = (*liq1) + liq2;
    enthalpy1 = (CONST_CPFWICE * (*liq1) + CONST_CPICE * (*ice1)) *
                (*temp1 - CONST_TKFRZ) + *liq1 * CONST_LATICE;
    enthalpy2 = (CONST_CPFWICE * liq2 + CONST_CPICE * ice2) *
                (temp2 - CONST_TKFRZ) + liq2 * CONST_LATICE;
    comb_enthalpy = enthalpy1 + enthalpy2;
    if (comb_enthalpy < 0.) {
        comb_temp = CONST_TKFRZ + comb_enthalpy / (CONST_CPFWICE * comb_liq +
                                            CONST_CPICE * comb_ice);
    }
    else if (comb_enthalpy < comb_ice * CONST_RHOFW) {
        comb_temp = CONST_TKFRZ;
    }
    else {
        comb_temp = CONST_TKFRZ + (comb_enthalpy - comb_liq * CONST_LATICE) / 
                        (CONST_CPFWICE * comb_liq + CONST_CPICE * comb_ice);
    }
    (*depth1) = comb_depth;
    (*liq1) = comb_liq;
    (*ice1) = comb_ice;
    (*temp1) = comb_temp;

}

/******************************************************************************
* @brief    This subroutine update snow water and temperature for 
*           combined snowpack layer
******************************************************************************/
void
distribute_snow_state(bool              IS_GLAC,
                      double            step_dt,
                      snow_data_struct *snow)
{
    extern option_struct    options;
    size_t i;
    double *dz_snow = snow->dz_snow;
    double *Zsum_snow = snow->Zsum_snow;
    double *zc_snow = snow->zc_snow;
    double *pack_T = snow->pack_T;
    double *pack_ice = snow->pack_ice;
    double *pack_liq = snow->pack_liq;
    if (snow->snow_depth <= 0.025) {
        snow->Nsnow = 0;
        pack_T[0] = min(250., CONST_TKFRZ);   // K
        pack_ice[0] = snow->swq;
        pack_liq[0] = 0.;
    }
    else if (snow->snow_depth <= 0.05) {
        snow->Nsnow = 1;
        dz_snow[0] = snow->snow_depth;
        zc_snow[0] = -dz_snow[0] / 2.0;
        Zsum_snow[0] = -dz_snow[0];
        pack_T[0] = min(250., CONST_TKFRZ);   // K
        pack_ice[0] = snow->swq;
        pack_liq[0] = 0.;
    }
    else if (snow->snow_depth <= 0.20) {
        snow->Nsnow = 2;
        dz_snow[0] = 0.05;
        dz_snow[1] = snow->snow_depth - 0.05;
        zc_snow[0] = -dz_snow[0] / 2.0;
        zc_snow[1] = -(dz_snow[0] + dz_snow[1] / 2.0);
        Zsum_snow[0] = -dz_snow[0];
        Zsum_snow[1] = -(dz_snow[0] + dz_snow[1]);
        pack_T[0] = min(250., CONST_TKFRZ);   // K
        pack_T[1] = min(250., CONST_TKFRZ);   // K
        pack_ice[0] = snow->swq * (dz_snow[0] / snow->snow_depth);
        pack_ice[1] = snow->swq * (dz_snow[1] / snow->snow_depth);
        pack_liq[0] = 0.;
        pack_liq[1] = 0.;

    }
    else {
        snow->Nsnow = 3;
        dz_snow[0] = 0.05;
        dz_snow[1] = 0.15;
        dz_snow[2] = snow->snow_depth - 0.20;
        zc_snow[0] = -dz_snow[0] / 2.0;
        zc_snow[1] = -(dz_snow[0] + dz_snow[1] / 2.0);
        zc_snow[2] = -(dz_snow[0] + dz_snow[1] + dz_snow[2] / 2.0);
        Zsum_snow[0] = -dz_snow[0];
        Zsum_snow[1] = -(dz_snow[0] + dz_snow[1]);
        Zsum_snow[2] = -(dz_snow[0] + dz_snow[1] + dz_snow[2]);
        pack_T[0] = min(250., CONST_TKFRZ);   // K
        pack_T[1] = min(250., CONST_TKFRZ);   // K
        pack_T[2] = min(250., CONST_TKFRZ);   // K
        pack_ice[0] = snow->swq * (dz_snow[0] / snow->snow_depth);
        pack_ice[1] = snow->swq * (dz_snow[1] / snow->snow_depth);
        pack_ice[2] = snow->swq * (dz_snow[2] / snow->snow_depth);
        pack_liq[0] = 0.;
        pack_liq[1] = 0.;
        pack_liq[2] = 0.;
    }
}
