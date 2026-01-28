/******************************************************************************
 * @section DESCRIPTION
 *
 * This set of functions calculate basic physical quantities that are
 * frequently calculated throughout the model
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    Compute outgoing longwave using the Stefan-Boltzman Law.
 *****************************************************************************/
double
calc_longwave(double  temp,
              double  coef_lw)
{
    double lwout;

    lwout = coef_lw * pow(temp, 4);

    return(lwout);
}

/******************************************************************************
 * @brief    Compute the sensible heat flux.
 *****************************************************************************/
double
calc_sensible_heat(double  Tgrnd,
                   double  Tair,
                   double  coef_sensible)
{
    double sensible;

    sensible = coef_sensible * (Tgrnd - Tair);

    return(sensible);
}

/******************************************************************************
 * @brief    Computes the latent heat from the glacier ground.
 *****************************************************************************/
double
calc_latent_heat(double  esat,
                 double  rel_humid,
                 double  vp,
                 double  coef_latent)
{
    double LatentHeat;

    LatentHeat = coef_latent * (esat * rel_humid - vp);

    return(LatentHeat);
}

/******************************************************************************
 * @brief    Computes the latent heat from the glacier ground.
 *****************************************************************************/
double
calc_ground_heat(double  TGrnd,
                 double  pack_temp,
                 double  coef_ground)

{
    double GroundFlux;

    GroundFlux = coef_ground * (TGrnd - pack_temp);

    return(GroundFlux);
}

