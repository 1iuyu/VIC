/******************************************************************************
 * @section DESCRIPTION
 *
 * Compute snow effective grain size (radius) based on SNICAR scheme 
 * (Flanner et al. (2021) GMD)
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    Determines from the air temperature what fraction of incoming
 *           precipitation is frozen and unfrozen (snow and rain).
 *****************************************************************************/
double
snow_aging(double            step_dt,
           double            swq,
           double           *old_swq,
           double            Tgrnd,
           double            snowfall,
           double           *SnowAge)
{
    extern option_struct     options;
    extern parameters_struct param;

    double      t_factor;
    double      snow_vapor;
    double      snow_refrz;
    double      snow_soot;
    double      delta_aging;
    double      new_snow_fact;
    double      tmp_SnowAge;
    double      SnowAge_fact;

    if (options.SNOW_AGING == SNICAR) {
/*      size_t      i;
        double      temp_upper;
        double      temp_lower;
        double      grad_temp;
        double      temp_index;
        double      grad_index;
        double      dens_index;
        double      tau;
        double      kappa;
        double      drdt0;
        double      delta_new_radius;
        double      part;
        double      exponent;
        double      delta_radius;
        double      water_frac;
        double      delta_wet_radius;
        double      new_snow;
        double      refrz_snow;
        double      new_snow_frac;
        double      refrz_frac;
        double      old_snow_frac;

        for (i = 0; i < options.Nsnows; i++) {
            snow_mass = snow->pack_liq[i] + snow->pack_ice[i];
            if (i == options.Nsnows) {
                temp_upper = snow->pack_temp[i];
                temp_lower = snow->pack_temp[i - 1] * snow->depth[i] +
                             snow->pack_temp[i] * snow->depth[i - 1] /
                             (snow->depth[i - 1] + snow->depth[i]);
            }
            else {
                temp_upper = snow->pack_temp[i + 1] * snow->depth[i] +
                             snow->pack_temp[i] * snow->depth[i + 1] /
                             (snow->depth[i] + snow->depth[i + 1]);
                temp_lower = snow->pack_temp[i - 1] * snow->depth[i] +
                             snow->pack_temp[i] * snow->depth[i - 1] /
                             (snow->depth[i - 1] + snow->depth[i]);
            }
            grad_temp[i] = fabs(temp_upper - temp_lower) / snow->depth[i];

            snow_density = snow_mass / snow->depth[i];
            if (snow_density < param.SNOW_NEW_SNOW_DENSITY) {
                snow_density = param.SNOW_NEW_SNOW_DENSITY;
            }

            temp_index = nint((snow->pack_temp[i] - 223.15) / 5) + 1.;
            grad_index = mint(grad_temp[i] / 10) + 1.;
            dens_index = mint((snow_density - 50.) / 50) + 1.;

            if (temp_index < param.SNOW_TEMP_MIN) {
                temp_index = param.SNOW_TEMP_MIN;
            }
            if (grad_index < param.SNOW_GRAD_MIN) {
                grad_index = param.SNOW_GRAD_MIN;
            }
            if (dens_index < param.SNOW_DENS_MIN) {
                dens_index = param.SNOW_DENS_MIN;
            }

            tau = tau_table[temp_index, grad_index, dens_index];
            kappa = kappa_table[temp_index, grad_index, dens_index];
            drdt = drdt_table[temp_index, grad_index, dens_index];

            delta_new_radius = snow->radius - param.SNOW_RADIUS_MIN;
            part = tau / (delta_new_radius + tau);

            if (part < 0.0) {
                part = 0.0; 
            }
            exponent = 1.0 / kappa;
            delta_radius = (drdt0 * pow(part, exponent));

            water_frac = min(0.1, (pack_water / snow_mass));
            delta_wet_radius = 1.0e18 * (step_dt * param.SNOW_WET_C1 + 
                               param.SNOW_WET_C2 * pow(water_frac, 3.0) /
                               (4.0 * CONST_PI * pow(snow->radius, 2.0)));
            delta_radius += delta_wet_radius;

            delta_radius *= param.SNOW_AGE_SCALE_F;

            new_snow = snowfall / MM_PER_M;
            refrz_snow = snowfreeze / MM_PER_M;

            refrz_frac = refrz_snow / snow_mass;

            if (lindex == options.Nsnows) {
                new_snow_frac = new_snow / snow_mass;
            }
            else {
                new_snow_frac = 0.;
            }

            if (new_snow_frac + refrz_frac > 1.) {
                refrz_frac = refrz_frac / (refrz_frac + new_snow_frac);
                new_snow_frac = 1.0 - refrz_frac;
                old_snow_frac = 0.0;
            }
            else {
                old_snow_frac = 1. - new_snow_frac - refrz_frac;
            }

            snow->radius = (snow->radius + delta_radius) * old_snow_frac +
                                        delta_new_radius * new_snow_frac + 
                                        delta_wet_radius * refrz_frac;


            if (snow->radius < param.SNOW_RADIUS_MIN) {
                snow->radius = param.SNOW_RADIUS_MIN;
            }
            if (snow->radius > param.SNOW_RADIUS_MAX) {
                snow->radius = param.SNOW_RADIUS_MAX;
            }
        }

        if (options.Nsnows == 0. && snowfall > 0. || 
                    snow->swq > 0. || snow->snow_depth > 0.) {
            snow->radius[0] = param.SNOW_NEW_RADIUS;
        }   */
        SnowAge_fact = 0.0;
        
        return (SnowAge_fact);
    }
    else if (options.SNOW_AGING == BATS) {

        SnowAge_fact = 0.0;
        if (swq <= 0.0) {
            *SnowAge = 0.;
        }
        else {
            /* No snow falling or present */
            t_factor = step_dt / param.SNOW_AGE_FACT;
            snow_vapor = exp(param.SNOW_AGE_VAPF * (1. / CONST_TKFRZ - 1. / Tgrnd));
            snow_refrz = exp(min(0., param.SNOW_AGE_FRZF * snow_vapor));
            snow_soot = param.SNOW_AGE_SOTF;
            delta_aging = t_factor * (snow_vapor + snow_refrz + snow_soot);
            new_snow_fact = max(0., swq - *old_swq) / param.SNOW_NEW_SNOW_COVER;
            tmp_SnowAge = (*SnowAge + delta_aging) * (1. - new_snow_fact);
            if (tmp_SnowAge < 0.) {
                *SnowAge = 0.;
            }
            else {
                *SnowAge = tmp_SnowAge;
            }
            SnowAge_fact = *SnowAge / (*SnowAge + 1.0);
        }
        *old_swq = swq;
        return (SnowAge_fact);
    }
}



