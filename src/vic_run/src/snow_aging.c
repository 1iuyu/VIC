/******************************************************************************
 * @section DESCRIPTION
 *
 * Compute snow effective grain size (radius) based on SNICAR scheme 
 * (Flanner et al. (2021) GMD)
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    Determines from the air temperature what fraction of incoming
 *           precipitation is frozen and unfrozen (snow and rain).
 *****************************************************************************/
int
snow_aging(double            step_dt,
           double            Tgrnd,
           double            air_temp,
           double            snowfall,
           snow_data_struct *snow)
{
    extern option_struct     options;
    extern parameters_struct param;

    if (options.SNOW_AGING == SNICAR) {
        size_t i, Nsnow;
        double snow_mass = 0.0;
        double grad_temp = 0.0;
        double snow_density = 0.0;
        double *tau_table;
        double *kappa_table;
        double *drdt_table;
        double *pack_T = snow->pack_T;
        double *radius = snow->radius;
        double *dz_snow = snow->dz_snow;
        double *pack_frze = snow->pack_frze;
        double *pack_liq = snow->pack_liq;
        double *pack_ice = snow->pack_ice;
        double temp_upper;
        double temp_lower;
        double grad_temp;
        double delta_new_radius;
        double part;
        double exponent;
        double delta_radius = 0.0;
        double frac_liq;
        double delta_wet_radius;
        double new_snow;
        double refrz_snow;
        double new_snow_frac;
        double refrz_frac;
        double old_snow_frac;
        Nsnow = snow->Nsnow;
        for (i = 0; i < Nsnow; i++) {
            snow_mass = pack_liq[i] + pack_ice[i];
            if (i == 0) {
                temp_upper = pack_T[i];
                temp_lower = pack_T[i+1] * dz_snow[i] + pack_T[i] * 
                        dz_snow[i+1] / (dz_snow[i+1] + dz_snow[i]);
            }
            else {
                temp_upper = pack_T[i-1] * dz_snow[i] + pack_T[i] * 
                            dz_snow[i-1] / (dz_snow[i] + dz_snow[i-1]);
                temp_lower = pack_T[i+1] * dz_snow[i] + pack_T[i] * 
                            dz_snow[i+1] / (dz_snow[i+1] + dz_snow[i]);
            }
            grad_temp = fabs(temp_upper - temp_lower) / dz_snow[i];

            snow_density = snow_mass / dz_snow[i];
            if (snow_density < param.SNOW_NEW_SNOW_DENSITY) {
                snow_density = param.SNOW_NEW_SNOW_DENSITY;
            }
            // best-fit table indices
            size_t temp_idx = nint((pack_T[i] - 223.15) / 5.0) + 1;
            size_t grad_idx = mint(grad_temp / 10.0) + 1;
            size_t dens_idx = mint((snow_density - 50.0) / 50.0) + 1;
            // boundary checks
            temp_idx = min(max(temp_idx, 1), 11);
            grad_idx = min(max(grad_idx, 1), 31);
            dens_idx = min(max(dens_idx, 1), 8);
            // best-fit parameters
            double best_tau = tau_table[temp_idx, grad_idx, dens_idx];
            double best_kappa = kappa_table[temp_idx, grad_idx, dens_idx];
            double best_drdt = drdt_table[temp_idx, grad_idx, dens_idx];
            radius[i] = max(radius[i], param.SNOW_RADIUS_MIN);
            delta_new_radius = radius[i] - param.SNOW_RADIUS_MIN;
            part = best_tau / (delta_new_radius + best_tau);
            if (part < 0.0) {
                part = 0.0; 
            }
            exponent = 1.0 / best_kappa;
            delta_radius = (best_drdt * pow(part, exponent)) * (step_dt / SEC_PER_HOUR);

            frac_liq = min(0.1, (pack_liq[i] / snow_mass));
            delta_wet_radius = 1.0e18 * (step_dt * param.SNOW_WET_C2 * pow(frac_liq, 3.0) /
                                                (4.0 * CONST_PI * pow(radius[i], 2.0)));
            delta_radius += delta_wet_radius;

            delta_radius *= param.SNOW_AGE_SCALE_F;
            // new snowfall [mm]
            new_snow = max(snowfall * step_dt, 0.0);
            // snow that has re-frozen [mm]
            refrz_snow = max(pack_frze[i] * step_dt, 0.0);
            refrz_frac = refrz_snow / snow_mass;
            // fraction of layer mass that is new snow
            if (i == 0) {
                new_snow_frac = new_snow / snow_mass;
            }
            else {
                new_snow_frac = 0.0;
            }

            if (new_snow_frac + refrz_frac > 1.0) {
                refrz_frac = refrz_frac / (refrz_frac + new_snow_frac);
                new_snow_frac = 1.0 - refrz_frac;
                old_snow_frac = 0.0;
            }
            else {
                old_snow_frac = 1.0 - new_snow_frac - refrz_frac;
            }
            // temperature dependent fresh grain size
            double fresh_radius = new_snow_radius(air_temp);

            radius[i] = (radius[i] + delta_radius) * old_snow_frac +
                    fresh_radius * new_snow_frac + param.SNOW_REFRZF * refrz_frac;

            if (radius[i] < param.SNOW_RADIUS_MIN) {
                radius[i] = param.SNOW_RADIUS_MIN;
            }
            if (radius[i] > param.SNOW_RADIUS_MAX) {
                radius[i] = param.SNOW_RADIUS_MAX;
            }
        }
        // set to fresh snow grain size
        if (Nsnow == 0.0 && snowfall > 0.0 || 
                    snow->swq > 0.0 || snow->snow_depth > 0.0) {
            radius[0] = param.SNOW_NEW_RADIUS;
        }
    }
    else if (options.SNOW_AGING == BATS) {
        // initialize
        double t_factor;
        double snow_vapor;
        double snow_refrz;
        double snow_soot;
        double delta_aging;
        double new_snow_fact;
        double tmp_SnowAge;
        double snowage = 0.0;
        double swq = snow->swq;
        double last_swq = snow->last_swq;

        if (swq <= 0.0) {
            snowage = 0.;
        }
        else {
            /* No snow falling or present */
            t_factor = step_dt / param.SNOW_AGE_FACT;
            snow_vapor = exp(param.SNOW_AGE_VAPF * (1. / CONST_TKFRZ - 1. / Tgrnd));
            snow_refrz = exp(min(0., param.SNOW_AGE_FRZF * snow_vapor));
            snow_soot = param.SNOW_AGE_SOTF;
            delta_aging = t_factor * (snow_vapor + snow_refrz + snow_soot);
            new_snow_fact = max(0., swq - last_swq) / param.SNOW_NEW_SNOW_COVER;
            tmp_SnowAge = (snowage + delta_aging) * (1. - new_snow_fact);
            if (tmp_SnowAge < 0.) {
                snowage = 0.;
            }
            else {
                snowage = tmp_SnowAge;
            }
        }
        snow->snowage = snowage;
    }

    return (0);
}