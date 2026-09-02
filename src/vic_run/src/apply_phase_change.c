/******************************************************************************
 * @section DESCRIPTION
 *   
 * This routine computes phase-change heat flux for next iteration.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    This routine computes phase-change heat flux for next iteration.
 *****************************************************************************/
void
apply_phase_change(double            step_dt,
                   snow_data_struct *snow)
{
    double phase_mass;
    double max_freeze;
    double max_melt;
    double coverage = snow->coverage;
    double *pack_ice = snow->pack_ice;
    double *pack_liq = snow->pack_liq;
    double *dz_snow = snow->dz_snow;
    double *porosity = snow->porosity;
    double *theta_ice = snow->theta_ice;
    double *theta_liq = snow->theta_liq;
    double *pack_frze = snow->pack_frze;
    double *pack_melt = snow->pack_melt;
    double *phase_snow = snow->phase_snow;

    for (size_t i = 0; i < snow->Nsnow; i++) {
        // Initialize phase-change fluxes for this snow layer
        pack_frze[i] = 0.0;
        pack_melt[i] = 0.0;
        /* No phase change */
        if (fabs(phase_snow[i]) <= 1.0e-12) {
            phase_snow[i] = 0.0;
            continue;
        }
        /* Convert phase-change energy to mass */
        phase_mass = fabs(phase_snow[i]) * step_dt * coverage / CONST_LATICE;
        /* Freezing: liquid water -> ice */
        if (phase_snow[i] > 0.0) {
            max_freeze = max(0.0, pack_liq[i]);
            phase_mass = min(phase_mass, max_freeze);
            pack_liq[i] -= phase_mass;
            pack_ice[i] += phase_mass;
            pack_frze[i] += phase_mass;
        }
        else if (phase_snow[i] < 0.0) {
            max_melt = max(0.0, pack_ice[i]);
            phase_mass = min(phase_mass, max_melt);
            pack_ice[i] -= phase_mass;
            pack_liq[i] += phase_mass;
            pack_melt[i] += phase_mass;
        }

        pack_ice[i] = max(0.0, pack_ice[i]);
        pack_liq[i] = max(0.0, pack_liq[i]);

        if (dz_snow[i] > 0.0 && coverage > 0.0) {
            theta_ice[i] = pack_ice[i] / (dz_snow[i] * coverage * CONST_RHOICE);
            theta_liq[i] = pack_liq[i] / (dz_snow[i] * coverage * CONST_RHOFW);
            theta_ice[i] = min(1.0, max(0.0, theta_ice[i]));
            theta_liq[i] = max(0.0, theta_liq[i]);
            porosity[i] = max(0.0, 1.0 - theta_ice[i]);
        }
        else {
            theta_ice[i] = 0.0;
            theta_liq[i] = 0.0;
            porosity[i] = 0.0;
        }
        // Clear phase-change flux after applying it.
        phase_snow[i] = 0.0;
    }
}