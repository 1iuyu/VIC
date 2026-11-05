/******************************************************************************
 * @brief
 *
 * Update variables from current time step for the next implicit calculation.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    This routine should be called after a successful sub-step convergence.
 *****************************************************************************/
void
update_last_state(energy_bal_struct *energy,
                  cell_data_struct  *cell,
                  snow_data_struct  *snow)
{
    // initialize variables
    size_t i;
    size_t Nnode = cell->Nnode;
    size_t Nsoil = cell->Nsoil;
    size_t Nsnow = snow->Nsnow;
    double *T = energy->T;
    double *last_T = energy->last_T;
    double *ice = cell->ice;
    double *liq = cell->liq;
    double *matric = cell->matric;
    double *last_ice = cell->last_ice;
    double *last_liq = cell->last_liq;
    double *moist = cell->moist;
    double *last_moist = cell->last_moist;
    double *last_matric = cell->last_matric;
    double *theta_ice = snow->theta_ice;
    double *theta_liq = snow->theta_liq;
    double *enthalpy = energy->enthalpy;
    double *last_enthalpy= energy->last_enthalpy;
    double *last_thtice = snow->last_thtice;
    double *last_thtliq = snow->last_thtliq;
    double *theta_moist = snow->theta_moist;
    double *last_thtmoist = snow->last_thtmoist;
    // Set to the default value
    energy->energy_flag = false;
    energy->moist_flag = false;
    energy->Esignchg_count = 0;
    energy->Msignchg_count = 0;
    energy->energy_error = 0.0;
    energy->moist_error = 0.0;
    
    // Update thermal states
    for (i = 0; i < Nnode; i++) {
        last_T[i] = T[i];
        last_enthalpy[i] = enthalpy[i];
    }

    // Update soil water states
    for (i = 0; i < Nsoil; i++) {
        last_ice[i] = ice[i];
        last_liq[i] = liq[i];
        last_moist[i] = moist[i];
        last_matric[i] = matric[i];
    }

    // Update snow states
    for (i = 0; i < Nsnow; i++) {
        last_thtice[i] = theta_ice[i];
        last_thtliq[i] = theta_liq[i];
        last_thtmoist[i] = theta_moist[i];
    }
}