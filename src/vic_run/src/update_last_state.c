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
    size_t i;
    size_t Nnode = cell->Nnode;
    size_t Nsoil = cell->Nsoil;
    size_t Nsnow = snow->Nsnow;
    double *T = energy->T;
    double *Cs_node = energy->Cs_node;
    double *last_T = energy->last_T;
    double *last_Cs = energy->last_Cs;
    double *ice = cell->ice;
    double *liq = cell->liq;
    double *matric = cell->matric;
    double *last_ice = cell->last_ice;
    double *last_liq = cell->last_liq;
    double *last_matric = cell->last_matric;
    double *theta_ice = snow->theta_ice;
    double *theta_liq = snow->theta_liq;
    double *snow_frac = snow->snow_frac;
    double *last_packice = snow->last_packice;
    double *last_packliq = snow->last_packliq;
    double *last_snowfrac = snow->last_snowfrac;
    // 
    energy->energy_flag = false;
    energy->moist_flag = false;
    energy->Esignchg_count = 0;
    energy->Msignchg_count = 0;
    energy->energy_error = 0.0;
    energy->moist_error = 0.0;
    
    // Update thermal states
    for (i = 0; i < Nnode; i++) {
        last_T[i] = T[i];
        last_Cs[i] = Cs_node[i];
    }

    // Update soil water states
    for (i = 0; i < Nsoil; i++) {
        last_ice[i] = ice[i];
        last_liq[i] = liq[i];
        last_matric[i] = matric[i];
    }

    // Update snow states
    for (i = 0; i < Nsnow; i++) {
        last_packice[i] = theta_ice[i];
        last_packliq[i] = theta_liq[i];
        last_snowfrac[i] = snow_frac[i];
    }
}