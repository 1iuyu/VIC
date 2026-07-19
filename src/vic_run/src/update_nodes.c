/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine update the nodes properties in the vertical direction.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    This routine update the nodes properties in the vertical direction.
 *****************************************************************************/
void
update_nodes(double             pressure,
             energy_bal_struct *energy,
             cell_data_struct  *cell,
             snow_data_struct  *snow,
             soil_con_struct   *soil_con)
{
    // initialize variables
    size_t i, lidx;
    size_t Nsnow = snow->Nsnow;
    size_t tmp_Nsnow = snow->Nsnow;
    size_t tmp_Nnode = soil_con->Nbedrock;
    size_t last_Nsnow = snow->last_Nsnow;
    double *T = energy->T;
    double *density = snow->density;
    double *last_T = energy->last_T;
    double *soil_T = cell->soil_T;
    double *Cs_node = energy->Cs_node;
    double *last_Cs = energy->last_Cs;
    double *pack_T = snow->pack_T;
    double *dz_snow = snow->dz_snow; 
    double *pack_ice = snow->pack_ice;
    double *pack_liq = snow->pack_liq;
    double *porosity = snow->porosity;
    double *snow_frac = snow->snow_frac;
    double *theta_ice = snow->theta_ice;
    double *theta_liq = snow->theta_liq;
    double *kappa_node = energy->kappa_node;
    double *kappa_int = energy->kappa_int;
    double *last_snowfrac = snow->last_snowfrac;
    double *pack_outflow = snow->pack_outflow;
    double *last_packice = snow->last_packice;
    double *last_packliq = snow->last_packliq;
    // update the number of nodes
    if (cell->h2osfc > 0.0) {
        tmp_Nnode++;
        tmp_Nsnow++;
    }
    if (Nsnow > 0) {
        tmp_Nnode += Nsnow;
    }
    // update snow layer properties
    cell->Nnode = tmp_Nnode;
    for (i = 0; i < Nsnow; i++) {
        theta_ice[i] = min(1.0, pack_ice[i] / (dz_snow[i] * CONST_RHOICE));
        porosity[i] = 1.0 - theta_ice[i];
        theta_liq[i] = min(porosity[i], pack_liq[i] / (dz_snow[i] * CONST_RHOFW));
        double SnowMass = pack_ice[i] + pack_liq[i];
        snow_frac[i] = pack_ice[i] / SnowMass;
        density[i] = pack_ice[i] / dz_snow[i];
    }

    // update new snow layer properties
    if (Nsnow != last_Nsnow) {
        double new_T[MAX_NODES];
        for (i = 0; i < MAX_NODES; i++) {
            new_T[i] = 0.0;
        }
        for (i = 0; i < cell->Nnode; i++) {
            if (i < Nsnow) {
                new_T[i] = pack_T[i];
            }
            else if (i == Nsnow && cell->h2osfc > 0.0) {
                new_T[i] = cell->h2osfc_T;
            }
            else {
                lidx = i - tmp_Nsnow;
                new_T[i] = soil_T[lidx];
            }
        }
        for (i = 0; i < cell->Nnode; i++) {
            T[i] = new_T[i];
        }
        // update the Cs_node for the new snow layer
        prepare_full_energy(pressure, 
                            cell, energy,
                            snow, soil_con);
    }

    // update the last time step values
    snow->last_Nsnow = snow->Nsnow;
    snow->last_swq = snow->swq;
    for (i = 0; i < tmp_Nnode; i++) {
        last_T[i] = T[i];
        last_Cs[i] = Cs_node[i];
    }
    // Update snow states
    for (i = 0; i < Nsnow; i++) {
        last_packice[i] = theta_ice[i];
        last_packliq[i] = theta_liq[i];
        last_snowfrac[i] = snow_frac[i];
    }

    /* remove old snow layers */
    for(i = Nsnow; i < last_Nsnow; i++) {
        density[i] = 0.0;
        porosity[i] = 0.0;
        theta_ice[i] = 0.0;
        theta_liq[i] = 0.0;
        last_packice[i] = 0.0;
        last_packliq[i] = 0.0;
        pack_outflow[i] = 0.0;
    }
    if (last_Nsnow > Nsnow) {
        lidx = last_Nsnow - Nsnow + tmp_Nnode;
        for (i = tmp_Nnode; i < lidx; i++) {
            T[i] = 0.0;
            last_T[i] = 0.0;
            last_Cs[i] = 0.0;
            Cs_node[i] = 0.0;
            kappa_int[i] = 0.0;
            kappa_node[i] = 0.0;
        }
    }
}