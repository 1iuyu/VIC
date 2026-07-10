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
    size_t Nsoil = cell->Nsoil;
    size_t tmp_Nnode = soil_con->Nbedrock;
    size_t last_Nsnow = snow->last_Nsnow;
    double *T = energy->T;
    double *last_T = energy->last_T;
    double *Cs_node = energy->Cs_node;
    double *last_Cs = energy->last_Cs;
    double *ice = cell->ice;
    double *liq = cell->liq;
    double *last_ice = cell->last_ice;
    double *last_liq = cell->last_liq;
    double *matric = cell->matric;
    double *last_matric = cell->last_matric;
    double *pack_T = snow->pack_T;
    double *theta_ice = snow->theta_ice;
    double *theta_liq = snow->theta_liq;
    double *last_packice = snow->last_packice;
    double *last_packliq = snow->last_packliq;
    // update the number of nodes
    if (cell->h2osfc > 0.0) {
        tmp_Nnode++;
    }
    if (Nsnow > 0) {
        tmp_Nnode += Nsnow;
    }
    cell->Nnode = tmp_Nnode;
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
                lidx = i - Nsnow;
                new_T[i] = T[lidx];
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
    for (i = 0; i < tmp_Nnode; i++) {
        last_T[i] = T[i];
        last_Cs[i] = Cs_node[i];
    }
    for (i = 0; i < Nsoil; i++) {
        last_ice[i] = ice[i];
        last_liq[i] = liq[i];
        last_matric[i] = matric[i];
    }
    for (i = 0; i < Nsnow; i++) {
        last_packice[i] = theta_ice[i];
        last_packliq[i] = theta_liq[i];
    }
}