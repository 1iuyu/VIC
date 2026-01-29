/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate the aerodynamic resistances.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    Calculate the aerodynamic resistance for each vegetation layer,
 *           based on Monin-Obukhov (M-O) Similarity Theory (MOST).
 *****************************************************************************/
void
ActiveLayer(energy_bal_struct *energy,
            cell_data_struct  *cell,
            soil_con_struct   *soil_con) 
{
    extern option_struct options;

    size_t Nnode= options.Nnode;
    size_t idx = 0;
    size_t i;
    bool thawlayer = false;
    double *soil_T = cell->soil_T;
    double Nthaw = energy->Nthaw;
    double *Zsum_soil = soil_con->Zsum_soil;

    if (soil_T[Nnode-1] > CONST_TKFRZ) {
        Nthaw = Nnode - 1;
    }
    else {
        for (i = Nnode - 1; i >= 0; i--) {
            if (soil_T[i] > CONST_TKFRZ && thawlayer == false) {
               idx = i;
               thawlayer = true;
            }
        }
    }
    if (Nthaw > 0) {
        energy->tdepth = Zsum_soil[idx] + (soil_T[idx] - CONST_TKFRZ) * 
                            (Zsum_soil[idx+1] - Zsum_soil[idx]) / 
                                (soil_T[idx] - soil_T[idx + 1]);
        energy->Nthaw = Nthaw;
    }
    else {
        energy->tdepth = 0.0;
        energy->Nthaw = 0;
    }
}