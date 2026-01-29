/******************************************************************************
 * @section DESCRIPTION
 *
 * Compute spatial average water table position (zwt).
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    Compute spatial average water table position (zwt).  Water table
 *           position is measured in cm and is negative below the soil surface.
 *****************************************************************************/
double
compute_zwt(double           dt,
            double           max_frost,
            cell_data_struct *cell,
            soil_con_struct  *soil_con)
{
    extern option_struct options;
    extern parameters_struct param;

    size_t  i, layer;
    double  decay_fator;
    double  discharge;
    double  coef_baseflow;
    double  sat_degree;
    double  psi_freeze;
    double  aquifer_conductivity;
    double  waterhead;
    double  recharge;
    double  storage_soil;
    double  storage_aqf;
    double  waterfill;
    double  WatConduct;
    double  min_moist;
    double  soil_excess;
    double  tmp_liq[MAX_NODES];
    double  tmp_conductivity[MAX_NODES];
    double *moist = cell->moist;
    double *zc_soil = soil_con->zc_soil;
    double *dz_soil = soil_con->dz_soil;
    double *Zsum_soil = soil_con->Zsum_soil;
    double *ice = cell->ice;
    double *liq = cell->liq;
    double *conductivity = soil_con->Ksat_node;
    double *Wsat_node = soil_con->Wsat_node;
    double *porosity = soil_con->porosity_node;

    size_t Nnode = options.Nnode;
    double zwt = cell->zwt;

    for (i = 0; i < Nnode; i++) {
        moist[i] = ice[i] + liq[i];
        tmp_liq[i] = liq[i] * dz_soil[i] * MM_PER_M;
        porosity[i] = max(0.01, Wsat_node[i] - ice[i]);
        tmp_conductivity[i] = conductivity[i] * MM_PER_M;  // convert from m/s to mm/s
    }

    layer = Nnode - 1;
    for (i = 1; i < Nnode; i++) {
        if (zwt <= Zsum_soil[i + 1]) {
            layer = i - 1;
            break;
        }
    }
    decay_fator = soil_con->bexp_node[layer] / 3.0;
    coef_baseflow = tmp_conductivity[layer] * 1.0e3 * exp(3.0);
    discharge = (1.0 - max_frost) * coef_baseflow * exp(-param.TOPOINDEX) * exp(-decay_fator * zwt);
    sat_degree = min(1.0, moist[layer] / Wsat_node[layer]);
    sat_degree = max(sat_degree, 0.01);
    if (sat_degree < 0.01) {
        sat_degree = 0.01;
    }
    psi_freeze = -soil_con->psi_sat_node[layer] * MM_PER_M * pow(sat_degree, soil_con->bexp_node[layer]);
    psi_freeze = max(psi_freeze * 0.2, -120000.0);  // cm
    aquifer_conductivity = 2.0 * (tmp_conductivity[layer] * conductivity[layer] * MM_PER_M) / 
                                    (tmp_conductivity[layer] + conductivity[layer] * MM_PER_M); // mm/s
    waterhead = psi_freeze - zc_soil[layer] * MM_PER_M;
    recharge = -aquifer_conductivity * (-zwt * MM_PER_M - waterhead) / ((zwt - zc_soil[layer]) * MM_PER_M);
    recharge = max(-10.0 / dt, min(10.0 / dt, recharge));
    storage_soil = (recharge - discharge) * dt; // mm
    if (layer == options.Nnode - 1) {
        storage_aqf += (recharge - discharge) * dt;
        storage_soil = storage_aqf;
        zwt = (Zsum_soil[layer] + 25.0) - storage_aqf / MM_PER_M / 0.2;
        liq[Nnode] -= recharge * dt;
        liq[Nnode] = max(0.0, storage_aqf - 5000.0);
        storage_aqf = min(5000.0, storage_aqf);
    }
    else {
        if (layer == Nnode - 2) {
            zwt = Zsum_soil[layer] - (storage_soil - 0.2 * MM_PER_M * 25.0) / porosity[layer + 1] / MM_PER_M;
        }
        else {
            waterfill = 0.0;
            for (i = layer + 2; i < Nnode; i++) {
                waterfill += porosity[i] * dz_soil[i] * MM_PER_M;
            }
            zwt = Zsum_soil[layer + 1] - (storage_soil - waterfill - 0.2 * MM_PER_M * 25.0) / porosity[layer + 1] / MM_PER_M;
        }
        WatConduct = 0.0;
        for (i = 0; i < Nnode; i++) {
            WatConduct += tmp_conductivity[i] * dz_soil[i] * MM_PER_M;
        }
        for (i = 0; i < Nnode; i++) {
            tmp_liq[i] -= discharge * dt * tmp_conductivity[i] * dz_soil[i] * MM_PER_M / WatConduct;
        }
    }
    cell->zwt = min(1.5, zwt);

    min_moist = 0.01;
    for (i = 0; i < Nnode - 1; i++) {
        if (tmp_liq[i] < 0.0) {
            soil_excess = min_moist - tmp_liq[i];
        }
        else {
            soil_excess = 0.0;
        }
        tmp_liq[i] += soil_excess;
        tmp_liq[i + 1] -= soil_excess;
    }

    if (tmp_liq[Nnode - 1] < min_moist) {
        soil_excess = min_moist - tmp_liq[Nnode - 1];
    }
    else {
        soil_excess = 0.0;
    }
    tmp_liq[Nnode - 1] += soil_excess;
    storage_aqf -= soil_excess;
    storage_soil -= soil_excess;

    for (i = 0; i < Nnode; i++) {
        liq[i] = tmp_liq[i] / (dz_soil[i] * MM_PER_M);
        liq[i] = max(liq[i], min_moist);
    }

    return(discharge);
}

/******************************************************************************
 * @brief    Function to compute spatial average water table position (zwt) for
 *           individual layers as well as various total-column versions of zwt.
 *           Water table position is measured in cm and is negative below the
 *           soil surface.
 *****************************************************************************/
void
wrap_compute_zwt(soil_con_struct  *soil_con,
                 cell_data_struct *cell)
{
    //
}
