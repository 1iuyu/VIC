/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine computes the state variables (energy balance, water balance,
 * and snow components) that are derived from the variables that are stored in
 * state files.
 *****************************************************************************/

#include "vic_driver_shared_all.h"

/******************************************************************************
 * @brief    Compute the state variables (energy balance, water balance,
 *           and snow components) that are derived from the variables that
 *           are stored in state files.
 *****************************************************************************/
void
compute_derived_state_vars(all_vars_struct *all_vars,
                           soil_con_struct *soil_con,
                           veg_con_struct  *veg_con)
{
    extern global_param_struct global_param;
    extern parameters_struct   param;
    extern option_struct       options;
    size_t veg, Nveg;
    size_t Nsnow;
    size_t i, lidx;
    double Cv;
    double dt_thresh;
    cell_data_struct *cell;
    snow_data_struct *snow;
    energy_bal_struct *energy;
    cell = all_vars->cell;
    snow = all_vars->snow;
    energy = all_vars->energy;
    Nveg = veg_con[0].vegetat_type_num;

    if (options.INIT_STATE) {
        // compute snow, soil and h2osfc_T temperatures
        for (veg = 0; veg <= Nveg; veg++) {
            if (veg_con[veg].Cv > 0) {
                Nsnow = snow[veg].Nsnow;
                // Initialize snow node temperatures
                for (i = 0; i < Nsnow; i++) {
                    snow[veg].pack_T[i] = energy[veg].T[i];
                }
                if (cell[veg].h2osfc > param.TOL_A) {
                    cell[veg].h2osfc_T = energy[veg].T[Nsnow];
                    Nsnow++;
                }
                /* Initialize soil node temperatures */
                for (i = Nsnow; i < cell[veg].Nnode; i++) {
                    lidx = i - Nsnow;
                    cell[veg].soil_T[lidx] = energy[veg].T[i];
                }
            }
        }
        // 计算累计厚度和中心位置
        for (veg = 0; veg <= Nveg; veg++) {
            double cum_depth = 0.0;
            for (i = 0; i < snow[veg].Nsnow; i++) {
                cum_depth += snow[veg].dz_snow[i];
                snow[veg].Zsum_snow[i] = cum_depth;
                snow[veg].zc_snow[i] = cum_depth - snow[veg].dz_snow[i] / 2.0;
                snow[veg].theta_ice[i] = min(1.0, snow[veg].pack_ice[i] / 
                            (snow[veg].dz_snow[i] * snow[veg].coverage * CONST_RHOICE));
                snow[veg].porosity[i] = 1.0 - snow[veg].theta_ice[i];
                snow[veg].theta_liq[i] = min(snow[veg].porosity[i], snow[veg].pack_liq[i] / 
                                (snow[veg].dz_snow[i] * snow[veg].coverage * CONST_RHOFW));
                double SnowMass = snow[veg].pack_ice[i] + snow[veg].pack_liq[i];
                snow[veg].snow_frac[i] = snow[veg].pack_ice[i] / SnowMass;
                snow[veg].density[i] = snow[veg].pack_ice[i] / 
                                      (snow[veg].dz_snow[i] * snow[veg].coverage);
                double CP_snow = 92.96 + 7.37 * snow[veg].pack_T[i];
                snow[veg].enthalpy[i] = (snow[veg].pack_ice[i] * CP_snow + 
                                snow[veg].pack_liq[i] * CONST_CPFWICE) * (snow[veg].pack_T[i] - 
                                CONST_TKFRZ) - snow[veg].pack_ice[i] * CONST_LATICE;
                snow[veg].enthalpy[i] /= snow[veg].coverage;
            }
        }
    }

    /******************************************
       Compute maximum daylight duration
    ******************************************/
    double max_daylen = calc_max_daylength(soil_con->lat);
    for (veg = 0; veg <= Nveg; veg++) {
        if (veg_con[veg].Cv > 0.0) {
            cell[veg].max_daylen = max_daylen;
        }
    }

    /******************************************
       Compute soil thermal node properties
    ******************************************/
    for (veg = 0; veg <= Nveg; veg++) {
        // Initialize soil for existing vegetation types
        Cv = veg_con[veg].Cv;
        if (Cv > 0) {
            // Check node spacing v time step
            // (note this is only approximate since heat capacity and
            // conductivity can change considerably during the
            // simulation depending on soil moisture and ice content)
            // in seconds
            dt_thresh = 0.5 * energy[veg].Cs_node[0] /
                        energy[veg].kappa_node[0] *
                        pow((soil_con->dz_soil[0]), 2.0);

            if (global_param.step_dt > dt_thresh) {
                log_warn("Crank-Nicolson stability condition is violated. "
                            "The current dimensionless parameter "
                            "exceeds the critical threshold of 0.5. "
                            "This may cause spurious oscillations or unphysical "
                            "decay in the numerical solution. "
                            "Now switch to the implicit method "
                            "Current settings: time step = %f s, "
                            "first layer node spacing = %f m, "
                            "soil heat capacity = %f J/m3/K, "
                            "soil thermal conductivity = %f J/m/s/K. "
                            "Recommended maximum stable time step: %f s. "
                            "To resolve this: Reduce the time step to <= %f s,",
                            global_param.step_dt,
                            soil_con->dz_soil[0],
                            energy[veg].Cs_node[0],
                            energy[veg].kappa_node[0],
                            dt_thresh, dt_thresh);
            }
        }
    }
}
