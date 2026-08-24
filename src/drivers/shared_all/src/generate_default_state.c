/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes the model state (energy balance, water balance, and
 * snow components) to default values.
 *****************************************************************************/

#include "vic_driver_shared_all.h"

/******************************************************************************
 * @brief    Initialize the model state (energy balance, water balance, and
 *           snow components) to default values.
 *****************************************************************************/
void
generate_default_state(force_data_struct *force,
                       all_vars_struct   *all_vars,
                       soil_con_struct   *soil_con,
                       veg_con_struct    *veg_con,
                       veg_lib_struct    *veg_lib)
{
    extern option_struct       options;
    extern global_param_struct global_param;
    int         ErrorFlag;
    size_t      Nveg;
    size_t      veg;
    size_t      band;
    size_t      Nsnow;
    size_t      lidx;
    size_t      i, k;
    size_t      veg_class;
    double      Canopy_Upper;
    cell_data_struct *cell;
    snow_data_struct *snow;
    veg_var_struct *veg_var;
    energy_bal_struct *energy;
    cell = all_vars->cell;
    snow = all_vars->snow;
    veg_var = all_vars->veg_var;
    energy = all_vars->energy;
    Nveg = veg_con[0].vegetat_type_num;
    double *zc_soil = soil_con->zc_soil;
    double step_dt = global_param.step_dt;
    double wind = force->wind[NR];
    double air_temp = force->air_temp[NR];
    double pressure = force->pressure[NR];

    /************************************
      Initialize snowpack properties
    ************************************/
    for (veg = 0; veg <= Nveg; veg++) {
        if (veg_con[veg].Cv > 0) {
            snow[veg].swq = 15.0;   // [mm] or [kg/m^2]
            double snow_density = 200.0;   // [kg/m^3]
            if (snow[veg].swq > 0.0) {
                // calculate snow coverage
                double GridSize = sqrt(soil_con->cell_area * veg_con[veg].Cv);
                double gridScalePara = min((max(GridSize, 500.0) / M_PER_KM), 36.0);
                double SNOW_MeltFac = 0.9713 + tanh(0.7436 * gridScalePara);
                double SNOW_CoverFac = 0.0062 * sinh(0.0555 * gridScalePara) + 0.0555;
                double MeltFac = pow(snow_density / 100, SNOW_MeltFac);
                // 迭代求解一个合理的coverage
                double coverage = 1.0;
                double A = snow[veg].swq / (snow_density * SNOW_CoverFac * MeltFac);
                for (int iter = 0; iter < 10; iter++) {
                    double A_over_c = A / coverage;
                    double new_coverage = tanh(A_over_c);
                    double sech2_val = 1.0 - new_coverage * new_coverage;  // sech²(x) = 1 - tanh²(x)
                    
                    double f = coverage - new_coverage;
                    double df = 1.0 + sech2_val * A / (coverage * coverage);
                    
                    double delta = f / df;
                    coverage -= delta;
                    if (fabs(delta) < 1.0e-4) {
                        break;
                    }
                }
                snow[veg].snow_depth = snow[veg].swq / (coverage * snow_density);
                snow[veg].coverage = coverage;
            }
            // set snow layer properties
            distribute_snow_state(step_dt, wind, 
                                  air_temp,
                                  pressure, &snow[veg]);
        }
    }

    /*********************************
      Initialize Surface Water
    *********************************/
    for (veg = 0; veg <= Nveg; veg++) {
        if (cell[veg].IS_WET) {
            cell[veg].h2osfc = 0.0;
            cell[veg].frac_h2o = 0.0;
        }
        else if (cell[veg].IS_GLAC) {
            cell[veg].frac_h2o = veg_con[veg].Cv;
            double glac_volume = 0.0365 * pow(soil_con->cell_area * cell[veg].frac_h2o, 1.375);
            cell[veg].h2osfc = glac_volume / soil_con->cell_area * 
                                MM_PER_M * (CONST_RHOICE / CONST_RHOFW);
        }
    }

    /*********************************
      Initialize Nnode_parameters
    *********************************/
    size_t Nbedrock = soil_con->Nbedrock;
    size_t Nsoil = Nbedrock - 1;
    for (veg = 0; veg <= Nveg; veg++) {
        Nsnow = snow[veg].Nsnow;
        if (veg_con[veg].Cv > 0) {
            cell[veg].Nsoil = Nsoil;
            cell[veg].Nnode = Nbedrock + Nsnow;
            // If surface water exists, add an additional node for the surface water layer
            if (cell[veg].h2osfc > 0.0) {
                cell[veg].Nnode += 1;
            }
        }
    }

    /*********************************
       Initialize water table depth
    *********************************/
    for (veg = 0; veg <= Nveg; veg++) {
        if (veg_con[veg].Cv > 0) {
            cell[veg].zwt = soil_con->init_zwt;
        }
    }

    /*********************************
       Initialize node temperatures
    *********************************/
    for (veg = 0; veg <= Nveg; veg++) {
        if (veg_con[veg].Cv > 0) {
            band = veg_con[veg].BandIndex;
            Nsnow = snow[veg].Nsnow;
            double T_soil = 0.0;
            double surf_temp = air_temp + soil_con->Tfactor[band];
            // Initialize snow node temperatures
            for (i = 0; i < Nsnow; i++) {
                energy[veg].T[i] = snow[veg].pack_T[i];
                energy[veg].last_T[i] = snow[veg].pack_T[i];
            }
            /* Initialize soil node temperatures */
            for (k = Nsnow; k < cell[veg].Nnode; k++) {
                lidx = k - Nsnow;
                if (zc_soil[lidx] <= 5.0) {
                    T_soil = surf_temp + (soil_con->avg_temp - 
                                surf_temp) * (zc_soil[lidx] / 5.0);
                }
                else if (zc_soil[lidx] <= 10.0) {
                    T_soil = soil_con->avg_temp;
                }
                else {
                    T_soil = soil_con->avg_temp + 0.05 * (zc_soil[lidx] - 10.0);
                }
                energy[veg].T[k] = T_soil;
                energy[veg].last_T[k] = T_soil;
                cell[veg].soil_T[lidx] = T_soil;
            }
            /* Initial estimate of temperatures */
            energy[veg].Tfoliage = surf_temp;
            energy[veg].Tstem = surf_temp;
            energy[veg].Tcanopy = surf_temp;
            energy[veg].Tgrnd = surf_temp;
            energy[veg].Tsurf = surf_temp;
        }
    }

    /******************************
       Initialize soil moistures 
    ******************************/
    for (veg = 0; veg <= Nveg; veg++) {
        if (veg_con[veg].Cv > 0) {            
            /* Initialize soil moistures */
            for (lidx = 0; lidx < Nsoil; lidx++) {
                // 温度大于0，地下水位以上设为田间持水量，地下水位以下设为饱和含水量
                if (soil_con->Zsum_soil[lidx] <= cell[veg].zwt) {
                    if (cell[veg].soil_T[lidx] >= CONST_TKFRZ) {
                        cell[veg].moist[lidx] = soil_con->Wsat_node[lidx] * 0.7;
                    }
                    else {
                        cell[veg].moist[lidx] = soil_con->Wsat_node[lidx] * 0.9;
                    }
                }
                else {
                    cell[veg].moist[lidx] = soil_con->Wsat_node[lidx];
                }
                cell[veg].liq[lidx] = cell[veg].moist[lidx]; // 将liq初始化为moist

                if (cell[veg].moist[lidx] >
                    soil_con->Wsat_node[lidx]) {
                    cell[veg].moist[lidx] =
                            soil_con->Wsat_node[lidx];
                }
            }
        }
    }

    /******************************************
       Compute soil thermal node properties
    ******************************************/
    for (veg = 0; veg <= Nveg; veg++) {
        // Initialize soil for existing vegetation types
        if (veg_con[veg].Cv > 0) {
            // set soil moisture properties for all soil thermal nodes
            ErrorFlag =
                    distribute_node_moisture_properties(
                        &cell[veg],
                        soil_con);
            if (ErrorFlag == ERROR) {
                log_err("Error setting physical properties for "
                        "soil thermal nodes"); 
            }
            // initialize last step soil layer properties
            for (i = 0; i < cell[veg].Nsoil; i++) {
                cell[veg].last_ice[i] = cell[veg].ice[i];
                cell[veg].last_liq[i] = cell[veg].liq[i];
                cell[veg].last_matric[i] = cell[veg].matric[i];
            }
            /* Soil and ice thermal properties for the layer */
            prepare_full_energy(pressure, &cell[veg], 
                                &energy[veg],
                                &snow[veg], soil_con);
            for (i = 0; i <= snow[veg].Nsnow; i++) {
                snow[veg].last_enthalpy[i] = snow[veg].enthalpy[i];
            }
        }
    }

    /******************************************
       initialize vegetation potential and actual water stress
    ******************************************/
    double root_psi = 0.0;
    double total_root = 0.0;
    for (veg = 0; veg <= Nveg; veg++) {
        if (veg_con[veg].Cv > 0.0 && !cell[veg].IS_GLAC) { // 冰川不计算植被水势
            veg_class = veg_con[veg].veg_class;
            Canopy_Upper = veg_lib[veg_class].Canopy_Upper;
            
            for (i = 0; i < veg_var[veg].Nroot; i++) {
                root_psi += cell[veg].root[i] * (cell[veg].matric[i] - soil_con->zc_soil[i]);
                total_root += cell[veg].root[i];
            }
            if (total_root > 0.0) {
                root_psi /= total_root;  // 加权平均根区水势
            } 
            else {
                root_psi = cell[veg].matric[0] - soil_con->zc_soil[0];
            }
            veg_var[veg].mat_VEG[3] = root_psi; // 根区水势(root)
            veg_var[veg].mat_VEG[2] = root_psi - Canopy_Upper; // 木质部水势(xyl)
            veg_var[veg].mat_VEG[1] = veg_var[veg].mat_VEG[2];  /* sunlit */
            veg_var[veg].mat_VEG[0] = veg_var[veg].mat_VEG[2];  /* shaded */
        }
    }
}