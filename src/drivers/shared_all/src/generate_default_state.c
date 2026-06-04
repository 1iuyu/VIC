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
    extern option_struct     options;
    int         ErrorFlag;
    size_t      Nveg;
    size_t      veg;
    size_t      band;
    size_t      Nsnow;
    size_t      lidx;
    size_t      i, k;
    double      Cv;
    double      init_temp;
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
    double *gravel_node = soil_con->gravel_node;
    double *bulk_dens_node = soil_con->bulk_dens_node;
    double air_temp = force->air_temp[NR];
    double pressure = force->pressure[NR];
    /*********************************
      Initialize snowpack temperatures
    *********************************/
    for (veg = 0; veg <= Nveg; veg++) {
        Cv = veg_con[veg].Cv;
        if (Cv > 0) {
            snow[veg].swq = 0.0;   // [mm]
            snow[veg].last_swq = 0.0;   // [mm]
            // set snow layer properties
            distribute_snow_state(&snow[veg]);
        }
    }
    /*********************************
      Initialize Nnode_parameters
    *********************************/
    size_t Nbedrock = soil_con->Nbedrock;
    size_t Nsoil = Nbedrock - 1;
    for (veg = 0; veg <= Nveg; veg++) {
        Cv = veg_con[veg].Cv;
        Nsnow = snow[veg].Nsnow;
        if (Cv > 0) {
            // 温度节点[0 到 Nsnow-1]
            for (i = 0; i < Nsnow; i++) {
                cell[veg].dz_node[i] = snow[veg].dz_snow[i];
            }
            // 土壤节点+基岩节点
            for (i = 0; i < Nbedrock; i++) {
                lidx = Nsnow + i;
                cell[veg].dz_node[lidx] = soil_con->dz_soil[i];
            }
            cell[veg].Nsoil = Nsoil;
            cell[veg].Nnode = Nsoil + Nsnow + 1;
            /* Compute node depths */
            if (Nsnow > 0) {
                for (k = 0; k < Nsnow; k++) {
                    cell[veg].Zsum_node[k] = snow[veg].Zsum_snow[k];
                }
            }
            cell[veg].Zsum_node[Nsnow] = 0.0;
            double sum_dz = 0.;
            for (k = Nsnow; k < cell[veg].Nnode; k++) {
                sum_dz += cell[veg].dz_node[k];
                cell[veg].Zsum_node[k+1] = sum_dz;
            }
            // 节点中心坐标
            if (Nsnow > 0) {
                for (k = 0; k < Nsnow; k++) {
                    cell[veg].zc_node[k] = snow[veg].zc_snow[k];
                }
            }
            for (k = Nsnow; k < cell[veg].Nnode; k++) {
                cell[veg].zc_node[k] = cell[veg].Zsum_node[k+1] - 
                                        cell[veg].dz_node[k] / 2.0;
            }
        }
    }

    /******************************************************************************
       Initialize soil moistures
       currently setting moist to porosity as default and eliminate init_moist 
    ******************************************************************************/
    for (veg = 0; veg <= Nveg; veg++) {
        if (veg_con[veg].Cv > 0) {
            // 初始化土壤水力学参数
            if (options.PARAM_FROM_SOIL) {
                for (lidx = 0; lidx < Nsoil; lidx++) {
                    bulk_dens_node[lidx] = (bulk_dens_node[lidx] * 
                        (1.0 - gravel_node[lidx]) + gravel_node[lidx] * 2650);
                }
            }
            else {
                // 土壤参数从PedoTransfer函数中计算得到
                PedoTransfer(soil_con); 
            }
            
            /* Initialize soil moistures */
            for (lidx = 0; lidx < Nsoil; lidx++) {
                cell[veg].moist[lidx] =
                        soil_con->Wsat_node[lidx] * 0.6;
                cell[veg].liq[lidx] = cell[veg].moist[lidx]; // 将liq初始化为moist
                if (cell[veg].moist[lidx] >
                    soil_con->Wsat_node[lidx]) {
                    cell[veg].moist[lidx] =
                            soil_con->Wsat_node[lidx];
                }
            }
        }
    }

    /*********************************
       Initialize node temperatures
    *********************************/
    for (veg = 0; veg <= Nveg; veg++) {
        if (veg_con[veg].Cv > 0) {
            band = veg_con[veg].BandIndex;
            Nsnow = snow[veg].Nsnow;
            init_temp = air_temp + soil_con->Tfactor[band];
            // Initialize snow node temperatures
            for (i = 0; i < Nsnow; i++) {
                energy[veg].T[i] = snow[veg].pack_T[i];
            }
            /* Initialize soil node temperatures */
            for (k = Nsnow; k < cell[veg].Nnode; k++) {
                lidx = k - Nsnow;
                energy[veg].T[k] = air_temp;
                cell[veg].soil_T[lidx] = air_temp;
            }
            /* Initial estimate of temperatures */
            energy[veg].Tfoliage = init_temp;
            energy[veg].Tstem = init_temp;
            energy[veg].Tcanopy = init_temp;
            energy[veg].Tgrnd = init_temp;
            energy[veg].Tsurf = init_temp;
        }
    }
    /******************************************
       Compute soil thermal node properties
    ******************************************/
    for (veg = 0; veg <= Nveg; veg++) {
        // Initialize soil for existing vegetation types
        cell[veg].IS_GLAC = veg_con[veg].IS_GLAC;
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
            /* Soil and ice thermal properties for the layer */
            prepare_full_energy(pressure, &cell[veg], 
                                &energy[veg],
                                &snow[veg], soil_con);
        }
    }
    /******************************************
       Compute maximum daylight duration
    ******************************************/
    double max_daylen = calc_max_daylength(soil_con->lat);
    for (veg = 0; veg <= Nveg; veg++) {
        Cv = veg_con[veg].Cv;
        if (Cv > 0.0) {
            cell[veg].max_daylen = max_daylen;
        }
    }
    /******************************************
       initialize vegetation potential and actual water stress
    ******************************************/
    double root_psi = 0.0;
    double total_root = 0.0;
    for (veg = 0; veg <= Nveg; veg++) {
        if (veg_con[veg].Cv > 0.0 && !veg_con[veg].IS_GLAC) { // 冰川不计算植被水势
            veg_class = veg_con[veg].veg_class;
            Canopy_Upper = veg_lib[veg_class].Canopy_Upper;
            
            for (i = 0; i < veg_con[veg].Nroot; i++) {
                root_psi += veg_con[veg].root[i] * (cell[veg].matric[i] - soil_con->zc_soil[i]);
                total_root += veg_con[veg].root[i];
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
