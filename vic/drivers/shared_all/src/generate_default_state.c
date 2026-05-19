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
generate_default_state(all_vars_struct *all_vars,
                       soil_con_struct *soil_con,
                       veg_con_struct  *veg_con)
{
    extern option_struct     options;

    size_t                   Nveg;
    size_t                   veg;
    size_t                   band;
    size_t                   Nsnow;
    size_t                   lidx;
    size_t                   i, k;
    double                   Cv;

    cell_data_struct        *cell;
    snow_data_struct        *snow;
    energy_bal_struct       *energy;
    cell = all_vars->cell;
    snow = all_vars->snow;
    energy = all_vars->energy;
    Nveg = veg_con[0].vegetat_type_num;
    double *gravel_node = soil_con->gravel_node;
    double *bulk_dens_node = soil_con->bulk_dens_node;

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
       (require user to use a state file if they want control over initial moist)
    ******************************************************************************/
    for (veg = 0; veg <= Nveg; veg++) {
        Cv = veg_con[veg].Cv;
        if (Cv > 0) {
            // 初始化土壤水力学参数
            if (options.DENSITY_FROM_SOIL) {
                for (lidx = 0; lidx < Nsoil; lidx++) {
                    bulk_dens_node[lidx] = (bulk_dens_node[lidx] * 
                        (1.0 - gravel_node[lidx]) + gravel_node[lidx] * 2650);
                }
            }
            else {
                PedoTransfer(soil_con); // 土壤参数从PedoTransfer函数中计算得到
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
        Cv = veg_con[veg].Cv;
        if (Cv > 0) {
            band = veg_con[veg].BandIndex;
            Nsnow = snow[veg].Nsnow;
            // Initialize snow node temperatures
            for (i = 0; i < Nsnow; i++) {
                energy[veg].T[i] = snow[veg].pack_T[i];
            }
            /* Initialize soil node temperatures */
            for (k = Nsnow; k < cell[veg].Nnode; k++) {
                lidx = k - Nsnow;
                if (options.FROZEN_SOIL) {
                    energy[veg].T[k] = 260.0;
                    cell[veg].soil_T[lidx] = 260.0;
                }
                else {
                    energy[veg].T[k] = CONST_TKFRZ;
                    cell[veg].soil_T[lidx] = CONST_TKFRZ;
                }
            }
            /* Initial estimate of temperatures */
            energy[veg].Tfoliage = 260.0 + soil_con->Tfactor[band];
            energy[veg].Tstem = 260.0 + soil_con->Tfactor[band];
            energy[veg].Tcanopy = 260.0 + soil_con->Tfactor[band];
            energy[veg].Tgrnd = 260.0 + soil_con->Tfactor[band];
            energy[veg].Tsurf = 260.0 + soil_con->Tfactor[band];
        }
    }  
}
