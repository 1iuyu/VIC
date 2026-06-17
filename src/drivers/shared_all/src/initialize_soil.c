/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes the soil variable arrays for each new grid cell.
 *****************************************************************************/

#include "vic_driver_shared_all.h"

/******************************************************************************
 * @brief    This routine initializes the soil variable arrays for each new
 *           grid cell.
 *****************************************************************************/
void
initialize_soil(cell_data_struct *cell,
                size_t            veg_num)
{
    size_t               i, veg, lidx;

    for (veg = 0; veg <= veg_num; veg++) {
        // Prognostic states
        for (i = 0; i < 3; i++) {
            cell[veg].Ra_over[i] = 0.0;
            cell[veg].Ra_sub[i] = 0.0;
            cell[veg].Ra_grnd[i] = 0.0;
        }
        cell[veg].IS_GLAC = false;
        cell[veg].IS_VEG = false;
        cell[veg].Nsoil = 0;
        cell[veg].Nroot = 0;
        cell[veg].Nnode = 0;
        cell[veg].Ra_evap = 0.0;
        cell[veg].Ra_leaf = 0.0;
        cell[veg].Ra_stem = 0.0;
        cell[veg].Qair_grnd = 0.0;
        cell[veg].Qair_soil = 0.0;
        cell[veg].Qair_snow = 0.0;
        cell[veg].Qair_deriv = 0.0;
        for (i = 0; i < 3; i++) {
            cell[veg].Z0m_grnd[i] = 0.0;
            cell[veg].Z0m_sub[i] = 0.0;
            cell[veg].ref_height[i] = 0.0;
        }        
        for (i = 0; i < 2; i++) {
            cell[veg].displacement[i] = 0.0;
        } 
        cell[veg].asat = 0.0;
        cell[veg].zwt = 0.0;  // initial water table depth (m)
        for (lidx = 0; lidx < MAX_SOILS; lidx++) {
            cell[veg].ice[lidx] = 0.0;
            cell[veg].liq[lidx] = 0.0;
            cell[veg].last_ice[lidx] = 0.0;
            cell[veg].last_liq[lidx] = 0.0;
            cell[veg].lateral_flow[lidx] = 0.0;
            cell[veg].deriv_vapor[lidx] = 0.0;
            cell[veg].moist[lidx] = 0.0;
            cell[veg].soil_T[lidx] = 0.0;
            cell[veg].porosity[lidx] = 0.0;
            cell[veg].matric[lidx] = 0.0;
            cell[veg].last_matric[lidx] = 0.0;
            cell[veg].soil_imped[lidx] = 0.0;
            cell[veg].conductivity[lidx] = 0.0;
            cell[veg].conduct_int[lidx] = 0.0;
            cell[veg].transp_sink[lidx] = 0.0;
            cell[veg].root[lidx] = 0.0;
        }
        cell[veg].rootmoist = 0.0;
        // Fluxes
        cell[veg].baseflow = 0.0;
        cell[veg].runoff = 0.0;
        cell[veg].soil_inflow = 0.0;
        cell[veg].recharge = 0.0;
        cell[veg].storage_aqf = 5.0; // m
        cell[veg].evap = 0.0;
        cell[veg].snowfrost = 0.0;
        cell[veg].snow_sublim = 0.0;
        cell[veg].Nthaw = 0;
        cell[veg].Nfrost = 0;
        cell[veg].fdepth = 0.0;
        cell[veg].tdepth = 0.0;
        cell[veg].esoil = 0.0;
        cell[veg].dewsoil = 0.0;
        cell[veg].soil_excess = 0.0;
        cell[veg].transp_fact = 0.0;
        // Canopy terms
        cell[veg].transp = 0.0;
        cell[veg].canopyevap = 0.0;
        cell[veg].canopydew = 0.0;
        cell[veg].canopyfrost = 0.0;
        cell[veg].canopy_sublim = 0.0;
        cell[veg].canopy_vapor = 0.0;
    }
}

