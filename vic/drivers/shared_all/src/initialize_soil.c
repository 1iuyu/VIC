/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes the soil variable arrays for each new grid cell.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

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
        cell[veg].Ra_evap = 0.0;
        cell[veg].Ra_leaf = 0.0;

        cell[veg].asat = 0.0;
        cell[veg].zwt = 2.5;  // initial water table depth (m)
        cell[veg].zwt_lumped = 0.0;
        cell[veg].VPcanopy = 0.0;
        cell[veg].CLitter = 0.0;
        cell[veg].CInter = 0.0;
        cell[veg].CSlow = 0.0;
        for (lidx = 0; lidx < MAX_NODES; lidx++) {
            cell[veg].ice[lidx] = 0.0;
            cell[veg].liq[lidx] = 0.0;
            cell[veg].moist[lidx] = 0.0;
            cell[veg].soil_T[lidx] = 0.0;
            cell[veg].SnowFrac[lidx] = 0.0; // ????
            cell[veg].soil_frost[lidx] = 0.0;
            cell[veg].diffusivity[lidx] = 0.0;
            cell[veg].conductivity[lidx] = 0.0;
            cell[veg].soil_transp[lidx] = 0.0;
            cell[veg].transp_fact[lidx] = 0.0;
        }
        for (lidx = 0; lidx < MAX_SNOWS + MAX_NODES; lidx++) {
            cell[veg].PhaseChange[i] = 0;
            cell[veg].dz_node[lidx] = 0.0;
            cell[veg].Zsum_node[lidx] = 0.0;
            cell[veg].zc_node[lidx] = 0.0;
        }
        cell[veg].rootmoist = 0.0;
        cell[veg].wetness = 0.0;
        // Fluxes
        cell[veg].baseflow = 0.0;
        cell[veg].runoff = 0.0;
        cell[veg].soil_inflow = 0.0;
        cell[veg].evap = 0.0;
        cell[veg].vapor_grnd = 0.0;
        cell[veg].conden_grnd = 0.0;
        cell[veg].snowfrost = 0.0;
        cell[veg].snow_sublim = 0.0;
        // Soil terms
        cell[veg].esoil = 0.0;
        cell[veg].dewsoil = 0.0;
        cell[veg].soil_excess = 0.0;
        cell[veg].total_transp = 0.0;
        // Canopy terms
        cell[veg].transp = 0.0;
        cell[veg].canopyevap = 0.0;
        cell[veg].canopydew = 0.0;
        cell[veg].canopyfrost = 0.0;
        cell[veg].canopy_sublim = 0.0;
        cell[veg].canopy_vapor = 0.0;
        // Carbon terms
        cell[veg].RhLitter = 0.0;
        cell[veg].RhLitter2Atm = 0.0;
        cell[veg].RhInter = 0.0;
        cell[veg].RhSlow = 0.0;
        cell[veg].RhTot = 0.0;
    }
}
