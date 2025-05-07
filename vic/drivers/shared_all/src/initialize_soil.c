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
    extern option_struct options;

    size_t               veg, lindex, frost_area;

    for (veg = 0; veg <= veg_num; veg++) {
            // Prognostic states
            cell[veg].aero_resist[0] = 0.0;
            cell[veg].aero_resist[1] = 0.0;
            cell[veg].CLitter = 0.0;
            cell[veg].CInter = 0.0;
            cell[veg].CSlow = 0.0;
            for (lindex = 0; lindex < options.Nlayer; lindex++) {
                cell[veg].layer[lindex].Cs = 0.0;
                cell[veg].layer[lindex].T = 0.0;
                for (frost_area = 0; frost_area < options.Nfrost;
                     frost_area++) {
                    cell[veg].layer[lindex].ice[frost_area] = 0.0;
                }
                cell[veg].layer[lindex].kappa = 0.0;
                cell[veg].layer[lindex].moist = 0.0;
                cell[veg].layer[lindex].phi = 0.0;
            }
            cell[veg].rootmoist = 0.0;
            cell[veg].wetness = 0.0;
            // Fluxes
            cell[veg].pot_evap = 0.0;
            cell[veg].baseflow = 0.0;
            cell[veg].runoff = 0.0;
            cell[veg].inflow = 0.0;
            cell[veg].RhLitter = 0.0;
            cell[veg].RhLitter2Atm = 0.0;
            cell[veg].RhInter = 0.0;
            cell[veg].RhSlow = 0.0;
            cell[veg].RhTot = 0.0;
            for (lindex = 0; lindex < options.Nlayer; lindex++) {
                cell[veg].layer[lindex].esoil = 0.0;
                cell[veg].layer[lindex].transp = 0.0;
                cell[veg].layer[lindex].evap = 0.0;
            }
    }
}
