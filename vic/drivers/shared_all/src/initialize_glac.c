 /******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initailizes the glacier variable array.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    This routine initailizes the glacier variable array.
 *****************************************************************************/ 
void
initialize_glac(glac_data_struct *glacier,
               size_t             veg_num)
{

    size_t               i;

    for (i = 0; i <= veg_num; i++) {
        // Prognostic states
        glacier[i].cold_content = 0.0;
        glacier[i].surf_temp = 0.0;
        glacier[i].surf_temp_fbcount = 0.0;
        glacier[i].surf_temp_fbflag = false;
        glacier[i].Qnet = 0.0;
        glacier[i].mass_balance = 0.0;
        glacier[i].ice_mass_balance = 0.0;
        glacier[i].vapor_flux = 0.0;
        glacier[i].water_storage = 0.0;
        glacier[i].outflow = 0.0;
        glacier[i].outflow_coef = 0.0;
        glacier[i].inflow = 0.0;
    }
}
