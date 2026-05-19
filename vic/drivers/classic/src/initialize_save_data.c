
/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes the save data structure for each grid cell.
 * It calls the put_data function to populate the save data storage terms.
 * It also initializes the output data structure to zero.
 *****************************************************************************/

#include "vic_driver_classic.h"

/******************************************************************************
 * @brief    Initialize the save data structure.
 *****************************************************************************/
void
initialize_save_data(grid_cell_struct *grid,
                     veg_lib_struct   *veg_lib,
                     double          **out_data,
                     timer_struct     *timer)
{
    all_vars_struct   *all_vars;
    force_data_struct *force;
    soil_con_struct   *soil_con;
    veg_con_struct    *veg_con;
    save_data_struct  *save_data;

    all_vars = grid->all_vars;
    force = grid->force;
    soil_con = grid->soil_con;
    veg_con = grid->veg_con;
    save_data = grid->save_data;
    
    // Calling put data will populate the save data storage terms
    put_data(all_vars, force, soil_con, veg_con, veg_lib,
             out_data, save_data, timer);

    zero_output_list(out_data);
}