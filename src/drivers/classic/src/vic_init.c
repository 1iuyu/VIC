/******************************************************************************
 * @section DESCRIPTION
 *
 * Initialize model parameters
 *****************************************************************************/

#include "vic_driver_classic.h"

/******************************************************************************
 * @brief    Initialize grid cell model data structures.
 *****************************************************************************/
void
vic_init(size_t            Nveg_type,
         grid_cell_struct *grid,
         filep_struct     *filep,
         dmy_struct       *dmy)
{
    extern global_param_struct global_param;

    /** Read Grid Cell Vegetation Parameters **/
    grid->veg_con = read_vegparam(filep->vegparam, grid->gridcel,
                            Nveg_type);
    // 计算根区分数
    calc_root_fractions(grid);
    
    /** Read Elevation Band Data if Used **/
    read_snowband(filep->snowband, grid->soil_con);

    /** Make Top-level Control Structure **/
    grid->all_vars = make_all_vars(grid->veg_con[0].vegetat_type_num);

    /** allocate memory for the veg_hist_struct **/
    alloc_veg_hist(global_param.nrecs,
                   grid->veg_con[0].vegetat_type_num,
                   &grid->veg_hist);

    /**************************************************
        Initialize Meteological Forcing Values That
        Have not Been Specifically Set
    **************************************************/
    alloc_force(grid->force);

    vic_force(grid, dmy, filep->forcing);

    // populate model state, either using a cold start or from a restart file
    vic_populate_model_state(grid, filep, grid->gridcel);

}