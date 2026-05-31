/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine handles the startup tasks for the image driver.
 *****************************************************************************/

#include "vic_driver_image.h"

/******************************************************************************
 * @brief    Wrapper function for VIC startup tasks.
 *****************************************************************************/
void
vic_image_start(void)
{
    extern filep_struct     filep;
    extern filenames_struct filenames;
    extern int              mpi_rank;
    extern veg_lib_struct  *veg_lib;
    size_t       Nveg_type;

    // Initialize structures
    initialize_global_structures();

    if (mpi_rank == VIC_MPI_ROOT) {
        // Read the global parameter file
        filep.globalparam = open_file(filenames.global, "r");
        get_global_param(filep.globalparam);
        
        /** Read Vegetation Library File **/
        filep.veglib = open_file(filenames.veglib, "r");
        veg_lib = read_veglib(filep.veglib, &Nveg_type);
    }

    // initialize image mode structures and settings
    vic_start();

}
