/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine frees all components of the veg_con structure.
 *****************************************************************************/

#include "vic_driver_shared_all.h"

/******************************************************************************
 * @brief    This subroutine frees all components of the veg_con structure.
 *****************************************************************************/
void
free_vegcon(veg_con_struct **veg_con)
{
    
    free((char *) veg_con[0]);
}
