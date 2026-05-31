/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine makes an array of vegitation variables for each vegitation type.
 *****************************************************************************/

#include "vic_driver_shared_all.h"

/******************************************************************************
 * @brief    Make an array of vegitation variables for each vegitation type.
 *****************************************************************************/
veg_var_struct *
make_veg_var(size_t veg_type_num)
{
    veg_var_struct      *temp = NULL;

    temp = calloc(veg_type_num, sizeof(*temp));
    check_alloc_status(temp, "Memory allocation error.");

    return temp;
}
