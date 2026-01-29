/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine makes an array of frozen soil data structures, one for each
 * vegetation type and bare soil.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    This routine makes an array of frozen soil data structures, one
 *           for each vegetation type and bare soil.
 *****************************************************************************/
energy_bal_struct *
make_energy_bal(size_t nveg)
{

    energy_bal_struct  *temp = NULL;

    temp = calloc(nveg, sizeof(*temp));
    check_alloc_status(temp, "Memory allocation error.");

    return temp;
}
