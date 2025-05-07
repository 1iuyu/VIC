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

    size_t               i;
    energy_bal_struct  *temp = NULL;

    temp = calloc(nveg, sizeof(*temp));
    check_alloc_status(temp, "Memory allocation error.");

    /** Initialize all records to unfrozen conditions */
    for (i = 0; i < nveg; i++) {
        temp[i].frozen = false;
    }

    return temp;
}
