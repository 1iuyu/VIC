/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine makes an array of glacier cover data structures.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    Make an array of glacier cover data structures.
 *****************************************************************************/
glac_data_struct *
make_glac_data(size_t nveg)
{

    glac_data_struct   *temp = NULL;

    temp = calloc(nveg, sizeof(*temp));

    check_alloc_status(temp, "Memory allocation error.");

    return temp;
}
