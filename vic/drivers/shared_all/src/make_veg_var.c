/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine makes an array of vegitation variables for each vegitation type.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    Make an array of vegitation variables for each vegitation type.
 *****************************************************************************/
veg_var_struct *
make_veg_var(size_t veg_type_num)
{
    extern option_struct options;

    size_t               i;
    veg_var_struct     *temp = NULL;

    temp = calloc(veg_type_num, sizeof(*temp));
    check_alloc_status(temp, "Memory allocation error.");

    if (options.CARBON) {
        for (i = 0; i < veg_type_num; i++) {
            temp[i].NscaleFactor = calloc(options.Ncanopy,
                                          sizeof(*(temp[i].
                                                   NscaleFactor)));
            check_alloc_status(temp[i].NscaleFactor,
                               "Memory allocation error.");
            temp[i].aPARLayer = calloc(options.Ncanopy,
                                       sizeof(*(temp[i].aPARLayer)));
            check_alloc_status(temp[i].aPARLayer,
                               "Memory allocation error.");
            temp[i].CiLayer = calloc(options.Ncanopy,
                                     sizeof(*(temp[i].CiLayer)));
            check_alloc_status(temp[i].CiLayer,
                               "Memory allocation error.");
            temp[i].rsLayer = calloc(options.Ncanopy,
                                     sizeof(*(temp[i].rsLayer)));
            check_alloc_status(temp[i].rsLayer,
                               "Memory allocation error.");
        }
    }

    return temp;
}
