/******************************************************************************
 * @section DESCRIPTION
 *
 * Allocate and free memory for the veg_hist data struct
 *****************************************************************************/

 #include "vic_driver_shared_image.h"

/******************************************************************************
 * @brief
 *****************************************************************************/
void
alloc_veg_hist(veg_hist_struct *veg_hist)
{
    veg_hist->fcanopy = calloc(NR + 1, sizeof(*(veg_hist->fcanopy)));
    check_alloc_status(veg_hist->fcanopy, "Memory allocation error.");

    veg_hist->LAI = calloc(NR + 1, sizeof(*(veg_hist->LAI)));
    check_alloc_status(veg_hist->LAI, "Memory allocation error.");

    veg_hist->SAI = calloc(NR + 1, sizeof(*(veg_hist->SAI)));
    check_alloc_status(veg_hist->SAI, "Memory allocation error.");
}

/******************************************************************************
 * @brief    Free veg hist structure.
 *****************************************************************************/
void
free_veg_hist(veg_hist_struct *veg_hist)
{
    if (veg_hist == NULL) {
        return;
    }
    free(veg_hist->fcanopy);
    free(veg_hist->LAI);
    free(veg_hist->SAI);
}
