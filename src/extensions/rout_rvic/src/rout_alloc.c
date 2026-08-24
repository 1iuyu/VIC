/******************************************************************************
 * @section DESCRIPTION
 *
 * Allocate memory for Routing structures.
 *****************************************************************************/

#include "rout.h"

/******************************************************************************
 * @brief    Allocate memory for Routing structures.
 *****************************************************************************/
void
rout_alloc(void)
{
    extern domain_struct local_domain;
    extern rout_struct  *rout;

    // Allocate memory in rout param_struct
    rout = malloc(local_domain.ncells_active * sizeof(*rout));
    check_alloc_status(rout, "Memory allocation error.");

    // 初始化
    memset(rout, 0, local_domain.ncells_active * sizeof(*rout));
}