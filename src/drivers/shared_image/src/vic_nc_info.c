/******************************************************************************
 * @section DESCRIPTION
 *
 * Setup netCDF output variables.
 *****************************************************************************/

#include "vic_driver_shared_image.h"

/******************************************************************************
 * @brief    Setup netCDF output variables.
 *****************************************************************************/
void
set_nc_var_info(unsigned int       varid,
                unsigned short int dtype,
                unsigned short int domain,
                nc_file_struct    *nc_hist_file,
                nc_var_struct     *nc_var)
{
    size_t i;
    size_t dim = 0;
    extern metadata_struct out_metadata[N_OUTVAR_TYPES];

    // set datatype
    nc_var->nc_type = get_nc_dtype(dtype);

    for (i = 0; i < MAXDIMS; i++) {
        nc_var->nc_dimids[i] = -1;
        nc_var->nc_counts[i] = 0;
    }
    /* Time dimension is always the first dimension */
    nc_var->nc_counts[dim++] = 1;

    /* HRU dimension, if requested */
    if (domain == OUT_DOMAIN_HRU) {
        nc_var->nc_counts[dim++] = nc_hist_file->hru_size;
    }

    // Set the number of dimensions and the count sizes
    switch (out_metadata[varid].elem_type) {
    case OUT_ELEM_TURBUL:
        nc_var->nc_counts[dim++] = nc_hist_file->turbul_size;
        break;

    case OUT_ELEM_VEGMAT:
        nc_var->nc_counts[dim++] = nc_hist_file->vegmat_size;
        break;

    case OUT_ELEM_SOIL:
        nc_var->nc_counts[dim++] = nc_hist_file->soil_size;
        break;

    case OUT_ELEM_SNOW:
        nc_var->nc_counts[dim++] = nc_hist_file->snow_size;
        break;
    
    case OUT_ELEM_NODE:
        nc_var->nc_counts[dim++] = nc_hist_file->node_size;
        break;

    case OUT_ELEM_DEFAULT:
        break;

    default:
        log_err("Unknown output element type for variable %u", varid);
    }
    /* Spatial dimensions */
    nc_var->nc_counts[dim++] = nc_hist_file->nj_size;
    nc_var->nc_counts[dim++] = nc_hist_file->ni_size;

    nc_var->nc_dims = dim;
}

/******************************************************************************
 * @brief    Set netcdf dim ids.
 *****************************************************************************/
void
set_nc_var_dimids(unsigned int    varid,
                  unsigned short int domain,
                  nc_file_struct *nc_hist_file,
                  nc_var_struct  *nc_var)
{
    size_t i;
    size_t dim = 0;

    extern metadata_struct out_metadata[N_OUTVAR_TYPES];

    for (i = 0; i < MAXDIMS; i++) {
        nc_var->nc_dimids[i] = -1;
    }

    /* Time */
    nc_var->nc_dimids[dim++] = nc_hist_file->time_dimid;

    /* HRU */
    if (domain == OUT_DOMAIN_HRU) {
        nc_var->nc_dimids[dim++] = nc_hist_file->hru_dimid;
    }

    /* Element */
    switch (out_metadata[varid].elem_type) {
    case OUT_ELEM_TURBUL:
        nc_var->nc_dimids[dim++] = nc_hist_file->turbul_dimid;
        break;

    case OUT_ELEM_VEGMAT:
        nc_var->nc_dimids[dim++] = nc_hist_file->vegmat_dimid;
        break;

    case OUT_ELEM_SOIL:
        nc_var->nc_dimids[dim++] = nc_hist_file->soil_dimid;
        break;

    case OUT_ELEM_SNOW:
        nc_var->nc_dimids[dim++] = nc_hist_file->snow_dimid;
        break;

    case OUT_ELEM_NODE:
        nc_var->nc_dimids[dim++] = nc_hist_file->node_dimid;
        break;

    case OUT_ELEM_DEFAULT:
        break;

    default:
        log_err("Unknown output element type for variable %u", varid);
    }

    /* Spatial dimensions */
    nc_var->nc_dimids[dim++] = nc_hist_file->nj_dimid;
    nc_var->nc_dimids[dim++] = nc_hist_file->ni_dimid;

}

/******************************************************************************
 * @brief    Determine the netCDF file format
 *****************************************************************************/
int
get_nc_mode(unsigned short int format)
{
    int nc_format;

    switch (format) {
    case NETCDF3_CLASSIC:
        nc_format = NC_CLASSIC_MODEL;
        break;
    case NETCDF3_64BIT_OFFSET:
        nc_format = NC_64BIT_OFFSET;
        break;
    case NETCDF4_CLASSIC:
        nc_format = NC_CLASSIC_MODEL;
        break;
    case NETCDF4:
        nc_format = NC_NETCDF4;
        break;
    default:
        log_err("Unrecognized netCDF file format");
    }
    return nc_format;
}

/******************************************************************************
 * @brief    Determine the netCDF data type
 *****************************************************************************/
int
get_nc_dtype(unsigned short int dtype)
{
    int type;

    switch (dtype) {
    case OUT_TYPE_DEFAULT:
        type = OUT_TYPE_DOUBLE;
        break;
    case OUT_TYPE_CHAR:
        type = NC_CHAR;
        break;
    case OUT_TYPE_SINT:
        type = NC_SHORT;
        break;
    case OUT_TYPE_USINT:
        type = NC_UINT;
        break;
    case OUT_TYPE_INT:
        type = NC_INT;
        break;
    case OUT_TYPE_FLOAT:
        type = NC_FLOAT;
        break;
    case OUT_TYPE_DOUBLE:
        type = NC_DOUBLE;
        break;
    default:
        log_err("Unrecognized netCDF variable datatype: %hu", dtype);
    }
    return type;
}
