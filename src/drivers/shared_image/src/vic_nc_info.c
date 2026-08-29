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
                nc_file_struct    *nc_hist_file,
                nc_var_struct     *nc_var)
{
    size_t i;

    // set datatype
    nc_var->nc_type = get_nc_dtype(dtype);

    for (i = 0; i < MAXDIMS; i++) {
        nc_var->nc_dimids[i] = -1;
        nc_var->nc_counts[i] = 0;
    }

    // Set the number of dimensions and the count sizes
    switch (varid) {
    case OUT_RA_OVER:
    case OUT_RA_SUB:
    case OUT_RA_GRND:
        nc_var->nc_dims = 4;
        nc_var->nc_counts[1] = nc_hist_file->turbul_size;
        nc_var->nc_counts[2] = nc_hist_file->nj_size;
        nc_var->nc_counts[3] = nc_hist_file->ni_size;
        break;
    case OUT_VEG_MATRIC:
        nc_var->nc_dims = 4;
        nc_var->nc_counts[1] = nc_hist_file->vegmat_size;
        nc_var->nc_counts[2] = nc_hist_file->nj_size;
        nc_var->nc_counts[3] = nc_hist_file->ni_size;
        break;
    case OUT_MATRIC:
    case OUT_SOIL_ICE:
    case OUT_SOIL_LIQ:
    case OUT_SOIL_MOIST:
    case OUT_SOIL_TEMP:
        nc_var->nc_dims = 4;
        nc_var->nc_counts[1] = nc_hist_file->soil_size;
        nc_var->nc_counts[2] = nc_hist_file->nj_size;
        nc_var->nc_counts[3] = nc_hist_file->ni_size;
        break;
    case OUT_SNOW_DENSITY:
    case OUT_SNOW_PACK_ICE:
    case OUT_SNOW_PACK_LIQ:
    case OUT_SNOW_ICEFRAC:
    case OUT_SNOW_LIQFRAC:
    case OUT_SNOW_POROSITY:
    case OUT_SNOW_RADIUS:
    case OUT_PACK_OUTFLOW:
    case OUT_SNOW_MELT:
    case OUT_SNOW_FRZE:
    case OUT_SNOW_PACK_TEMP:
        nc_var->nc_dims = 4;
        nc_var->nc_counts[1] = nc_hist_file->snow_size;
        nc_var->nc_counts[2] = nc_hist_file->nj_size;
        nc_var->nc_counts[3] = nc_hist_file->ni_size;
        break;
    default:
        nc_var->nc_dims = 3;
        nc_var->nc_counts[1] = nc_hist_file->nj_size;
        nc_var->nc_counts[2] = nc_hist_file->ni_size;
    }
}

/******************************************************************************
 * @brief    Set netcdf dim ids.
 *****************************************************************************/
void
set_nc_var_dimids(unsigned int    varid,
                  nc_file_struct *nc_hist_file,
                  nc_var_struct  *nc_var)
{
    size_t i;

    for (i = 0; i < MAXDIMS; i++) {
        nc_var->nc_dimids[i] = -1;
    }

    // Set the non-default ones
    switch (varid) {
    case OUT_RA_OVER:
    case OUT_RA_SUB:
    case OUT_RA_GRND:
        nc_var->nc_dimids[0] = nc_hist_file->time_dimid;
        nc_var->nc_dimids[1] = nc_hist_file->turbul_dimid;
        nc_var->nc_dimids[2] = nc_hist_file->nj_dimid;
        nc_var->nc_dimids[3] = nc_hist_file->ni_dimid;
        break;
    case OUT_VEG_MATRIC:
        nc_var->nc_dimids[0] = nc_hist_file->time_dimid;
        nc_var->nc_dimids[1] = nc_hist_file->vegmat_dimid;
        nc_var->nc_dimids[2] = nc_hist_file->nj_dimid;
        nc_var->nc_dimids[3] = nc_hist_file->ni_dimid;
        break;
    case OUT_MATRIC:
    case OUT_SOIL_ICE:
    case OUT_SOIL_LIQ:
    case OUT_SOIL_MOIST:
    case OUT_SOIL_TEMP:
        nc_var->nc_dimids[0] = nc_hist_file->time_dimid;
        nc_var->nc_dimids[1] = nc_hist_file->soil_dimid;
        nc_var->nc_dimids[2] = nc_hist_file->nj_dimid;
        nc_var->nc_dimids[3] = nc_hist_file->ni_dimid;
        break;
    case OUT_SNOW_DENSITY:
    case OUT_SNOW_PACK_ICE:
    case OUT_SNOW_PACK_LIQ:
    case OUT_SNOW_ICEFRAC:
    case OUT_SNOW_LIQFRAC:
    case OUT_SNOW_POROSITY:
    case OUT_SNOW_RADIUS:
    case OUT_PACK_OUTFLOW:
    case OUT_SNOW_MELT:
    case OUT_SNOW_FRZE:
    case OUT_SNOW_PACK_TEMP:
        nc_var->nc_dimids[0] = nc_hist_file->time_dimid;
        nc_var->nc_dimids[1] = nc_hist_file->snow_dimid;
        nc_var->nc_dimids[2] = nc_hist_file->nj_dimid;
        nc_var->nc_dimids[3] = nc_hist_file->ni_dimid;
        break;
    default:
        nc_var->nc_dimids[0] = nc_hist_file->time_dimid;
        nc_var->nc_dimids[1] = nc_hist_file->nj_dimid;
        nc_var->nc_dimids[2] = nc_hist_file->ni_dimid;
    }
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
