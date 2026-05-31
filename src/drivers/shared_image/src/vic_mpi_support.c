/******************************************************************************
 * @section DESCRIPTION
 *
 * MPI support routines for VIC
 *****************************************************************************/

#include "vic_driver_shared_image.h"

/******************************************************************************
* @brief   Print MPI Error String to LOG_DEST, this function is used by loggers
******************************************************************************/
void
print_mpi_error_str(int error_code)
{
    extern int mpi_rank;

    char       error_string[MAXSTRING];
    int        length_of_error_string;

    if (error_code != MPI_SUCCESS) {
        MPI_Error_string(error_code, error_string, &length_of_error_string);
        fprintf(LOG_DEST, "MPI Error in rank %d\n%s\n", mpi_rank, error_string);
    }
}

/******************************************************************************
 * @brief   Initialize MPI functionality
 *****************************************************************************/
void
initialize_mpi(void)
{
    extern MPI_Datatype mpi_global_struct_type;
    extern MPI_Datatype mpi_filenames_struct_type;
    extern MPI_Datatype mpi_location_struct_type;
    extern MPI_Datatype mpi_option_struct_type;
    extern MPI_Datatype mpi_param_struct_type;
    extern MPI_Datatype mpi_alarm_struct_type;
    extern MPI_Comm     MPI_COMM_VIC;
    extern int          mpi_rank;
    extern int          mpi_size;
    int                 status;

    // get MPI mpi_rank and mpi_size
    status = MPI_Comm_rank(MPI_COMM_VIC, &mpi_rank);
    check_mpi_status(status, "MPI Error");

    // set mpi error handling
    MPI_Comm_set_errhandler(MPI_COMM_VIC, MPI_ERRORS_RETURN);

    status = MPI_Comm_size(MPI_COMM_VIC, &mpi_size);
    check_mpi_status(status, "MPI Error");

    // initialize MPI data structures
    create_MPI_global_struct_type(&mpi_global_struct_type);
    create_MPI_filenames_struct_type(&mpi_filenames_struct_type);
    create_MPI_location_struct_type(&mpi_location_struct_type);
    create_MPI_alarm_struct_type(&mpi_alarm_struct_type);
    create_MPI_option_struct_type(&mpi_option_struct_type);
    create_MPI_param_struct_type(&mpi_param_struct_type);
}

/******************************************************************************
 * @brief   Create an MPI_Datatype that represents the global_param_struct
 * @details This allows MPI operations in which the entire global_param_struct
 *          can be treated as an MPI_Datatype. NOTE: This function needs to be
 *          kept in-sync with the global_param_struct data type in vic_def.h.
 *
 * @param mpi_type MPI_Datatype that can be used in MPI operations
 *****************************************************************************/
void
create_MPI_global_struct_type(MPI_Datatype *mpi_type)
{
    extern MPI_Comm MPI_COMM_VIC;
    int             nitems; // number of elements in struct
    int             status;
    int            *blocklengths;
    size_t          i;
    MPI_Aint       *offsets;
    MPI_Datatype   *mpi_types;

    // nitems has to equal the number of elements in global_param_struct
    nitems = 26;
    blocklengths = malloc(nitems * sizeof(*blocklengths));
    check_alloc_status(blocklengths, "Memory allocation error.");

    offsets = malloc(nitems * sizeof(*offsets));
    check_alloc_status(offsets, "Memory allocation error.");

    mpi_types = malloc(nitems * sizeof(*mpi_types));
    check_alloc_status(mpi_types, "Memory allocation error.")

    // most of the elements in global_param_struct are not arrays. Use 1 as
    // the default block length and reset as needed
    for (i = 0; i < (size_t) nitems; i++) {
        blocklengths[i] = 1;
    }

    // reset i
    i = 0;

    // double resolution;
    offsets[i] = offsetof(global_param_struct, resolution);
    mpi_types[i++] = MPI_DOUBLE;

    // double dt;
    offsets[i] = offsetof(global_param_struct, step_dt);
    mpi_types[i++] = MPI_DOUBLE;

    // size_t model_steps_per_day;
    offsets[i] = offsetof(global_param_struct, model_steps_per_day);
    mpi_types[i++] = MPI_AINT;

    // unsigned short endday;
    offsets[i] = offsetof(global_param_struct, endday);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short endmonth;
    offsets[i] = offsetof(global_param_struct, endmonth);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short endyear;
    offsets[i] = offsetof(global_param_struct, endyear);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short forceday[2];
    offsets[i] = offsetof(global_param_struct, forceday);
    blocklengths[i] = 2;
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned int forcesec[2];
    offsets[i] = offsetof(global_param_struct, forcesec);
    blocklengths[i] = 2;
    mpi_types[i++] = MPI_UNSIGNED;

    // unsigned short forcemonth[2];
    offsets[i] = offsetof(global_param_struct, forcemonth);
    blocklengths[i] = 2;
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short forceoffset[2];
    offsets[i] = offsetof(global_param_struct, forceoffset);
    blocklengths[i] = 2;
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned int forceskip[2];
    offsets[i] = offsetof(global_param_struct, forceskip);
    blocklengths[i] = 2;
    mpi_types[i++] = MPI_UNSIGNED;

    // unsigned short int forceyear[2];
    offsets[i] = offsetof(global_param_struct, forceyear);
    blocklengths[i] = 2;
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // size_t nrecs;
    offsets[i] = offsetof(global_param_struct, nrecs);
    mpi_types[i++] = MPI_AINT;

    // unsigned short int startday;
    offsets[i] = offsetof(global_param_struct, startday);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned int startsec;
    offsets[i] = offsetof(global_param_struct, startsec);
    mpi_types[i++] = MPI_UNSIGNED;

    // unsigned short int startmonth;
    offsets[i] = offsetof(global_param_struct, startmonth);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short int startyear;
    offsets[i] = offsetof(global_param_struct, startyear);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short int stateday;
    offsets[i] = offsetof(global_param_struct, stateday);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short int statemonth;
    offsets[i] = offsetof(global_param_struct, statemonth);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned int statesec;
    offsets[i] = offsetof(global_param_struct, statesec);
    mpi_types[i++] = MPI_UNSIGNED;

    // unsigned short int stateyear;
    offsets[i] = offsetof(global_param_struct, stateyear);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short int calendar;
    offsets[i] = offsetof(global_param_struct, calendar);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short int time_units;
    offsets[i] = offsetof(global_param_struct, time_units);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // double time_origin_num;
    offsets[i] = offsetof(global_param_struct, time_origin_num);
    mpi_types[i++] = MPI_DOUBLE;

    // char time_origin_str[MAXSTRING];
    offsets[i] = offsetof(global_param_struct, time_origin_str);
    blocklengths[i] = MAXSTRING;
    mpi_types[i++] = MPI_CHAR;

    // make sure that the we have the right number of elements
    if (i != (size_t) nitems) {
        log_err("Miscount: %zd not equal to %d.", i, nitems);
    }

    status = MPI_Type_create_struct(nitems, blocklengths, offsets, mpi_types,
                                    mpi_type);
    check_mpi_status(status, "MPI error.");

    status = MPI_Type_commit(mpi_type);
    check_mpi_status(status, "MPI error.");

    // cleanup
    free(blocklengths);
    free(offsets);
    free(mpi_types);
}

/******************************************************************************
 * @brief   Create an MPI_Datatype that represents the filenames_struct
 * @details This allows MPI operations in which the entire filenames_struct
 *          can be treated as an MPI_Datatype. NOTE: This function needs to be
 *          kept in-sync with the global_param_struct data type in vic_def.h.
 *
 * @param mpi_type MPI_Datatype that can be used in MPI operations
 *****************************************************************************/
void
create_MPI_filenames_struct_type(MPI_Datatype *mpi_type)
{
    extern MPI_Comm MPI_COMM_VIC;

    int             nitems; // number of elements in struct
    int             status;
    int            *blocklengths;
    size_t          i;
    MPI_Aint       *offsets;
    MPI_Datatype   *mpi_types;

    // nitems has to equal the number of elements in filenames_struct
    nitems = 10;
    blocklengths = malloc(nitems * sizeof(*blocklengths));
    check_alloc_status(blocklengths, "Memory allocation error.");

    offsets = malloc(nitems * sizeof(*offsets));
    check_alloc_status(offsets, "Memory allocation error.");

    mpi_types = malloc(nitems * sizeof(*mpi_types));
    check_alloc_status(mpi_types, "Memory allocation error.");

    // most of the elements in filenames_struct are character arrays. Use
    // MAXSTRING as the default block length and reset as needed
    for (i = 0; i < (size_t) nitems; i++) {
        blocklengths[i] = MAXSTRING;
    }

    // reset i
    i = 0;

    // char forcing[2][MAXSTRING];
    offsets[i] = offsetof(filenames_struct, forcing);
    blocklengths[i] *= 2;
    mpi_types[i++] = MPI_CHAR;

    // char f_path_pfx[2][MAXSTRING];
    offsets[i] = offsetof(filenames_struct, f_path_pfx);
    blocklengths[i] *= 2;
    mpi_types[i++] = MPI_CHAR;

    // char global[MAXSTRING];
    offsets[i] = offsetof(filenames_struct, global);
    mpi_types[i++] = MPI_CHAR;

    // char domain[MAXSTRING];
    offsets[i] = offsetof(filenames_struct, domain);
    mpi_types[i++] = MPI_CHAR;

    // char constants[MAXSTRING];
    offsets[i] = offsetof(filenames_struct, constants);
    mpi_types[i++] = MPI_CHAR;

    // char init_state[MAXSTRING];
    offsets[i] = offsetof(filenames_struct, init_state);
    mpi_types[i++] = MPI_CHAR;

    // char params[MAXSTRING];
    offsets[i] = offsetof(filenames_struct, params);
    mpi_types[i++] = MPI_CHAR;

    // char result_dir[MAXSTRING];
    offsets[i] = offsetof(filenames_struct, result_dir);
    mpi_types[i++] = MPI_CHAR;

    // char statefile[MAXSTRING];
    offsets[i] = offsetof(filenames_struct, statefile);
    mpi_types[i++] = MPI_CHAR;

    // char log_path[MAXSTRING];
    offsets[i] = offsetof(filenames_struct, log_path);
    mpi_types[i++] = MPI_CHAR;


    // make sure that the we have the right number of elements
    if (i != (size_t) nitems) {
        log_err("Miscount: %zd not equal to %d.", i, nitems);
    }

    status = MPI_Type_create_struct(nitems, blocklengths, offsets, mpi_types,
                                    mpi_type);
    check_mpi_status(status, "MPI error.");

    status = MPI_Type_commit(mpi_type);
    check_mpi_status(status, "MPI error.");

    // cleanup
    free(blocklengths);
    free(offsets);
    free(mpi_types);
}

/******************************************************************************
 * @brief   Create an MPI_Datatype that represents the location_struct
 * @details This allows MPI operations in which the entire location_struct
 *          can be treated as an MPI_Datatype. NOTE: This function needs to be
 *          kept in-sync with the location_struct data type in
 *          vic_image_driver.h.
 *
 * @param mpi_type MPI_Datatype that can be used in MPI operations
 *****************************************************************************/
void
create_MPI_location_struct_type(MPI_Datatype *mpi_type)
{
    extern MPI_Comm MPI_COMM_VIC;

    int             nitems; // number of elements in struct
    int             status;
    int            *blocklengths;
    size_t          i;
    MPI_Aint       *offsets;
    MPI_Datatype   *mpi_types;

    // nitems has to equal the number of elements in location_struct
    nitems = 9;
    blocklengths = malloc(nitems * sizeof(*blocklengths));
    check_alloc_status(blocklengths, "Memory allocation error.");

    offsets = malloc(nitems * sizeof(*offsets));
    check_alloc_status(offsets, "Memory allocation error.");

    mpi_types = malloc(nitems * sizeof(*mpi_types));
    check_alloc_status(mpi_types, "Memory allocation error.");

    // none of the elements in location_struct are arrays.
    for (i = 0; i < (size_t) nitems; i++) {
        blocklengths[i] = 1;
    }

    // reset i
    i = 0;

    // size_t run;
    offsets[i] = offsetof(location_struct, run);
    mpi_types[i++] = MPI_C_BOOL;

    // double latitude;
    offsets[i] = offsetof(location_struct, latitude);
    mpi_types[i++] = MPI_DOUBLE;

    // double longitude;
    offsets[i] = offsetof(location_struct, longitude);
    mpi_types[i++] = MPI_DOUBLE;

    // double area;
    offsets[i] = offsetof(location_struct, area);
    mpi_types[i++] = MPI_DOUBLE;

    // double frac;
    offsets[i] = offsetof(location_struct, frac);
    mpi_types[i++] = MPI_DOUBLE;

    // size_t nveg;
    offsets[i] = offsetof(location_struct, nveg);
    mpi_types[i++] = MPI_AINT;

    // size_t global_idx;
    offsets[i] = offsetof(location_struct, global_idx);
    mpi_types[i++] = MPI_AINT;

    // size_t io_idx;
    offsets[i] = offsetof(location_struct, io_idx);
    mpi_types[i++] = MPI_AINT;

    // size_t local_idx;
    offsets[i] = offsetof(location_struct, local_idx);
    mpi_types[i++] = MPI_AINT;

    // make sure that the we have the right number of elements
    if (i != (size_t) nitems) {
        log_err("Miscount: %zd not equal to %d.", i, nitems);
    }

    status = MPI_Type_create_struct(nitems, blocklengths, offsets, mpi_types,
                                    mpi_type);
    check_mpi_status(status, "MPI error.");

    status = MPI_Type_commit(mpi_type);
    check_mpi_status(status, "MPI error.");

    // cleanup
    free(blocklengths);
    free(offsets);
    free(mpi_types);
}

/******************************************************************************
 * @brief   Create an MPI_Datatype that represents the option_struct
 * @details This allows MPI operations in which the entire option_struct can
 *          be treated as an MPI_Datatype. NOTE: This function needs to be kept
 *          in-sync with the option_struct data type in vic_def.h.
 *
 * @param mpi_type MPI_Datatype that can be used in MPI operations
 *****************************************************************************/
void
create_MPI_option_struct_type(MPI_Datatype *mpi_type)
{
    extern MPI_Comm MPI_COMM_VIC;

    int             nitems; // number of elements in struct
    int             status;
    int            *blocklengths;
    size_t          i;
    MPI_Aint       *offsets;
    MPI_Datatype   *mpi_types;

    // nitems has to equal the number of elements in option_struct
    nitems = 29;
    blocklengths = malloc(nitems * sizeof(*blocklengths));
    check_alloc_status(blocklengths, "Memory allocation error.");

    offsets = malloc(nitems * sizeof(*offsets));
    check_alloc_status(offsets, "Memory allocation error.");

    mpi_types = malloc(nitems * sizeof(*mpi_types));
    check_alloc_status(mpi_types, "Memory allocation error.");

    // none of the elements in option_struct are arrays
    for (i = 0; i < (size_t) nitems; i++) {
        blocklengths[i] = 1;
    }

    // reset i
    i = 0;

    // unsigned short AERO_RESIST;
    offsets[i] = offsetof(option_struct, AERO_RESIST);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // bool CARBON;
    offsets[i] = offsetof(option_struct, CARBON);
    mpi_types[i++] = MPI_C_BOOL;

    // bool CONTINUEONERROR;
    offsets[i] = offsetof(option_struct, CONTINUEONERROR);
    mpi_types[i++] = MPI_C_BOOL;

    // bool CORRPREC;
    offsets[i] = offsetof(option_struct, CORRPREC);
    mpi_types[i++] = MPI_C_BOOL;

    // bool FROZEN_SOIL;
    offsets[i] = offsetof(option_struct, FROZEN_SOIL);
    mpi_types[i++] = MPI_C_BOOL;

    // size_t Nlayer;
    offsets[i] = offsetof(option_struct, Nlayer);
    mpi_types[i++] = MPI_AINT; // note there is no MPI_SIZE_T equivalent

    // size_t Nswband;
    offsets[i] = offsetof(option_struct, Nswband);
    mpi_types[i++] = MPI_AINT; // note there is no MPI_SIZE_T equivalent  

    // bool NOFLUX;
    offsets[i] = offsetof(option_struct, NOFLUX);
    mpi_types[i++] = MPI_C_BOOL;

    // size_t NVEGTYPES;
    offsets[i] = offsetof(option_struct, NVEGTYPES);
    mpi_types[i++] = MPI_AINT; // note there is no MPI_SIZE_T equivalent

    // unsigned short SNOW_DENSITY;
    offsets[i] = offsetof(option_struct, SNOW_DENSITY);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // size_t SNOW_BAND;
    offsets[i] = offsetof(option_struct, SNOW_BAND);
    mpi_types[i++] = MPI_AINT; // note there is no MPI_SIZE_T equivalent

    // size_t GLACIER_ID;
    offsets[i] = offsetof(option_struct, GLACIER_ID);
    mpi_types[i++] = MPI_AINT; // note there is no MPI_SIZE_T equivalent

    // bool TFALLBACK;
    offsets[i] = offsetof(option_struct, TFALLBACK);
    mpi_types[i++] = MPI_C_BOOL;

    // unsigned short INTERCEPTION;
    offsets[i] = offsetof(option_struct, CANOPY_INTERCEP);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short SNOW_AGING;
    offsets[i] = offsetof(option_struct, SNOW_AGING);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short GRID_DECIMAL;
    offsets[i] = offsetof(option_struct, GRID_DECIMAL);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // bool VEGLIB_FCAN;
    offsets[i] = offsetof(option_struct, VEGLIB_FCAN);
    mpi_types[i++] = MPI_C_BOOL;

    // bool VEGPARAM_FCAN;
    offsets[i] = offsetof(option_struct, VEGPARAM_FCAN);
    mpi_types[i++] = MPI_C_BOOL;

    // bool VEGPARAM_LAI;
    offsets[i] = offsetof(option_struct, VEGPARAM_LAI);
    mpi_types[i++] = MPI_C_BOOL;

    // unsigned short FCAN_SRC;
    offsets[i] = offsetof(option_struct, FCAN_SRC);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short LAI_SRC;
    offsets[i] = offsetof(option_struct, LAI_SRC);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // bool DENSITY_FROM_SOIL;
    offsets[i] = offsetof(option_struct, PARAM_FROM_SOIL);
    mpi_types[i++] = MPI_C_BOOL;

    // unsigned short STATE_FORMAT;
    offsets[i] = offsetof(option_struct, STATE_FORMAT);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // bool INIT_STATE;
    offsets[i] = offsetof(option_struct, INIT_STATE);
    mpi_types[i++] = MPI_C_BOOL;

    // bool SAVE_STATE;
    offsets[i] = offsetof(option_struct, SAVE_STATE);
    mpi_types[i++] = MPI_C_BOOL;

    // make sure that the we have the right number of elements
    if (i != (size_t) nitems) {
        log_err("Miscount: %zd not equal to %d.", i, nitems);
    }

    status = MPI_Type_create_struct(nitems, blocklengths, offsets, mpi_types,
                                    mpi_type);
    check_mpi_status(status, "MPI error.");
    status = MPI_Type_commit(mpi_type);
    check_mpi_status(status, "MPI error.");

    // cleanup
    free(blocklengths);
    free(offsets);
    free(mpi_types);
}

/******************************************************************************
 * @brief   Create an MPI_Datatype that represents the parameters_struct
 * @details This allows MPI operations in which the entire parameters_struct
 *          can be treated as an MPI_Datatype. NOTE: This function needs to be
 *          kept in-sync with the parameters_struct data type in vic_def.h.
 *
 * @param mpi_type MPI_Datatype that can be used in MPI operations
 *****************************************************************************/
void
create_MPI_param_struct_type(MPI_Datatype *mpi_type)
{
    extern MPI_Comm MPI_COMM_VIC;
    int             nitems; // number of elements in struct
    int             status;
    int            *blocklengths;
    size_t          i;
    MPI_Aint       *offsets;
    MPI_Datatype   *mpi_types;

    // nitems has to equal the number of elements in parameters_struct
    nitems = 186;
    blocklengths = malloc(nitems * sizeof(*blocklengths));
    check_alloc_status(blocklengths, "Memory allocation error.");

    offsets = malloc(nitems * sizeof(*offsets));
    check_alloc_status(offsets, "Memory allocation error.");

    mpi_types = malloc(nitems * sizeof(*mpi_types));
    check_alloc_status(mpi_types, "Memory allocation error.");

    // none of the elements in parameters_struct are arrays
    for (i = 0; i < (size_t) nitems; i++) {
        blocklengths[i] = 1;
    }

    // reset i
    i = 0;

    // double LAPSE_RATE;
    offsets[i] = offsetof(parameters_struct, LAPSE_RATE);
    mpi_types[i++] = MPI_DOUBLE;

    // double GAUGE_HEIGHT;
    offsets[i] = offsetof(parameters_struct, GAUGE_HEIGHT);
    mpi_types[i++] = MPI_DOUBLE;

    // double HUGE_RESIST;
    offsets[i] = offsetof(parameters_struct, HUGE_RESIST);
    mpi_types[i++] = MPI_DOUBLE;

    // double REF_HEIGHT;
    offsets[i] = offsetof(parameters_struct, REF_HEIGHT);
    mpi_types[i++] = MPI_DOUBLE;

    // double REF_HEIGHT_WIND;
    offsets[i] = offsetof(parameters_struct, REF_HEIGHT_WIND);
    mpi_types[i++] = MPI_DOUBLE;

    // double ALBEDO_BARE_SOIL;
    offsets[i] = offsetof(parameters_struct, ALBEDO_BARE_SOIL);
    mpi_types[i++] = MPI_DOUBLE;

    // double EMISS_GRND;
    offsets[i] = offsetof(parameters_struct, EMISS_GRND);
    mpi_types[i++] = MPI_DOUBLE;

    // double EMISS_VEG
    offsets[i] = offsetof(parameters_struct, EMISS_VEG);
    mpi_types[i++] = MPI_DOUBLE;

    // double EMISS_ICE
    offsets[i] = offsetof(parameters_struct, EMISS_ICE);
    mpi_types[i++] = MPI_DOUBLE;

    // double EMISS_SNOW
    offsets[i] = offsetof(parameters_struct, EMISS_SNOW);
    mpi_types[i++] = MPI_DOUBLE;

    // double EMISS_H2O
    offsets[i] = offsetof(parameters_struct, EMISS_H2O);
    mpi_types[i++] = MPI_DOUBLE;

    // double VEG_LAI_SNOW_MULTIPLIER
    offsets[i] = offsetof(parameters_struct, VEG_LAI_SNOW_MULTIPLIER);
    mpi_types[i++] = MPI_DOUBLE;

    // double VEG_LAI_WATER_FACTOR
    offsets[i] = offsetof(parameters_struct, VEG_LAI_WATER_FACTOR);
    mpi_types[i++] = MPI_DOUBLE;

    // double CANOPY_CLOSURE
    offsets[i] = offsetof(parameters_struct, CANOPY_CLOSURE);
    mpi_types[i++] = MPI_DOUBLE;

    // double CANOPY_RSMAX
    offsets[i] = offsetof(parameters_struct, CANOPY_RSMAX);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_A0
    offsets[i] = offsetof(parameters_struct, SVP_A0);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_A1
    offsets[i] = offsetof(parameters_struct, SVP_A1);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_A2
    offsets[i] = offsetof(parameters_struct, SVP_A2);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_A3
    offsets[i] = offsetof(parameters_struct, SVP_A3);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_A4
    offsets[i] = offsetof(parameters_struct, SVP_A4);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_A5
    offsets[i] = offsetof(parameters_struct, SVP_A5);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_A6
    offsets[i] = offsetof(parameters_struct, SVP_A6);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_A7
    offsets[i] = offsetof(parameters_struct, SVP_A7);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_A8
    offsets[i] = offsetof(parameters_struct, SVP_A8);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_B0
    offsets[i] = offsetof(parameters_struct, SVP_B0);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_B1
    offsets[i] = offsetof(parameters_struct, SVP_B1);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_B2
    offsets[i] = offsetof(parameters_struct, SVP_B2);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_B3
    offsets[i] = offsetof(parameters_struct, SVP_B3);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_B4
    offsets[i] = offsetof(parameters_struct, SVP_B4);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_B5
    offsets[i] = offsetof(parameters_struct, SVP_B5);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_B6
    offsets[i] = offsetof(parameters_struct, SVP_B6);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_B7
    offsets[i] = offsetof(parameters_struct, SVP_B7);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_B8
    offsets[i] = offsetof(parameters_struct, SVP_B8);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_C0
    offsets[i] = offsetof(parameters_struct, SVP_C0);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_C1
    offsets[i] = offsetof(parameters_struct, SVP_C1);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_C2
    offsets[i] = offsetof(parameters_struct, SVP_C2);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_C3
    offsets[i] = offsetof(parameters_struct, SVP_C3);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_C4
    offsets[i] = offsetof(parameters_struct, SVP_C4);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_C5
    offsets[i] = offsetof(parameters_struct, SVP_C5);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_C6
    offsets[i] = offsetof(parameters_struct, SVP_C6);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_C7
    offsets[i] = offsetof(parameters_struct, SVP_C7);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_C8
    offsets[i] = offsetof(parameters_struct, SVP_C8);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_D0
    offsets[i] = offsetof(parameters_struct, SVP_D0);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_D1
    offsets[i] = offsetof(parameters_struct, SVP_D1);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_D2
    offsets[i] = offsetof(parameters_struct, SVP_D2);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_D3
    offsets[i] = offsetof(parameters_struct, SVP_D3);
    mpi_types[i++] = MPI_DOUBLE;
    
    // double SVP_D4
    offsets[i] = offsetof(parameters_struct, SVP_D4);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_D5
    offsets[i] = offsetof(parameters_struct, SVP_D5);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_D6
    offsets[i] = offsetof(parameters_struct, SVP_D6);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_D7
    offsets[i] = offsetof(parameters_struct, SVP_D7);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_D8
    offsets[i] = offsetof(parameters_struct, SVP_D8);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_FRZ
    offsets[i] = offsetof(parameters_struct, SVP_FRZ);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_RDAIR
    offsets[i] = offsetof(parameters_struct, SVP_RDAIR);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_OMEGA
    offsets[i] = offsetof(parameters_struct, PHOTO_OMEGA);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_LAIMAX
    offsets[i] = offsetof(parameters_struct, PHOTO_LAIMAX);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_LAILIMIT
    offsets[i] = offsetof(parameters_struct, PHOTO_LAILIMIT);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_LAIMIN
    offsets[i] = offsetof(parameters_struct, PHOTO_LAIMIN);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_EPAR
    offsets[i] = offsetof(parameters_struct, PHOTO_EPAR);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_FCMAX
    offsets[i] = offsetof(parameters_struct, PHOTO_FCMAX);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_FCMIN
    offsets[i] = offsetof(parameters_struct, PHOTO_FCMIN);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_ALBSOIPARMIN
    offsets[i] = offsetof(parameters_struct, PHOTO_ALBSOIPARMIN);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_MINMAXETRANS
    offsets[i] = offsetof(parameters_struct, PHOTO_MINMAXETRANS);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_FCI1C3
    offsets[i] = offsetof(parameters_struct, PHOTO_FCI1C3);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_FCI1C4
    offsets[i] = offsetof(parameters_struct, PHOTO_FCI1C4);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_OX
    offsets[i] = offsetof(parameters_struct, PHOTO_OX);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_KCPHOTO_KC
    offsets[i] = offsetof(parameters_struct, PHOTO_KC);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_KO
    offsets[i] = offsetof(parameters_struct, PHOTO_KO);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_EC
    offsets[i] = offsetof(parameters_struct, PHOTO_EC);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_EO
    offsets[i] = offsetof(parameters_struct, PHOTO_EO);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_EV
    offsets[i] = offsetof(parameters_struct, PHOTO_EV);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_ER
    offsets[i] = offsetof(parameters_struct, PHOTO_ER);
    mpi_types[i++] = MPI_DOUBLE;

    // double ROUGH3
    offsets[i] = offsetof(parameters_struct, ROUGH3);
    mpi_types[i++] = MPI_DOUBLE;

    // double ROUGH_BETA
    offsets[i] = offsetof(parameters_struct, ROUGH_BETA);
    mpi_types[i++] = MPI_DOUBLE;

    // double ROUGH_NU
    offsets[i] = offsetof(parameters_struct, ROUGH_NU);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_ROUGH
    offsets[i] = offsetof(parameters_struct, SNOW_ROUGH);
    mpi_types[i++] = MPI_DOUBLE;

    // double SOIL_ROUGH
    offsets[i] = offsetof(parameters_struct, SOIL_ROUGH);
    mpi_types[i++] = MPI_DOUBLE;

    // double GLAC_ROUGH
    offsets[i] = offsetof(parameters_struct, GLAC_ROUGH);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_MAX_SURFACE_SWE
    offsets[i] = offsetof(parameters_struct, SNOW_MAX_SURFACE_SWE);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_LIQUID_WATER_CAPACITY
    offsets[i] = offsetof(parameters_struct, SNOW_LIQUID_WATER_CAPACITY);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_NEW_SNOW_DENSITY
    offsets[i] = offsetof(parameters_struct, SNOW_NEW_SNOW_DENSITY);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_MAX_LIQUID_FRAC
    offsets[i] = offsetof(parameters_struct, SNOW_MAX_LIQUID_FRAC);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_RELEASE_FAC
    offsets[i] = offsetof(parameters_struct, SNOW_RELEASE_FAC);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_NEW_SNOW_DENS_MAX
    offsets[i] = offsetof(parameters_struct, SNOW_NEW_SNOW_DENS_MAX);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_NEW_SNT_C1
    offsets[i] = offsetof(parameters_struct, SNOW_NEW_SNT_C1);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_NEW_SNT_C2
    offsets[i] = offsetof(parameters_struct, SNOW_NEW_SNT_C2);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_NEW_SNT_C3
    offsets[i] = offsetof(parameters_struct, SNOW_NEW_SNT_C3);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_NEW_BRAS_DENOM
    offsets[i] = offsetof(parameters_struct, SNOW_NEW_BRAS_DENOM);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_COMPACT_A
    offsets[i] = offsetof(parameters_struct, SNOW_COMPACT_A);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_COMPACT_B
    offsets[i] = offsetof(parameters_struct, SNOW_COMPACT_B);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_COMPACT_C
    offsets[i] = offsetof(parameters_struct, SNOW_COMPACT_C);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_COMPACT_P
    offsets[i] = offsetof(parameters_struct, SNOW_COMPACT_P);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_COMPACT_DM
    offsets[i] = offsetof(parameters_struct, SNOW_COMPACT_DM);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_COMPACT_ETA
    offsets[i] = offsetof(parameters_struct, SNOW_COMPACT_ETA);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_RADIUS_MIN
    offsets[i] = offsetof(parameters_struct, SNOW_RADIUS_MIN);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_RADIUS_MAX
    offsets[i] = offsetof(parameters_struct, SNOW_RADIUS_MAX);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_NEW_RADIUS
    offsets[i] = offsetof(parameters_struct, SNOW_NEW_RADIUS);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_AGE_FACT
    offsets[i] = offsetof(parameters_struct, SNOW_AGE_FACT);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_AGE_VAPF
    offsets[i] = offsetof(parameters_struct, SNOW_AGE_VAPF);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_AGE_FRZF
    offsets[i] = offsetof(parameters_struct, SNOW_AGE_FRZF);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_AGE_SOTF
    offsets[i] = offsetof(parameters_struct, SNOW_AGE_SOTF);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_NEW_SNOW_COVER
    offsets[i] = offsetof(parameters_struct, SNOW_NEW_SNOW_COVER);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_COSZEN_B
    offsets[i] = offsetof(parameters_struct, SNOW_COSZEN_B);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_AGE_DIR_VIS
    offsets[i] = offsetof(parameters_struct, SNOW_AGE_DIR_VIS);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_AGE_DIR_NIR
    offsets[i] = offsetof(parameters_struct, SNOW_AGE_DIR_NIR);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_AGE_DFS_VIS
    offsets[i] = offsetof(parameters_struct, SNOW_AGE_DFS_VIS);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_AGE_DFS_NIR
    offsets[i] = offsetof(parameters_struct, SNOW_AGE_DFS_NIR);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_NEW_SNOW_VIS
    offsets[i] = offsetof(parameters_struct, SNOW_NEW_SNOW_VIS);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_NEW_SNOW_NIR
    offsets[i] = offsetof(parameters_struct, SNOW_NEW_SNOW_NIR);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_CONDUCT
    offsets[i] = offsetof(parameters_struct, SNOW_CONDUCT);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_BETADS
    offsets[i] = offsetof(parameters_struct, SNOW_BETADS);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_BETAIS
    offsets[i] = offsetof(parameters_struct, SNOW_BETAIS);
    mpi_types[i++] = MPI_DOUBLE;

    // double RAD_DIR_F
    offsets[i] = offsetof(parameters_struct, RAD_DIR_F);
    mpi_types[i++] = MPI_DOUBLE;

    // double RAD_VIS_F
    offsets[i] = offsetof(parameters_struct, RAD_VIS_F);
    mpi_types[i++] = MPI_DOUBLE;

    // double TOL_A
    offsets[i] = offsetof(parameters_struct, TOL_A);
    mpi_types[i++] = MPI_DOUBLE;

    // double TOL_B
    offsets[i] = offsetof(parameters_struct, TOL_B);
    mpi_types[i++] = MPI_DOUBLE;

    // double MAX_LIMIT
    offsets[i] = offsetof(parameters_struct, MAX_LIMIT);
    mpi_types[i++] = MPI_DOUBLE;

    // double TOL_WETBULB
    offsets[i] = offsetof(parameters_struct, TOL_WETBULB);
    mpi_types[i++] = MPI_DOUBLE;

    // double CN_FACTOR
    offsets[i] = offsetof(parameters_struct, CN_FACTOR);
    mpi_types[i++] = MPI_DOUBLE;

    // int MAX_ITER_WETBULB
    offsets[i] = offsetof(parameters_struct, MAX_ITER_WETBULB);
    mpi_types[i++] = MPI_INT;

    // int MAX_ITER_OVER
    offsets[i] = offsetof(parameters_struct, MAX_ITER_OVER);
    mpi_types[i++] = MPI_INT;

    // int MAX_ITER_MOST
    offsets[i] = offsetof(parameters_struct, MAX_ITER_MOST);
    mpi_types[i++] = MPI_INT;

    // double SNOW_OMEGAS[MAX_SWBANDS]
    blocklengths[i] = MAX_SWBANDS;
    offsets[i] = offsetof(parameters_struct, SNOW_OMEGAS);
    mpi_types[i++] = MPI_DOUBLE;

    // double GLAC_ALBEDO[MAX_SWBANDS]
    blocklengths[i] = MAX_SWBANDS;
    offsets[i] = offsetof(parameters_struct, GLAC_ALBEDO);
    mpi_types[i++] = MPI_DOUBLE;

    // make sure that the we have the right number of elements
    if (i != (size_t) nitems) {
        log_err("Miscount: %zd not equal to %d.", i, nitems);
    }

    status = MPI_Type_create_struct(nitems, blocklengths, offsets, mpi_types,
                                    mpi_type);
    check_mpi_status(status, "MPI error.");
    status = MPI_Type_commit(mpi_type);
    check_mpi_status(status, "MPI error.");

    // cleanup
    free(blocklengths);
    free(offsets);
    free(mpi_types);
}

/******************************************************************************
 * @brief   Create an MPI_Datatype that represents the dmy_struct
 * @details This allows MPI operations in which the entire dmy_struct
 *          can be treated as an MPI_Datatype.
 * @param mpi_type MPI_Datatype that can be used in MPI operations
 *****************************************************************************/
void
create_MPI_dmy_struct_type(MPI_Datatype *mpi_type)
{
    extern MPI_Comm MPI_COMM_VIC;

    int             nitems; // number of elements in struct
    int             status;
    int            *blocklengths;
    size_t          i;
    MPI_Aint       *offsets;
    MPI_Datatype   *mpi_types;

    // nitems has to equal the number of elements in dmy_struct
    nitems = 5;
    blocklengths = malloc(nitems * sizeof(*blocklengths));
    check_alloc_status(blocklengths, "Memory allocation error.");

    offsets = malloc(nitems * sizeof(*offsets));
    check_alloc_status(offsets, "Memory allocation error.");

    mpi_types = malloc(nitems * sizeof(*mpi_types));
    check_alloc_status(mpi_types, "Memory allocation error.");

    // none of the elements in location_struct are arrays.
    for (i = 0; i < (size_t) nitems; i++) {
        blocklengths[i] = 1;
    }

    // reset i
    i = 0;

    // unsigned short int day;
    offsets[i] = offsetof(dmy_struct, day);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short int day_in_year;
    offsets[i] = offsetof(dmy_struct, day_in_year);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short int month;
    offsets[i] = offsetof(dmy_struct, month);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // int year;
    offsets[i] = offsetof(dmy_struct, year);
    mpi_types[i++] = MPI_INT;

    // unsigned int dayseconds;
    offsets[i] = offsetof(dmy_struct, dayseconds);
    mpi_types[i++] = MPI_UNSIGNED;

    // make sure that the we have the right number of elements
    if (i != (size_t) nitems) {
        log_err("Miscount: %zd not equal to %d.", i, nitems);
    }

    status = MPI_Type_create_struct(nitems, blocklengths, offsets, mpi_types,
                                    mpi_type);
    check_mpi_status(status, "MPI error.");

    status = MPI_Type_commit(mpi_type);
    check_mpi_status(status, "MPI error.");

    // cleanup
    free(blocklengths);
    free(offsets);
    free(mpi_types);
}

/******************************************************************************
 * @brief   Create an MPI_Datatype that represents the alarm_struct
 * @details This allows MPI operations in which the entire alarm_struct
 *          can be treated as an MPI_Datatype.
 * @param mpi_type MPI_Datatype that can be used in MPI operations
 *****************************************************************************/
void
create_MPI_alarm_struct_type(MPI_Datatype *mpi_type)
{
    extern MPI_Comm MPI_COMM_VIC;

    int             nitems; // number of elements in struct
    int             status;
    int            *blocklengths;
    size_t          i;
    MPI_Aint       *offsets;
    MPI_Datatype   *mpi_types;
    MPI_Datatype    mpi_dmy_type;

    // nitems has to equal the number of elements in alarm_struct
    nitems = 6;
    blocklengths = malloc(nitems * sizeof(*blocklengths));
    check_alloc_status(blocklengths, "Memory allocation error.");

    offsets = malloc(nitems * sizeof(*offsets));
    check_alloc_status(offsets, "Memory allocation error.");

    mpi_types = malloc(nitems * sizeof(*mpi_types));
    check_alloc_status(mpi_types, "Memory allocation error.");

    // none of the elements in location_struct are arrays.
    for (i = 0; i < (size_t) nitems; i++) {
        blocklengths[i] = 1;
    }

    // reset i
    i = 0;

    // unsigned int count;
    offsets[i] = offsetof(alarm_struct, count);
    mpi_types[i++] = MPI_UNSIGNED;

    // int next_count;
    offsets[i] = offsetof(alarm_struct, next_count);
    mpi_types[i++] = MPI_INT;

    // dmy_struct next_dmy;
    offsets[i] = offsetof(alarm_struct, next_dmy);
    create_MPI_dmy_struct_type(&mpi_dmy_type);
    mpi_types[i++] = mpi_dmy_type;

    // unsigned int freq;
    offsets[i] = offsetof(alarm_struct, freq);
    mpi_types[i++] = MPI_UNSIGNED;

    // int n;
    offsets[i] = offsetof(alarm_struct, n);
    mpi_types[i++] = MPI_INT;

    // bool is_subdaily;
    offsets[i] = offsetof(alarm_struct, is_subdaily);
    mpi_types[i++] = MPI_C_BOOL;

    // make sure that the we have the right number of elements
    if (i != (size_t) nitems) {
        log_err("Miscount: %zd not equal to %d.", i, nitems);
    }

    status = MPI_Type_create_struct(nitems, blocklengths, offsets, mpi_types,
                                    mpi_type);
    check_mpi_status(status, "MPI error.");

    status = MPI_Type_commit(mpi_type);
    check_mpi_status(status, "MPI error.");

    // cleanup
    free(blocklengths);
    free(offsets);
    free(mpi_types);
    MPI_Type_free(&mpi_dmy_type);
}

/******************************************************************************
 * @brief   Type-agnostic mapping function
 * @details Reorders the elements in 'from' to 'to' according to the ordering
 *          specified in 'map'.
 *          Note that this function can also be used for filtering, i.e. you
 *          can use a smaller number of elements in 'map' and 'to' than in
 *          'from' to get only a subset of the elements.
 *
 *          to[to_map[i]] = from[from_map[i]]
 *
 * @param size size of the datatype of 'from' and 'to', e.g. sizeof(int)
 * @param n number of elements in 'from_map' and 'to_map'
 * @param from_map array of length n with 'from' indices, if from_map == NULL,
 *        then the 'from' indices are sequential
 * @param to_map array of length n with 'to' indices, if to_map == NULL, then
 *        the 'to' indices are sequential
 * @param from array of with entries of size 'size' (unchanged)
 * @param to array of with entries of size 'size' (changed)
 *****************************************************************************/
void
map(size_t  size,
    size_t  n,
    size_t *from_map,
    size_t *to_map,
    void   *from,
    void   *to)
{
    size_t i;

    if (to_map == NULL && from_map == NULL) {
        for (i = 0; i < n; i++) {
            // type-agnostic version of to[i] = from[i];
            memcpy((void *)((char *)to + i * size),
                   (void *)((char *)from + i * size), size);
        }
    }
    if (to_map == NULL) {
        for (i = 0; i < n; i++) {
            // type-agnostic version of to[i] = from[from_map[i]];
            memcpy((void *)((char *)to + i * size),
                   (void *)((char *)from + from_map[i] * size), size);
        }
    }
    else if (from_map == NULL) {
        for (i = 0; i < n; i++) {
            // type-agnostic version of to[to_map[i]] = from[i];
            memcpy((void *)((char *)to + to_map[i] * size),
                   (void *)((char *)from + i * size), size);
        }
    }
    else {
        for (i = 0; i < n; i++) {
            // type-agnostic version of to[to_map[i]] = from[from_map[i]];
            memcpy((void *)((char *)to + to_map[i] * size),
                   (void *)((char *)from + from_map[i] * size), size);
        }
    }
}

/******************************************************************************
 * @brief   Decompose the domain for MPI operations
 * @details This function sets up the arrays needed to scatter and gather
 *          data from and to the master process to the individual mpi
 *          processes.
 *
 * @param ncells total number of cells
 * @param mpi_size number of mpi processes
 * @param mpi_map_local_array_sizes address of integer array with number of
 *        cells assigned to each node (MPI_Scatterv:sendcounts and
 *        MPI_Gatherv:recvcounts)
 * @param mpi_map_global_array_offsets address of integer array with offsets
 *        for sending and receiving data (MPI_Scatterv:displs and
 *        MPI_Gatherv:displs)
 * @param mpi_map_mapping_array address of size_t array with indices to prepare
 *        an array on the master process for MPI_Scatterv or map back after
 *        MPI_Gatherv
 *****************************************************************************/
void
mpi_map_decomp_domain(size_t   ncells,
                      size_t   mpi_size,
                      int    **mpi_map_local_array_sizes,
                      int    **mpi_map_global_array_offsets,
                      size_t **mpi_map_mapping_array)
{
    size_t i;
    size_t j;
    size_t k;
    size_t n;

    *mpi_map_local_array_sizes = calloc(mpi_size,
                                        sizeof(*(*mpi_map_local_array_sizes)));
    *mpi_map_global_array_offsets = calloc(mpi_size,
                                           sizeof(*(*
                                                    mpi_map_global_array_offsets)));
    *mpi_map_mapping_array = calloc(ncells, sizeof(*(*mpi_map_mapping_array)));

    // determine number of cells per node
    for (n = ncells, i = 0; n > 0; n--, i++) {
        if (i >= mpi_size) {
            i = 0;
        }
        (*mpi_map_local_array_sizes)[i] += 1;
    }

    // determine offsets to use for MPI_Scatterv and MPI_Gatherv
    for (i = 1; i < mpi_size; i++) {
        for (j = 0; j < i; j++) {
            (*mpi_map_global_array_offsets)[i] +=
                (*mpi_map_local_array_sizes)[j];
        }
    }

    // set mapping array
    for (i = 0, k = 0; i < (size_t) mpi_size; i++) {
        for (j = 0; j < (size_t) (*mpi_map_local_array_sizes)[i]; j++) {
            (*mpi_map_mapping_array)[k++] = (size_t) (i + j * mpi_size);
        }
    }
}

/******************************************************************************
 * @brief   Gather double precision variable
 * @details Values are gathered to the master node
 *****************************************************************************/
void
gather_field_double(double  fillval,
                    double *dvar,
                    double *var)
{
    extern MPI_Comm      MPI_COMM_VIC;
    extern domain_struct global_domain;
    extern domain_struct local_domain;
    extern int           mpi_rank;
    extern int          *mpi_map_global_array_offsets;
    extern int          *mpi_map_local_array_sizes;
    extern size_t       *filter_active_cells;
    extern size_t       *mpi_map_mapping_array;
    int                  status;
    double              *dvar_gathered = NULL;
    double              *dvar_remapped = NULL;
    size_t               grid_size;
    size_t               i;

    if (mpi_rank == VIC_MPI_ROOT) {
        grid_size = global_domain.n_nx * global_domain.n_ny;
        for (i = 0; i < grid_size; i++) {
            dvar[i] = fillval;
        }
        dvar_gathered =
            malloc(global_domain.ncells_active * sizeof(*dvar_gathered));
        check_alloc_status(dvar_gathered, "Memory allocation error.");

        dvar_remapped =
            malloc(global_domain.ncells_active * sizeof(*dvar_remapped));
        check_alloc_status(dvar_remapped, "Memory allocation error.");
    }
    // Gather the results from the nodes, result for the local node is in the
    // array *var (which is a function argument)
    status = MPI_Gatherv(var, local_domain.ncells_active, MPI_DOUBLE,
                         dvar_gathered, mpi_map_local_array_sizes,
                         mpi_map_global_array_offsets, MPI_DOUBLE,
                         VIC_MPI_ROOT, MPI_COMM_VIC);
    check_mpi_status(status, "MPI error.");
    if (mpi_rank == VIC_MPI_ROOT) {
        // remap the array
        map(sizeof(double), global_domain.ncells_active, NULL,
            mpi_map_mapping_array, dvar_gathered, dvar_remapped);
        // expand to full grid size
        map(sizeof(double), global_domain.ncells_active, NULL,
            filter_active_cells, dvar_remapped, dvar);
        // cleanup
        free(dvar_gathered);
        free(dvar_remapped);
    }
}

/******************************************************************************
 * @brief   Gather and write double precision NetCDF field
 * @details Values are gathered to the master node and then written from the
 *          master node
 *****************************************************************************/
void
gather_put_nc_field_double(int     nc_id,
                           int     var_id,
                           double  fillval,
                           size_t *start,
                           size_t *count,
                           double *var)
{
    extern int           mpi_rank;
    extern domain_struct global_domain;
    int                  status;
    size_t               grid_size;
    double              *dvar = NULL;

    // Allocate memory
    if (mpi_rank == VIC_MPI_ROOT) {
        grid_size = global_domain.n_nx * global_domain.n_ny;
        dvar = malloc(grid_size * sizeof(*dvar));
        check_alloc_status(dvar, "Memory allocation error.");
    }

    // Gather results from the nodes
    gather_field_double(fillval, dvar, var);

    // Write to netcdf
    if (mpi_rank == VIC_MPI_ROOT) {
        status = nc_put_vara_double(nc_id, var_id, start, count, dvar);
        check_nc_status(status, "Error writing values.");
        // cleanup
        free(dvar);
    }
}

/******************************************************************************
 * @brief   Gather and write double precision NetCDF field
 * @details Values are gathered to the master node and then written from the
 *          master node
 *****************************************************************************/
void
gather_put_nc_field_float(int     nc_id,
                          int     var_id,
                          float   fillval,
                          size_t *start,
                          size_t *count,
                          float  *var)
{
    extern MPI_Comm      MPI_COMM_VIC;
    extern domain_struct global_domain;
    extern domain_struct local_domain;
    extern int           mpi_rank;
    extern int          *mpi_map_global_array_offsets;
    extern int          *mpi_map_local_array_sizes;
    extern size_t       *filter_active_cells;
    extern size_t       *mpi_map_mapping_array;
    int                  status;
    float               *fvar = NULL;
    float               *fvar_gathered = NULL;
    float               *fvar_remapped = NULL;
    size_t               grid_size;
    size_t               i;

    if (mpi_rank == VIC_MPI_ROOT) {
        grid_size = global_domain.n_nx * global_domain.n_ny;
        fvar = malloc(grid_size * sizeof(*fvar));
        check_alloc_status(fvar, "Memory allocation error.");
        for (i = 0; i < grid_size; i++) {
            fvar[i] = fillval;
        }
        fvar_gathered =
            malloc(global_domain.ncells_active * sizeof(*fvar_gathered));
        check_alloc_status(fvar_gathered, "Memory allocation error.");

        fvar_remapped =
            malloc(global_domain.ncells_active * sizeof(*fvar_remapped));
        check_alloc_status(fvar_remapped, "Memory allocation error.");
    }
    // Gather the results from the nodes, result for the local node is in the
    // array *var (which is a function argument)
    status = MPI_Gatherv(var, local_domain.ncells_active, MPI_FLOAT,
                         fvar_gathered, mpi_map_local_array_sizes,
                         mpi_map_global_array_offsets, MPI_FLOAT,
                         VIC_MPI_ROOT, MPI_COMM_VIC);
    check_mpi_status(status, "Error with gather of floats");

    if (mpi_rank == VIC_MPI_ROOT) {
        // remap the array
        map(sizeof(float), global_domain.ncells_active, NULL,
            mpi_map_mapping_array, fvar_gathered, fvar_remapped);
        // expand to full grid size
        map(sizeof(float), global_domain.ncells_active, NULL,
            filter_active_cells, fvar_remapped, fvar);

        // write to file
        status = nc_put_vara_float(nc_id, var_id, start, count, fvar);
        check_nc_status(status, "Error writing values");

        // cleanup
        free(fvar);
        free(fvar_gathered);
        free(fvar_remapped);
    }
}

/******************************************************************************
 * @brief   Gather and write integer NetCDF field
 * @details Values are gathered to the master node and then written from the
 *          master node
 *****************************************************************************/
void
gather_put_nc_field_int(int     nc_id,
                        int     var_id,
                        int     fillval,
                        size_t *start,
                        size_t *count,
                        int    *var)
{
    extern MPI_Comm      MPI_COMM_VIC;
    extern domain_struct global_domain;
    extern domain_struct local_domain;
    extern int           mpi_rank;
    extern int          *mpi_map_global_array_offsets;
    extern int          *mpi_map_local_array_sizes;
    extern size_t       *filter_active_cells;
    extern size_t       *mpi_map_mapping_array;
    int                  status;
    int                 *ivar = NULL;
    int                 *ivar_gathered = NULL;
    int                 *ivar_remapped = NULL;
    size_t               grid_size;
    size_t               i;

    if (mpi_rank == VIC_MPI_ROOT) {
        grid_size = global_domain.n_nx * global_domain.n_ny;
        ivar = malloc(grid_size * sizeof(*ivar));
        check_alloc_status(ivar, "Memory allocation error.");

        for (i = 0; i < grid_size; i++) {
            ivar[i] = fillval;
        }

        ivar_gathered =
            malloc(global_domain.ncells_active * sizeof(*ivar_gathered));
        check_alloc_status(ivar_gathered, "Memory allocation error.");


        ivar_remapped =
            malloc(global_domain.ncells_active * sizeof(*ivar_remapped));
        check_alloc_status(ivar_remapped, "Memory allocation error.");
    }
    // Gather the results from the nodes, result for the local node is in the
    // array *var (which is a function argument)
    status = MPI_Gatherv(var, local_domain.ncells_active, MPI_INT,
                         ivar_gathered, mpi_map_local_array_sizes,
                         mpi_map_global_array_offsets, MPI_INT,
                         VIC_MPI_ROOT, MPI_COMM_VIC);
    check_mpi_status(status, "MPI error.");

    if (mpi_rank == VIC_MPI_ROOT) {
        // remap the array
        map(sizeof(int), global_domain.ncells_active, NULL,
            mpi_map_mapping_array,
            ivar_gathered, ivar_remapped);
        // expand to full grid size
        map(sizeof(int), global_domain.ncells_active, NULL, filter_active_cells,
            ivar_remapped, ivar);
        // write to file
        status = nc_put_vara_int(nc_id, var_id, start, count, ivar);
        check_nc_status(status, "Error writing values");

        // cleanup
        free(ivar);
        free(ivar_gathered);
        free(ivar_remapped);
    }
}

/******************************************************************************
 * @brief   Gather and write short integer NetCDF field
 * @details Values are gathered to the master node and then written from the
 *          master node
 *****************************************************************************/
void
gather_put_nc_field_short(int        nc_id,
                          int        var_id,
                          short int  fillval,
                          size_t    *start,
                          size_t    *count,
                          short int *var)
{
    extern MPI_Comm      MPI_COMM_VIC;
    extern domain_struct global_domain;
    extern domain_struct local_domain;
    extern int           mpi_rank;
    extern int          *mpi_map_global_array_offsets;
    extern int          *mpi_map_local_array_sizes;
    extern size_t       *filter_active_cells;
    extern size_t       *mpi_map_mapping_array;
    int                  status;
    short int           *svar = NULL;
    short int           *svar_gathered = NULL;
    short int           *svar_remapped = NULL;
    size_t               grid_size;
    size_t               i;

    if (mpi_rank == VIC_MPI_ROOT) {
        grid_size = global_domain.n_nx * global_domain.n_ny;
        svar = malloc(grid_size * sizeof(*svar));
        check_alloc_status(svar, "Memory allocation error.");

        for (i = 0; i < grid_size; i++) {
            svar[i] = fillval;
        }

        svar_gathered =
            malloc(global_domain.ncells_active * sizeof(*svar_gathered));
        check_alloc_status(svar_gathered, "Memory allocation error.");


        svar_remapped =
            malloc(global_domain.ncells_active * sizeof(*svar_remapped));
        check_alloc_status(svar_remapped, "Memory allocation error.");
    }
    // Gather the results from the nodes, result for the local node is in the
    // array *var (which is a function argument)
    status = MPI_Gatherv(var, local_domain.ncells_active, MPI_SHORT,
                         svar_gathered, mpi_map_local_array_sizes,
                         mpi_map_global_array_offsets, MPI_SHORT,
                         VIC_MPI_ROOT, MPI_COMM_VIC);
    check_mpi_status(status, "MPI error.");

    if (mpi_rank == VIC_MPI_ROOT) {
        // remap the array
        map(sizeof(short int), global_domain.ncells_active, NULL,
            mpi_map_mapping_array,
            svar_gathered, svar_remapped);
        // expand to full grid size
        map(sizeof(short int), global_domain.ncells_active, NULL,
            filter_active_cells, svar_remapped, svar);
        // write to file
        status = nc_put_vara_short(nc_id, var_id, start, count, svar);
        check_nc_status(status, "Error writing values");

        // cleanup
        free(svar);
        free(svar_gathered);
        free(svar_remapped);
    }
}

/******************************************************************************
 * @brief   Gather and write signed character NetCDF field
 * @details Values are gathered to the master node and then written from the
 *          master node
 *****************************************************************************/
void
gather_put_nc_field_schar(int     nc_id,
                          int     var_id,
                          char    fillval,
                          size_t *start,
                          size_t *count,
                          char   *var)
{
    extern MPI_Comm      MPI_COMM_VIC;
    extern domain_struct global_domain;
    extern domain_struct local_domain;
    extern int           mpi_rank;
    extern int          *mpi_map_global_array_offsets;
    extern int          *mpi_map_local_array_sizes;
    extern size_t       *filter_active_cells;
    extern size_t       *mpi_map_mapping_array;
    int                  status;
    signed char         *cvar = NULL;
    signed char         *cvar_gathered = NULL;
    signed char         *cvar_remapped = NULL;
    size_t               grid_size;
    size_t               i;

    if (mpi_rank == VIC_MPI_ROOT) {
        grid_size = global_domain.n_nx * global_domain.n_ny;
        cvar = malloc(grid_size * sizeof(*cvar));
        check_alloc_status(cvar, "Memory allocation error.");

        for (i = 0; i < grid_size; i++) {
            cvar[i] = fillval;
        }

        cvar_gathered =
            malloc(global_domain.ncells_active * sizeof(*cvar_gathered));
        check_alloc_status(cvar_gathered, "Memory allocation error.");


        cvar_remapped =
            malloc(global_domain.ncells_active * sizeof(*cvar_remapped));
        check_alloc_status(cvar_remapped, "Memory allocation error.");
    }
    // Gather the results from the nodes, result for the local node is in the
    // array *var (which is a function argument)
    status = MPI_Gatherv(var, local_domain.ncells_active, MPI_CHAR,
                         cvar_gathered, mpi_map_local_array_sizes,
                         mpi_map_global_array_offsets, MPI_CHAR,
                         VIC_MPI_ROOT, MPI_COMM_VIC);
    check_mpi_status(status, "MPI error.");

    if (mpi_rank == VIC_MPI_ROOT) {
        // remap the array
        map(sizeof(char), global_domain.ncells_active, NULL,
            mpi_map_mapping_array,
            cvar_gathered, cvar_remapped);
        // expand to full grid size
        map(sizeof(char), global_domain.ncells_active, NULL,
            filter_active_cells, cvar_remapped, cvar);
        // write to file
        status = nc_put_vara_schar(nc_id, var_id, start, count, cvar);
        check_nc_status(status, "Error writing values");

        // cleanup
        free(cvar);
        free(cvar_gathered);
        free(cvar_remapped);
    }
}

/******************************************************************************
 * @brief   Scatter double precision variable
 * @details values from master node are scattered to the local nodes
 *****************************************************************************/
void
scatter_field_double(double *dvar,
                     double *var)
{
    extern MPI_Comm      MPI_COMM_VIC;
    extern domain_struct global_domain;
    extern domain_struct local_domain;
    extern int           mpi_rank;
    extern int          *mpi_map_global_array_offsets;
    extern int          *mpi_map_local_array_sizes;
    extern size_t       *filter_active_cells;
    extern size_t       *mpi_map_mapping_array;
    int                  status;
    double              *dvar_filtered = NULL;
    double              *dvar_mapped = NULL;

    if (mpi_rank == VIC_MPI_ROOT) {
        dvar_filtered =
            malloc(global_domain.ncells_active * sizeof(*dvar_filtered));
        check_alloc_status(dvar_filtered, "Memory allocation error.");

        dvar_mapped =
            malloc(global_domain.ncells_active * sizeof(*dvar_mapped));
        check_alloc_status(dvar_mapped, "Memory allocation error.");

        // filter the active cells only
        map(sizeof(double), global_domain.ncells_active, filter_active_cells,
            NULL, dvar, dvar_filtered);
        // map to prepare for MPI_Scatterv
        map(sizeof(double), global_domain.ncells_active, mpi_map_mapping_array,
            NULL, dvar_filtered, dvar_mapped);
        free(dvar);
        free(dvar_filtered);
    }

    // Scatter the results to the nodes, result for the local node is in the
    // array *var (which is a function argument)
    status = MPI_Scatterv(dvar_mapped, mpi_map_local_array_sizes,
                          mpi_map_global_array_offsets, MPI_DOUBLE,
                          var, local_domain.ncells_active, MPI_DOUBLE,
                          VIC_MPI_ROOT, MPI_COMM_VIC);
    check_mpi_status(status, "MPI error.");

    if (mpi_rank == VIC_MPI_ROOT) {
        free(dvar_mapped);
    }
}

/******************************************************************************
 * @brief   Read double precision NetCDF field from file and scatter
 * @details Read happens on the master node and is then scattered to the local
 *          nodes
 *****************************************************************************/
void
get_scatter_nc_field_double(nameid_struct *nc_nameid,
                            char          *var_name,
                            size_t        *start,
                            size_t        *count,
                            double        *var)
{
    extern domain_struct global_domain;
    extern int           mpi_rank;
    double              *dvar = NULL;

    // Read variable from netcdf
    if (mpi_rank == VIC_MPI_ROOT) {
        dvar = malloc(global_domain.ncells_total * sizeof(*dvar));
        check_alloc_status(dvar, "Memory allocation error.");

        get_nc_field_double(nc_nameid, var_name, start, count, dvar);
    }

    // Scatter results to nodes
    scatter_field_double(dvar, var);
}

/******************************************************************************
 * @brief   Read single precision NetCDF field from file and scatter
 * @details Read happens on the master node and is then scattered to the local
 *          nodes
 *****************************************************************************/
void
get_scatter_nc_field_float(nameid_struct *nc_nameid,
                           char          *var_name,
                           size_t        *start,
                           size_t        *count,
                           float         *var)
{
    extern MPI_Comm      MPI_COMM_VIC;
    extern domain_struct global_domain;
    extern domain_struct local_domain;
    extern int           mpi_rank;
    extern int          *mpi_map_global_array_offsets;
    extern int          *mpi_map_local_array_sizes;
    extern size_t       *filter_active_cells;
    extern size_t       *mpi_map_mapping_array;
    int                  status;
    float               *fvar = NULL;
    float               *fvar_filtered = NULL;
    float               *fvar_mapped = NULL;

    if (mpi_rank == VIC_MPI_ROOT) {
        fvar = malloc(global_domain.ncells_total * sizeof(*fvar));
        check_alloc_status(fvar, "Memory allocation error.");

        fvar_filtered =
            malloc(global_domain.ncells_active * sizeof(*fvar_filtered));
        check_alloc_status(fvar_filtered, "Memory allocation error.");

        fvar_mapped =
            malloc(global_domain.ncells_active * sizeof(*fvar_mapped));
        check_alloc_status(fvar_mapped, "Memory allocation error.");

        get_nc_field_float(nc_nameid, var_name, start, count, fvar);
        // filter the active cells only
        map(sizeof(float), global_domain.ncells_active, filter_active_cells,
            NULL,
            fvar, fvar_filtered);
        // map to prepare for MPI_Scatterv
        map(sizeof(float), global_domain.ncells_active, mpi_map_mapping_array,
            NULL,
            fvar_filtered, fvar_mapped);
        free(fvar);
        free(fvar_filtered);
    }

    // Scatter the results to the nodes, result for the local node is in the
    // array *var (which is a function argument)
    status = MPI_Scatterv(fvar_mapped, mpi_map_local_array_sizes,
                          mpi_map_global_array_offsets, MPI_FLOAT,
                          var, local_domain.ncells_active, MPI_FLOAT,
                          VIC_MPI_ROOT, MPI_COMM_VIC);
    check_mpi_status(status, "MPI error.");

    if (mpi_rank == VIC_MPI_ROOT) {
        free(fvar_mapped);
    }
}

/******************************************************************************
 * @brief   Read integer NetCDF field from file and scatter
 * @details Read happens on the master node and is then scattered to the local
 *          nodes
 *****************************************************************************/
void
get_scatter_nc_field_int(nameid_struct *nc_nameid,
                         char          *var_name,
                         size_t        *start,
                         size_t        *count,
                         int           *var)
{
    extern MPI_Comm      MPI_COMM_VIC;
    extern domain_struct global_domain;
    extern domain_struct local_domain;
    extern int           mpi_rank;
    extern int          *mpi_map_global_array_offsets;
    extern int          *mpi_map_local_array_sizes;
    extern size_t       *filter_active_cells;
    extern size_t       *mpi_map_mapping_array;
    int                  status;
    int                 *ivar = NULL;
    int                 *ivar_filtered = NULL;
    int                 *ivar_mapped = NULL;

    if (mpi_rank == VIC_MPI_ROOT) {
        ivar = malloc(global_domain.ncells_total * sizeof(*ivar));
        check_alloc_status(ivar, "Memory allocation error.");
        ivar_filtered =
            malloc(global_domain.ncells_active * sizeof(*ivar_filtered));
        check_alloc_status(ivar_filtered, "Memory allocation error.");

        ivar_mapped =
            malloc(global_domain.ncells_active * sizeof(*ivar_mapped));
        check_alloc_status(ivar_mapped, "Memory allocation error.");

        get_nc_field_int(nc_nameid, var_name, start, count, ivar);
        // filter the active cells only
        map(sizeof(int), global_domain.ncells_active, filter_active_cells, NULL,
            ivar, ivar_filtered);
        // map to prepare for MPI_Scatterv
        map(sizeof(int), global_domain.ncells_active, mpi_map_mapping_array,
            NULL,
            ivar_filtered, ivar_mapped);
        free(ivar);
        free(ivar_filtered);
    }

    // Scatter the results to the nodes, result for the local node is in the
    // array *var (which is a function argument)
    status = MPI_Scatterv(ivar_mapped, mpi_map_local_array_sizes,
                          mpi_map_global_array_offsets, MPI_INT,
                          var, local_domain.ncells_active, MPI_INT,
                          VIC_MPI_ROOT, MPI_COMM_VIC);
    check_mpi_status(status, "MPI error.");

    if (mpi_rank == VIC_MPI_ROOT) {
        free(ivar_mapped);
    }
}
