/******************************************************************************
 * @section DESCRIPTION
 *
 * Initialize model parameters
 *****************************************************************************/

#include "vic_driver_classic.h"

/******************************************************************************
 * @brief    Initialize grid cell model data structures.
 *****************************************************************************/
void
vic_init_output(size_t            Ncells,
                grid_cell_struct *grid,
                filep_struct     *filep,
                dmy_struct       *dmy,
                veg_lib_struct   *veg_lib)
{
    extern option_struct    options;
    extern double        ***out_data;
    extern MPI_Datatype     mpi_alarm_struct_type;
    extern MPI_Comm         MPI_COMM_VIC;
    extern int              mpi_rank;
    extern stream_struct   *streams;
    extern nc_file_struct  *nc_hist_files;
    
    int                    status;
    size_t                 i;
    size_t                 streamnum;
    size_t                 nstream_vars[MAX_OUTPUT_STREAMS];
    bool                   default_outputs = false;
    timer_struct           timer;
    /** Set up output data structures **/
    set_output_met_data_info();

    // Allocate memory for out_data
    out_data = malloc(Ncells * sizeof(*out_data));
    check_alloc_status(out_data, "Memory allocation error.");
    // out_data is shape [ngridcells (1), N_OUTVAR_TYPES]
    alloc_out_data(Ncells, out_data);

    // initialize the save data structures
    for (i = 0; i < Ncells; i++) {
        initialize_save_data(&(grid[i]), veg_lib, 
                             out_data[i], &timer);
    }
    if (mpi_rank == VIC_MPI_ROOT) {
        // count the number of streams and variables in the global parameter file
        count_nstreams_nvars(filep->globalparam, &(options.Noutstreams),
                             nstream_vars);

        // If there weren't any output streams specified, get the defaults
        if (options.Noutstreams == 0) {
            default_outputs = true;
            get_default_nstreams_nvars(&(options.Noutstreams), nstream_vars);
        }
    }

    // broadcast Noutstreams and nstream_vars
    status = MPI_Bcast(&(options.Noutstreams), 1, MPI_AINT, VIC_MPI_ROOT,
                       MPI_COMM_VIC);
    check_mpi_status(status, "MPI Error.");
    status = MPI_Bcast(&(nstream_vars), MAX_OUTPUT_STREAMS, MPI_AINT,
                       VIC_MPI_ROOT,
                       MPI_COMM_VIC);
    check_mpi_status(status, "MPI Error.");

    // allocate output streams
    streams = calloc(options.Noutstreams, sizeof(*streams));
    check_alloc_status(streams, "Memory allocation error.");

    // allocate netcdf history files array
    nc_hist_files = calloc(options.Noutstreams, sizeof(*nc_hist_files));
    check_alloc_status(nc_hist_files, "Memory allocation error.");

    // allocate memory for streams, initialize to default/missing values
    for (streamnum = 0; streamnum < options.Noutstreams; streamnum++) {
        setup_stream(&(streams[streamnum]), nstream_vars[streamnum],
                     Ncells);
    }

    if (mpi_rank == VIC_MPI_ROOT) {
        if (default_outputs) {
            // determine which variables will be written to the history file
            set_output_defaults(&streams, dmy, NETCDF4_CLASSIC);
        }
        else {
            // set output defaults
            parse_output_info(filep->globalparam, &streams, dmy);
        }
    }

    // Now broadcast the arrays of shape nvars
    for (streamnum = 0; streamnum < options.Noutstreams; streamnum++) {
        // prefix
        status = MPI_Bcast(streams[streamnum].prefix,
                           MAXSTRING, MPI_CHAR, VIC_MPI_ROOT, MPI_COMM_VIC);
        check_mpi_status(status, "MPI error.");

        // filename
        status = MPI_Bcast(streams[streamnum].filename,
                           MAXSTRING, MPI_CHAR, VIC_MPI_ROOT, MPI_COMM_VIC);
        check_mpi_status(status, "MPI error.");

        // skip fh

        // file_format
        status = MPI_Bcast(&(streams[streamnum].file_format),
                           1, MPI_UNSIGNED_SHORT, VIC_MPI_ROOT, MPI_COMM_VIC);
        check_mpi_status(status, "MPI error.");

        // compress
        status = MPI_Bcast(&(streams[streamnum].compress),
                           1, MPI_SHORT, VIC_MPI_ROOT, MPI_COMM_VIC);
        check_mpi_status(status, "MPI error.");

        // type
        status = MPI_Bcast(streams[streamnum].type,
                           streams[streamnum].nvars,
                           MPI_UNSIGNED_SHORT, VIC_MPI_ROOT, MPI_COMM_VIC);
        check_mpi_status(status, "MPI error.");

        // mult
        status = MPI_Bcast(streams[streamnum].mult,
                           streams[streamnum].nvars,
                           MPI_DOUBLE, VIC_MPI_ROOT, MPI_COMM_VIC);
        check_mpi_status(status, "MPI error.");

        // format
        // skip broadcast

        // varid
        status = MPI_Bcast(streams[streamnum].varid,
                           streams[streamnum].nvars,
                           MPI_UNSIGNED, VIC_MPI_ROOT, MPI_COMM_VIC);
        check_mpi_status(status, "MPI error.");

        // aggtype
        status = MPI_Bcast(streams[streamnum].aggtype,
                           streams[streamnum].nvars, MPI_UNSIGNED_SHORT,
                           VIC_MPI_ROOT, MPI_COMM_VIC);
        check_mpi_status(status, "MPI error.");

        // skip agg data

        // Now brodcast the alarms
        status = MPI_Bcast(&(streams[streamnum].agg_alarm), 1,
                           mpi_alarm_struct_type, VIC_MPI_ROOT, MPI_COMM_VIC);
        check_mpi_status(status, "MPI error.");
        status = MPI_Bcast(&(streams[streamnum].write_alarm), 1,
                           mpi_alarm_struct_type, VIC_MPI_ROOT, MPI_COMM_VIC);
        check_mpi_status(status, "MPI error.");

        // allocate agg data
        alloc_aggdata(&(streams[streamnum]));

        // setup netcdf files
        initialize_nc_file(&(nc_hist_files[streamnum]),
                           streams[streamnum].nvars,
                           streams[streamnum].varid,
                           streams[streamnum].type);
    }
    // validate streams
    validate_streams(&streams);
}