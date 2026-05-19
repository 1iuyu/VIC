/******************************************************************************
 * @section DESCRIPTION
 *
 * Classic driver of the VIC model
 *****************************************************************************/

#include <vic_driver_classic.h>
#include <rout.h>   // Routing routine (extension)

// global variables
int                 flag;
size_t              NR; /* array index for atmos struct that indicates
                           the model step avarage or sum */
size_t              NF; /* array index loop counter limit for atmos
                           struct that indicates the SNOW_STEP values */
FILE *LOG_DEST = NULL;
char vic_run_ref_str[MAXSTRING] = "";
global_param_struct global_param;
veg_lib_struct     *veg_lib;
option_struct       options;
Error_struct        Error;
param_set_struct    param_set;
parameters_struct   param;
filenames_struct    filenames;
filep_struct        filep;
metadata_struct     out_metadata[N_OUTVAR_TYPES];

/******************************************************************************
 * @brief   Classic driver of the VIC model
 * @details The classic driver runs VIC for a single grid cell for all
 *          timesteps before moving on to the next grid cell.
 *
 * @param argc Argument count
 * @param argv Argument vector
 *****************************************************************************/
int
main(int   argc,
     char *argv[])
{
    /** Variable Declarations **/
    size_t             rec;
    size_t             Nveg_type;
    int                startrec;
    int                ErrorFlag;
    int                i, n;
    int                Ncells;
    size_t             streamnum;
    dmy_struct        *dmy;
    stream_struct     *streams = NULL;         // [nstreams]
    double          ***out_data = NULL;        // [ncells, nvars, nelem]
    save_data_struct  *save_data;              // [ncells]
    timer_struct       global_timers[N_TIMERS];
    timer_struct       cell_timer;
    // Extensions
    rout_struct        rout; // Routing routine (extension)

    // start vic all timer
    timer_start(&(global_timers[TIMER_VIC_ALL]));
    // start vic init timer
    timer_start(&(global_timers[TIMER_VIC_INIT]));

    // Initialize Log Destination
    initialize_log();

    /** Read Model Options **/
    cmd_proc(argc, argv, filenames.global);

    // Initialize global structures
    initialize_options();
    initialize_global();
    initialize_parameters();
    initialize_filenames();

    /* Initilize forcing file param structure */
    initialize_forcing_files();

    /** Read Global Control File **/
    filep.globalparam = open_file(filenames.global, "r");
    get_global_param(filep.globalparam);
    fclose(filep.globalparam);

    // Set Log Destination
    setup_logging(MISSING, filenames.log_path, &(filep.logfile));

    /** Set model constants **/
    if (strcmp(filenames.constants, "MISSING") != 0) {
        filep.constants = open_file(filenames.constants, "r");
        get_parameters(filep.constants);
    }
    // Check that model parameters are valid
    validate_parameters();

    /** Make Date Data Structure **/
    initialize_time();
    dmy = make_dmy(&global_param);

    filep.globalparam = open_file(filenames.global, "r");
    parse_output_info(filep.globalparam, &streams, &(dmy[0]));
    validate_streams(&streams);

    /** Check and Open Files **/
    check_files(&filep, &filenames);

    /** Read Vegetation Library File **/
    veg_lib = read_veglib(filep.veglib, &Nveg_type);

    /** Initial state **/
    startrec = 0;
    if (options.INIT_STATE) {
        filep.init_state = check_state_file(filenames.init_state,
                                            options.Nlayer, options.Nsoil,
                                            &startrec);
    }

    /** open state file if model state is to be saved **/
    if (options.SAVE_STATE && strcmp(filenames.statefile, "NONE") != 0) {
        filep.statefile = open_state_file(&global_param, filenames,
                                          options.Nlayer,
                                          options.Nsoil);
    }
    else {
        filep.statefile = NULL;
    }

    // stop init timer
    timer_stop(&(global_timers[TIMER_VIC_INIT]));
    // start vic run timer
    timer_start(&(global_timers[TIMER_VIC_RUN]));

    /*****************************************
     * Read soil for all "active" grid cells *
     *****************************************/
    grid_cell = read_grid_data(filep.soilparam, &Ncells);

    // allocate memory for routing
    rout_alloc();   // Routing routine (extension)

    /********************************************
     Initialize grid cell model data structures.
    *********************************************/
    for (i = 0; i < Ncells; i++) {

        vic_init(Nveg_type, &grid_cell[i],
                &filep, dmy);

    }   /* End Grid Loop */

    // initialize routing parameters from parameter files
    rout_init();    // Routing routine (extension)

    // initialize output structures
    vic_init_output(Ncells, grid_cell, &filep, &(dmy[0]), veg_lib);

    /******************************************
     Run Model in Grid Cell for all Time Steps
    ******************************************/
    for (rec = startrec; rec < global_param.nrecs; rec++) {

        timer_start(&cell_timer);

        for (i = 0; i < Ncells; i++) {

            /**************************************************
               Update data structures for current time step
            **************************************************/
            ErrorFlag = update_step_vars(grid_cell[i].all_vars, 
                                         grid_cell[i].veg_con,
                                         grid_cell[i].veg_hist[rec]);

            /**************************************************
               Compute cell physics for 1 timestep
            **************************************************/
            ErrorFlag = vic_run(&grid_cell[i].force[rec], 
                                grid_cell[i].all_vars,
                                &global_param,
                                grid_cell[i].soil_con, 
                                grid_cell[i].veg_con, veg_lib);

            if (ErrorFlag == ERROR) {
                if (options.CONTINUEONERROR) {
                    // Handle grid cell solution error
                    log_warn("ERROR: Grid cell %i failed in record %zu "
                                "so the simulation has not finished.  An "
                                "incomplete output file has been "
                                "generated, check your inputs before "
                                "rerunning the simulation.",
                                grid_cell[i].gridcel, rec);
                    break;
                }
                else {
                    // Else exit program on cell solution error as in previous versions
                    log_err("ERROR: Grid cell %i failed in record %zu "
                            "so the simulation has ended. Check your "
                            "inputs before rerunning the simulation.",
                            grid_cell[i].gridcel, rec);
                }
            }
            
        }
        timer_stop(&cell_timer);
        vic_write_output(&(dmy[0]));
        /************************************
            Save model state at assigned date
            (after the final time step of the assigned date)
        ************************************/
        if (filep.statefile != NULL &&
            check_save_state_flag(dmy, rec)) {
            for (i = 0; i < Ncells; i++) {
                write_model_state(grid_cell[i].all_vars, 
                                  grid_cell[i].veg_con,
                                  grid_cell[i].gridcel, &filep, 
                                  grid_cell[i].soil_con);
            }
        }
    } /* End Rec Loop */

    close_files(&filep, &streams);
    for (i = 0; i < Ncells; i++) {

        free_veg_hist(global_param.nrecs, 
                      grid_cell[i].veg_con[0].vegetat_type_num,
                      grid_cell[i].veg_hist);
        free_all_vars(grid_cell[i].all_vars, 
                      grid_cell[i].veg_con[0].vegetat_type_num);

        free_vegcon(grid_cell[i].veg_con);
        free(grid_cell[i].soil_con->AreaFract);
        free(grid_cell[i].soil_con->BandElev);
        free(grid_cell[i].soil_con->Tfactor);
        free(grid_cell[i].soil_con->Pfactor);
    }


    // stop vic run timer
    timer_stop(&(global_timers[TIMER_VIC_RUN]));
    // start vic final timer
    timer_start(&(global_timers[TIMER_VIC_FINAL]));

    /** cleanup **/
    free_dmy(&dmy);
    free_streams(&streams);
    free_out_data(Ncells, out_data);
    fclose(filep.soilparam);
    free_veglib(&veg_lib);
    fclose(filep.vegparam);
    fclose(filep.veglib);
    if (options.SNOW_BAND > 1) {
        fclose(filep.snowband);
    }
    if (options.INIT_STATE) {
        fclose(filep.init_state);
    }
    if (options.SAVE_STATE && strcmp(filenames.statefile, "NONE") != 0) {
        fclose(filep.statefile);
    }
    finalize_logging();
    log_info("Completed running VIC %s", VIC_DRIVER);

    // stop vic final timer
    timer_stop(&(global_timers[TIMER_VIC_FINAL]));
    // stop vic all timer
    timer_stop(&(global_timers[TIMER_VIC_ALL]));
    // write timing info
    write_vic_timing_table(global_timers);

    return EXIT_SUCCESS;
}   /* End Main Program */
