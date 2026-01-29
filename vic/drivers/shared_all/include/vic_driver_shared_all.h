/******************************************************************************
 * @section DESCRIPTION
 *
 * Header file for vic_driver_shared_all routines
 *****************************************************************************/

#ifndef VIC_DRIVER_SHARED_H
#define VIC_DRIVER_SHARED_H

#include <vic_run.h>
#include <vic_version.h>

// Define maximum array sizes for driver level objects
#define MAX_FORCE_FILES 2
#define MAX_OUTPUT_STREAMS 20

// Output compression setting
#define COMPRESSION_LVL_UNSET -1
#define COMPRESSION_LVL_DEFAULT 5

// Default ouput values
#define OUT_MULT_DEFAULT 0  // Why is this not 1?
#define OUT_ASCII_FORMAT_DEFAULT "%.4f"

// Default snow band setting
#define SNOW_BAND_TRUE_BUT_UNSET 99999

// Max counter for root distribution iteration
#define MAX_ROOT_ITER 9999

/******************************************************************************
 * @brief   File formats
 *****************************************************************************/
enum
{
    UNSET_FILE_FORMAT,
    ASCII,
    BINARY,
    NETCDF3_CLASSIC,
    NETCDF3_64BIT_OFFSET,
    NETCDF4_CLASSIC,
    NETCDF4
};

/******************************************************************************
 * @brief   endian flags
 *****************************************************************************/
enum
{
    LITTLE,  /**< little-endian flag */
    BIG      /**< big-endian flag */
};

/******************************************************************************
 * @brief   Veg param sources
 *****************************************************************************/
enum
{
    FROM_DEFAULT,
    FROM_VEGLIB,
    FROM_VEGPARAM,
    FROM_VEGHIST
};

/******************************************************************************
 * @brief   Forcing Variable Types
 *****************************************************************************/
enum
{
    AIR_TEMP,    /**< air temperature per time step [C] */
    ALBEDO,      /**< surface albedo [fraction] */
    CATM,        /**< atmospheric CO2 concentration [ppm] */
    CHANNEL_IN,  /**< incoming channel flow [m3] */
    FCANOPY,     /**< fractional area covered by plant canopy [fraction] */
    FDIR,        /**< fraction of incoming shortwave that is direct [fraction] */
    LAI,         /**< leaf area index [m2/m2] */
    LWDOWN,      /**< incoming longwave radiation [W/m2] */
    PAR,         /**< incoming photosynthetically active radiation [W/m2] */
    PREC,        /**< total precipitation (rain and snow) [mm] */
    RAINF,       /**< rain [mm] */
    SNOWF,       /**< snow [mm] */
    PRESSURE,    /**< atmospheric pressure [kPa] */
    QAIR,        /**< specific humidity [kg/kg] */
    REL_HUMID,   /**< rel_humid [%] */
    SWDOWN,      /**< incoming shortwave [W/m2] */
    WIND,        /**< wind speed [m/s] */
    SKIP,        /**< place holder for unused data columns */
    // Last value of enum - DO NOT ADD ANYTHING BELOW THIS LINE!!
    // used as a loop counter and must be >= the largest value in this enum
    N_FORCING_TYPES  /**< Number of forcing types */
};


/******************************************************************************
 * @brief   Output Variable Types
 *****************************************************************************/
enum
{
    // Water Balance Terms - state variables
    OUT_ASAT,             /**< Saturated Area Fraction */
    OUT_DELDEPTH,         /**< snow depth increasing rate */
    OUT_GLAC_EXCESS,      /**< glacier snow excess flow */
    OUT_ROOTMOIST,        /**< root zone soil moisture  [mm] */
    OUT_NEW_DENSITY,
    OUT_SMFROZFRAC,       /**< fraction of soil moisture (by mass) that is ice, for each soil layer */
    OUT_SMLIQFRAC,        /**< fraction of soil moisture (by mass) that is liquid, for each soil layer */
    OUT_SNOW_AGE,
    OUT_SNOW_CANOPY,      /**< snow interception storage in canopy  [mm] */
    OUT_SNOW_COVER,       /**< fractional area of snow cover [fraction] */
    OUT_SNOW_COMB,
    OUT_SNOW_DEPTH,       /**< depth of snow pack [cm] */
    OUT_SNOW_DENSITY,
    OUT_SNOW_PACK_ICE,
    OUT_SNOW_PACK_LIQ,
    OUT_SNOW_ICEFRAC,
    OUT_SNOW_LIQFRAC,
    OUT_SNOW_POROSITY,
    OUT_SOIL_ICE,         /**< soil ice content  [mm] for each soil layer */
    OUT_SOIL_LIQ,         /**< soil liquid content  [mm] for each soil layer */
    OUT_SOIL_MOIST,       /**< soil total moisture content  [mm] for each soil layer */
    OUT_SOIL_WET,         /**< vertical average of (soil moisture - wilting point)/(maximum soil moisture - wilting point) [mm/mm] */
    OUT_SURFSTOR,         /**< storage of liquid water and ice (not snow) on surface (ponding) [mm] */
    OUT_SURF_FROST_FRAC,  /**< fraction of soil surface that is frozen [fraction] */
    OUT_SWE,              /**< snow water equivalent in snow pack (including vegetation-intercepted snow)  [mm] */
    OUT_SNOW_TRANSP,
    OUT_WDEW,             /**< total moisture interception storage in canopy [mm] */
    OUT_ZWT,              /**< water table position [cm] (zwt within lowest unsaturated layer) */
    OUT_ZWT_LUMPED,       /**< lumped water table position [cm] (zwt of total moisture across all layers, lumped together) */
    // Water Balance Terms - fluxes
    OUT_BASEFLOW,         /**< baseflow out of the bottom layer  [mm] */
    OUT_CONDEN_GRND,
    OUT_DISCHARGE,        /**< river discharge [m3 s-1]) */
    OUT_DEW_CANOP,
    OUT_DEW_SOIL,
    OUT_DELSOILMOIST,     /**< change in soil water content  [mm] */
    OUT_EVAP,             /**< total net evaporation [mm] */
    OUT_EVAP_BARE,        /**< net evaporation from bare soil [mm] */
    OUT_EVAP_CANOP,       /**< net evaporation from canopy interception [mm] */
    OUT_INFLOW,           /**< moisture that reaches top of soil column [mm] */
    OUT_FROST_CANOP,
    OUT_FROST_SNOW,
    OUT_NET_EVAP,
    OUT_PREC,             /**< incoming precipitation [mm] */
    OUT_PACK_OUTFLOW,
    OUT_RAINF,            /**< rainfall  [mm] */
    OUT_RUNOFF,           /**< surface runoff [mm] */
    OUT_SNOW_MELT,        /**< snow melt  [mm] */
    OUT_SNOWF,            /**< snowfall  [mm] */
    OUT_SUB_CANOP,        /**< net sublimation from snow stored in canopy [mm] */
    OUT_SUB_BLOWING,      /**< net sublimation of blowing snow [mm] */
    OUT_SUB_SNOW,         /**< total net sublimation from snow pack (surface and blowing) [mm] */
    OUT_TRANSP_VEG,       /**< net transpiration from vegetation [mm] */
    OUT_VAPOR_CANOP,
    OUT_VAPOR_GRND,
    OUT_WATER_ERROR,      /**< water budget error [mm] */


    OUT_DELINTERCEPT,     /**< change in canopy interception storage  [mm] */
    
    OUT_DELSURFSTOR,      /**< change in surface liquid water storage  [mm] */
    OUT_DELSWE,           /**< change in snow water equivalent  [mm] */
    OUT_REFREEZE,         /**< refreezing of water in the snow  [mm] */
    OUT_SUB_SURFACE,      /**< net sublimation from snow pack surface [mm] */

    // Energy Balance Terms - state variables
    OUT_ALBEDO,           /**< average surface albedo [fraction] */
    OUT_BARESOILT,        /**< bare soil surface temperature [C] */
    OUT_FDEPTH,           /**< depth of freezing fronts [cm] */
    OUT_RAD_TEMP,         /**< average radiative surface temperature [K] */
    OUT_SALBEDO,          /**< snow pack albedo [fraction] */
    OUT_SNOW_PACK_TEMP,   /**< snow pack temperature [C] */
    OUT_SNOW_SURF_TEMP,   /**< snow surface temperature [C] */
    OUT_SOIL_TEMP,        /**< soil temperature [C] */
    OUT_SURF_TEMP,
    OUT_TRND_FBFLAG,      /**< surface temperature flag */
    OUT_TCAN_FBFLAG,      /**< Tcanopy flag */
    OUT_TDEPTH,           /**< depth of thawing fronts [cm] */
    OUT_VEGT,             /**< average vegetation canopy temperature [C] */
    OUT_VEGTAIR,
    // Energy Balance Terms - fluxes
    OUT_ADVECTION,        /**< advected energy [W/m2] */
    OUT_ADVECTSUB,
    OUT_ADVECTGRND,
    OUT_ADVECTOVER,
    OUT_DELTACC,          /**< rate of change in cold content in snow pack [W/m2] */
    OUT_DELTAH,           /**< rate of change in heat storage [W/m2] */
    OUT_ENERGY_ERROR,     /**< energy budget error [W/m2] */
    OUT_FUSION,           /**< net energy used to melt/freeze soil moisture [W/m2] */
    OUT_GRND_FLUX,        /**< net heat flux into ground [W/m2] */
    OUT_GRND_SUB,
    OUT_GRND_GRND,
    OUT_IN_LONG,          /**< incoming longwave at ground surface (under veg) [W/m2] */
    OUT_LATENT,           /**< net upward latent heat flux [W/m2] */
    OUT_LATENT_SUB,       /**< net upward latent heat flux from sublimation [W/m2] */
    OUT_LATENT_GRND,
    OUT_LATENT_CANOP,
    OUT_LATENT_TRANSP,
    OUT_MELT_ENERGY,      /**< energy of fusion (melting) in snowpack [W/m2] */
    OUT_LWNET,            /**< net downward longwave flux [W/m2] */
    OUT_LWSUB,
    OUT_LWGRND,
    OUT_LWOVER,
    OUT_SWNET,            /**< net downward shortwave flux [W/m2] */
    OUT_SWSUB,
    OUT_SWGRND,
    OUT_SENSIBLE,         /**< net upward sensible heat flux [W/m2] */
    OUT_SENSIBLE_SUB,
    OUT_SENSIBLE_GRND,
    OUT_SENSIBLE_OVER,
    // Miscellaneous Terms
    OUT_RA_EVAP,          /**< surface aerodynamic conductance [m/s] */
    OUT_RA_GRND,          /**< overstory aerodynamic conductance [m/s] */
    OUT_RA_LEAF,          /**< "scene"canopy aerodynamic resistance [s/m]*/
    OUT_RA_OVER,          /**< surface aerodynamic resistance [s/m] */
    OUT_RA_SUB,           /**< overstory aerodynamic resistance [s/m] */
    OUT_AIR_TEMP,         /**< air temperature [C] */
    OUT_CATM,             /**< atmospheric CO2 concentrtaion [ppm]*/
    OUT_DENSITY,          /**< near-surface atmospheric density [kg/m3]*/
    OUT_FCANOPY,          /**< fractional area covered by plant canopy [fraction] */
    OUT_FDIR,             /**< fraction of incoming shortwave that is direct [fraction]*/
    OUT_COSZEN,
    OUT_LAI,              /**< leaf area index [m2/m2] */
    OUT_LWDOWN,           /**< incoming longwave [W/m2] */
    OUT_PAR,              /**< incoming photosynthetically active radiation [W/m2] */
    OUT_PRESSURE,         /**< near surface atmospheric pressure [kPa] */
    OUT_QAIR,             /**< specific humidity [kg/kg] */
    OUT_REL_HUMID,        /**< relative humidity [%] */
    OUT_SWDOWN,           /**< incoming shortwave [W/m2] */
    OUT_SURF_COND,        /**< surface conductance [m/s] */
    OUT_VP,               /**< near surface vapor pressure [kPa] */
    OUT_VPD,              /**< near surface vapor pressure deficit [kPa] */
    OUT_WIND,             /**< near surface wind speed [m/s] */
    // Band-specific quantities
    OUT_ADVECTION_BAND,   /**< advected energy [W/m2] */
    OUT_GRND_FLUX_BAND,   /**< net heat flux into ground [W/m2] */
    OUT_LATENT_BAND,      /**< net upward latent heat flux [W/m2] */
    OUT_LATENT_SUB_BAND,  /**< net upward latent heat flux due to sublimation [W/m2] */
    OUT_LWNET_BAND,       /**< net downward longwave flux [W/m2] */
    OUT_SWNET_BAND,       /**< net downward shortwave flux [W/m2] */
    OUT_SENSIBLE_BAND,    /**< net upward sensible heat flux [W/m2] */
    OUT_SNOW_CANOPY_BAND, /**< snow interception storage in canopy [mm] */
    OUT_SNOW_COVER_BAND,  /**< fractional area of snow cover [fraction] */
    OUT_SNOW_DEPTH_BAND,  /**< depth of snow pack [cm] */
    OUT_SNOW_FLUX_BAND,   /**< energy flux through snow pack [W/m2] */
    OUT_SNOW_MELT_BAND,   /**< snow melt [mm] */
    OUT_SNOW_PACKT_BAND,  /**< snow pack temperature [C] */
    OUT_SWE_BAND,         /**< snow water equivalent in snow pack [mm] */
    // Carbon-Cycling Terms
    OUT_APAR,             /**< absorbed PAR [W/m2] */
    OUT_GPP,              /**< gross primary productivity [g C/m2d] */
    OUT_RAUT,             /**< autotrophic respiration [g C/m2d] */
    OUT_NPP,              /**< net primary productivity [g C/m2d] */
    OUT_LITTERFALL,       /**< flux of carbon from living biomass into soil [g C/m2d] */
    OUT_RHET,             /**< soil respiration (heterotrophic respiration) [g C/m2d] */
    OUT_NEE,              /**< net ecosystem exchange (=NPP-RHET) [g C/m2d] */
    OUT_CLITTER,          /**< Carbon density in litter pool [g C/m2] */
    OUT_CINTER,           /**< Carbon density in intermediate pool [g C/m2] */
    OUT_CSLOW,            /**< Carbon density in slow pool [g C/m2] */
    // Timing and Profiling Terms
    OUT_TIME_VICRUN_WALL, /**< Wall time spent inside vic_run [seconds] */
    OUT_TIME_VICRUN_CPU,  /**< Wall time spent inside vic_run [seconds] */
    
    //Glacier Water Balance Terms - state variables
    OUT_GLAC_WAT_STOR,     /* glacier water storage [mm] */
    OUT_GLAC_AREA,         /* glacier surface area fraction */

    //Glacier Water Balance Terms - fluxes
    OUT_GLAC_MBAL,         /* glacier mass balance [mm] */
    OUT_GLAC_IMBAL,        /* glacier ice mass balance [mm] */
    OUT_GLAC_ACCUM,        /* glacier ice accumulation from conversion of firn to ice [mm] */
    OUT_GLAC_MELT,         /* glacier ice melt [mm] */
    OUT_GLAC_INFLOW,       /* glacier water inflow from snow melt, ice melt and rainfall [mm] */
    OUT_GLAC_OUTFLOW,      /* glacier water outflow [mm] */

    // Last value of enum - DO NOT ADD ANYTHING BELOW THIS LINE!!
    // used as a loop counter and must be >= the largest value in this enum
    N_OUTVAR_TYPES        /**< used as a loop counter*/
};

/******************************************************************************
 * @brief   Output state variable.
 *****************************************************************************/
enum
{
    STATE_SOIL_MOISTURE,               /**<  total soil moisture */
    STATE_SOIL_ICE,                    /**<  ice content */
    STATE_CANOPY_WATER,                /**<  dew storage: tmpval = veg_var[veg][band].Wdew; */
    STATE_ANNUALNPP,                   /**<  cumulative NPP: tmpval = veg_var[veg][band].AnnualNPP; */
    STATE_ANNUALNPPPREV,               /**<  previous NPP: tmpval = veg_var[veg][band].AnnualNPPPrev; */
    STATE_CLITTER,                     /**<  litter carbon: tmpval = cell[veg][band].CLitter; */
    STATE_CINTER,                      /**<  intermediate carbon: tmpval = cell[veg][band].CInter; */
    STATE_CSLOW,                       /**<  slow carbon: tmpval = cell[veg][band].CSlow; */
    STATE_SNOW_AGE,                    /**<  snow age: snow[veg][band].last_snow */
    STATE_SNOW_MELT_STATE,             /**<  melting state: (int)snow[veg][band].MELTING */
    STATE_SNOW_COVERAGE,               /**<  snow covered fraction: snow[veg][band].coverage */
    STATE_SNOW_WATER_EQUIVALENT,       /**<  snow water equivalent: snow[veg][band].swq */
    STATE_SNOW_SURF_TEMP,              /**<  snow surface temperature: snow[veg][band].surf_temp */
    STATE_SNOW_SURF_WATER,             /**<  snow surface water: snow[veg][band].surf_water */
    STATE_SNOW_PACK_TEMP,              /**<  snow pack temperature: snow[veg][band].pack_temp */
    STATE_SNOW_PACK_WATER,             /**<  snow pack water: snow[veg][band].pack_water */
    STATE_SNOW_DENSITY,                /**<  snow density: snow[veg][band].density */
    STATE_SNOW_COLD_CONTENT,           /**<  snow cold content: snow[veg][band].coldcontent */
    STATE_SNOW_CANOPY,                 /**<  snow canopy storage: snow[veg][band].snow_canopy */
    STATE_SOIL_NODE_TEMP,              /**<  soil node temperatures: energy[veg][band].T[nidx] */
    STATE_FOLIAGE_TEMPERATURE,         /**<  Foliage temperature: energy[veg][band].Tfoliage */
    STATE_ENERGY_LONGUNDEROUT,         /**<  Outgoing longwave from understory: energy[veg][band].LongUnderOut */
    STATE_ENERGY_SNOW_FLUX,            /**<  Thermal flux through the snow pack: energy[veg][band].snow_flux */
    STATE_AVG_ALBEDO,                  /**<  gridcell-averaged albedo: gridcell_avg.avg_albedo */
    // Last value of enum - DO NOT ADD ANYTHING BELOW THIS LINE!!
    // used as a loop counter and must be >= the largest value in this enum
    N_STATE_VARS                       /**< used as a loop counter*/
};


/******************************************************************************
 * @brief   Output BINARY format types
 *****************************************************************************/
enum
{
    OUT_TYPE_DEFAULT,  /**< Default data type */
    OUT_TYPE_CHAR,     /**< char */
    OUT_TYPE_SINT,     /**< short int */
    OUT_TYPE_USINT,    /**< unsigned short int */
    OUT_TYPE_INT,      /**< int */
    OUT_TYPE_FLOAT,    /**< single-precision floating point */
    OUT_TYPE_DOUBLE    /**< double-precision floating point */
};

/******************************************************************************
 * @brief   Output aggregation method types
 *****************************************************************************/
enum
{
    AGG_TYPE_DEFAULT, /**< Default aggregation type */
    AGG_TYPE_AVG,     /**< average over agg interval */
    AGG_TYPE_BEG,     /**< value at beginning of agg interval */
    AGG_TYPE_END,     /**< value at end of agg interval */
    AGG_TYPE_MAX,     /**< maximum value over agg interval */
    AGG_TYPE_MIN,     /**< minimum value over agg interval */
    AGG_TYPE_SUM      /**< sum over agg interval */
};

/******************************************************************************
 * @brief   Frequency flags for raising alarms/flags
 *****************************************************************************/
enum
{
    FREQ_NEVER,      /**< Flag for never raising alarm */
    FREQ_NSTEPS,     /**< Flag for raising alarm every nsteps */
    FREQ_NSECONDS,   /**< Flag for raising alarm every nseconds */
    FREQ_NMINUTES,   /**< Flag for raising alarm every nminutes */
    FREQ_NHOURS,     /**< Flag for raising alarm every nhours */
    FREQ_NDAYS,      /**< Flag for raising alarm every ndays */
    FREQ_NMONTHS,    /**< Flag for raising alarm every nmonths */
    FREQ_NYEARS,     /**< Flag for raising alarm every nyears */
    FREQ_DATE,       /**< Flag for raising alarm on a specific date */
    FREQ_END         /**< Flag for raising alarm at the end of a simulation */
};

/******************************************************************************
 * @brief   Codes for displaying version information
 *****************************************************************************/
enum
{
    DISP_VERSION,
    DISP_COMPILE_TIME,
    DISP_ALL
};

/******************************************************************************
 * @brief   Codes for calendar option.
 *****************************************************************************/
enum calendars
{
    CALENDAR_STANDARD,
    CALENDAR_GREGORIAN,
    CALENDAR_PROLEPTIC_GREGORIAN,
    CALENDAR_NOLEAP,
    CALENDAR_365_DAY,
    CALENDAR_360_DAY,
    CALENDAR_JULIAN,
    CALENDAR_ALL_LEAP,
    CALENDAR_366_DAY
};

/******************************************************************************
 * @brief   Codes for time units option.
 *****************************************************************************/
enum time_units
{
    TIME_UNITS_SECONDS,
    TIME_UNITS_MINUTES,
    TIME_UNITS_HOURS,
    TIME_UNITS_DAYS
};

/******************************************************************************
 * @brief   Codes for timers
 *****************************************************************************/
enum timers
{
    TIMER_VIC_ALL,
    TIMER_VIC_INIT,
    TIMER_VIC_RUN,
    TIMER_VIC_FINAL,
    TIMER_VIC_FORCE,
    TIMER_VIC_WRITE,
    N_TIMERS
};

/******************************************************************************
 * @brief    Stores forcing file input information.
 *****************************************************************************/
typedef struct {
    size_t N_ELEM; /**< number of elements per record; for LAI and ALBEDO,
                        1 element per veg tile; for others N_ELEM = 1; */
    bool SIGNED;
    bool SUPPLIED;
    double multiplier;
    char varname[MAXSTRING];
} force_type_struct;

/******************************************************************************
 * @brief    This structure records the parameters set by the forcing file
             input routines.  Those filled, are used to estimate the paramters
             needed for the model run in initialize_atmos.c.
 *****************************************************************************/
typedef struct {
    force_type_struct TYPE[N_FORCING_TYPES];
    double FORCE_DT[2];    /**< forcing file time step */
    size_t force_steps_per_day[2];    /**< forcing file timesteps per day */
    unsigned short int FORCE_ENDIAN[2];  /**< endian-ness of input file, used for
                                            DAILY_BINARY format */
    int FORCE_FORMAT[2];            /**< ASCII or BINARY */
    int FORCE_INDEX[2][N_FORCING_TYPES];
    size_t N_TYPES[2];
} param_set_struct;

/******************************************************************************
 * @brief   This structure stores alarm information
 *****************************************************************************/
typedef struct {
    unsigned int count;  /**< current alarm count */
    dmy_struct next_dmy; /**< next dmy to raise alarm at */
    int next_count;      /**< next count to raise alarm at */
    unsigned int freq;   /**< enum value to describing alarm frequency */
    int n;               /**< variable that provides additional information with respect to alarm_freq */
    bool is_subdaily;    /**< flag denoting if alarm will be raised more than once per day */
} alarm_struct;

/******************************************************************************
 * @brief   This structure stores output information for one output stream.
 *****************************************************************************/
typedef struct {
    size_t nvars;                    /**< number of variables to store in the file */
    size_t ngridcells;               /**< number of grid cells in aggdata */
    dmy_struct time_bounds[2];       /**< timestep bounds of stream */
    char prefix[MAXSTRING];          /**< prefix of the file name, e.g. "fluxes" */
    char filename[MAXSTRING];        /**< complete file name */
    FILE *fh;                        /**< filehandle */
    unsigned short int file_format;  /**< output file format */
    short int compress;              /**< Compress output files in stream*/
    unsigned short int *type;        /**< type, when written to a binary file;
                                          OUT_TYPE_USINT  = unsigned short int
                                          OUT_TYPE_SINT   = short int
                                          OUT_TYPE_FLOAT  = single precision floating point
                                          OUT_TYPE_DOUBLE = double precision floating point */
    double *mult;                    /**< multiplier, when written to a binary file [shape=(nvars, )] */
    char **format;                    /**< format, when written to disk [shape=(nvars, )] */
    unsigned int *varid;             /**< id numbers of the variables to store in the file
                                          (a variable's id number is its index in the out_data array).
                                          The order of the id numbers in the varid array
                                          is the order in which the variables will be written. */
    unsigned short int *aggtype;     /**< type of aggregation to use [shape=(nvars, )] */
    double ****aggdata;              /**< array of aggregated data values [shape=(ngridcells, nvars, nelem, nbins)] */
    alarm_struct agg_alarm;          /**< alaram for stream aggregation */
    alarm_struct write_alarm;        /**< alaram for controlling stream write */
} stream_struct;

/******************************************************************************
 * @brief   This structure stores moisture state information for differencing
 *          with next time step.
 *****************************************************************************/
typedef struct {
    double total_moist_storage;   /**< total moisture storage [mm] */
    double total_soil_moist;      /**< total column soil moisture [mm] */
    double surfstor;              /**< surface water storage [mm] */
    double swe;                   /**< snow water equivalent [mm] */
    double wdew;                  /**< canopy interception [mm] */
} save_data_struct;

/******************************************************************************
 * @brief   This structure stores metadata for individual variables
 *****************************************************************************/
typedef struct {
    char varname[MAXSTRING];  /**< name of variable */
    char long_name[MAXSTRING];  /**< name of variable */
    char standard_name[MAXSTRING];  /**< cf long_name of variable */
    char units[MAXSTRING];  /**< units of variable */
    char description[MAXSTRING];  /**< descripition of variable */
    size_t nelem;          /**< number of data values */
} metadata_struct;

/******************************************************************************
 * @brief   This structure holds all variables needed for the error handling
 *          routines.
 *****************************************************************************/
typedef struct {
    force_data_struct *force;
    double dt;
    energy_bal_struct *energy;
    size_t rec;
    double **out_data;
    stream_struct *output_streams;
    snow_data_struct *snow;
    soil_con_struct soil_con;
    veg_con_struct *veg_con;
    veg_var_struct *veg_var;
} Error_struct;

/******************************************************************************
 * @brief   This structure holds timer information for profiling
 *****************************************************************************/
typedef struct {
    double start_wall;
    double start_cpu;
    double stop_wall;
    double stop_cpu;
    double delta_wall;
    double delta_cpu;
} timer_struct;

double air_density(double t, double p, double vp);
void agg_stream_data(stream_struct *stream, dmy_struct *dmy_current,
                     double ***out_data);
double all_30_day_from_dmy(dmy_struct *dmy);
double all_leap_from_dmy(dmy_struct *dmy);
void alloc_aggdata(stream_struct *stream);
void alloc_out_data(size_t ngridcells, double ***out_data);
double average(double *ar, size_t n);
double calc_energy_balance_error(double, double, double, double, double);
double calc_water_balance_error(double, double, double, double);
bool cell_method_from_agg_type(unsigned short int aggtype, char cell_method[]);
bool check_write_flag(int rec);
void collect_eb_terms(energy_bal_struct, snow_data_struct,
                      cell_data_struct, double, bool, double, 
                      int, double **);
void collect_wb_terms(cell_data_struct, veg_var_struct,
                      snow_data_struct, double, bool, double, double **);
void compute_derived_state_vars(all_vars_struct *, soil_con_struct *,
                                veg_con_struct *);
size_t count_force_vars(FILE *gp);
void count_nstreams_nvars(FILE *gp, size_t *nstreams, size_t nvars[]);
void cmd_proc(int argc, char **argv, char *globalfilename);
void compress_files(char string[], short int level);
stream_struct create_outstream(stream_struct *output_streams);
double get_cpu_time();
void get_current_datetime(char *cdt);
double get_wall_time();
double date2num(double origin, dmy_struct *date, double tzoffset,
                unsigned short int calendar, unsigned short int time_units);
void dmy_all_30_day(double julian, dmy_struct *dmy);
void dmy_all_leap(double julian, dmy_struct *dmy);
bool dmy_equal(dmy_struct *a, dmy_struct *b);
void dmy_julian_day(double julian, unsigned short int calendar,
                    dmy_struct *dmy);
void dmy_no_leap_day(double julian, dmy_struct *dmy);
void dt_seconds_to_time_units(unsigned short int time_units, double dt_seconds,
                              double *dt_time_units);
void display_current_settings(int);
double fractional_day_from_dmy(dmy_struct *dmy);
void free_all_vars(all_vars_struct *all_vars, int Nveg);
void free_dmy(dmy_struct **dmy);
void free_out_data(size_t ngridcells, double ***out_data);
void free_streams(stream_struct **streams);
void free_vegcon(veg_con_struct **veg_con);
void generate_default_state(all_vars_struct *, soil_con_struct *,
                            veg_con_struct *);
void get_default_nstreams_nvars(size_t *nstreams, size_t nvars[]);
void get_parameters(FILE *paramfile);
void init_output_list(double **out_data, int write, char *format, int type,
                      double mult);
void initialize_energy(energy_bal_struct *energy, size_t nveg);
void initialize_global(void);
void initialize_options(void);
void initialize_parameters(void);
void initialize_save_data(all_vars_struct *all_vars, force_data_struct *force,
                          soil_con_struct *soil_con, veg_con_struct *veg_con,
                          veg_lib_struct *veg_lib,
                          double **out_data, save_data_struct *save_data,
                          timer_struct *timer);
void initialize_snow(snow_data_struct *snow, size_t veg_num);
void initialize_soil(cell_data_struct *cell, size_t veg_num);
void initialize_time(void);
void initialize_veg(veg_var_struct *veg_var, size_t nveg);
double julian_day_from_dmy(dmy_struct *dmy, unsigned short int calendar);
bool leap_year(unsigned short int year, unsigned short int calendar);
all_vars_struct make_all_vars(size_t nveg);
cell_data_struct *make_cell_data(size_t veg_type_num);
dmy_struct *make_dmy(global_param_struct *global);
energy_bal_struct *make_energy_bal(size_t nveg);
void make_lastday(unsigned short int calendar, unsigned short int year,
                  unsigned short int lastday[]);
snow_data_struct *make_snow_data(size_t nveg);
veg_var_struct *make_veg_var(size_t veg_type_num);
double no_leap_day_from_dmy(dmy_struct *dmy);
void num2date(double origin, double time_value, double tzoffset,
              unsigned short int calendar, unsigned short int time_units,
              dmy_struct *date);
FILE *open_file(char string[], char type[]);
void parse_nc_time_units(char *nc_unit_chars, unsigned short int *units,
                         dmy_struct *dmy);
void put_data(all_vars_struct *, force_data_struct *, soil_con_struct *,
              veg_con_struct *, veg_lib_struct *veg_lib,
              double **out_data, save_data_struct *, timer_struct *timer);
void print_alarm(alarm_struct *alarm);
void print_cell_data(cell_data_struct *cell, size_t nlayers, size_t nfrost);
void print_dmy(dmy_struct *dmy);
void print_energy_bal(energy_bal_struct *eb, size_t nnodes, size_t nfronts);
void print_force_type(force_type_struct *force_type);
void print_global_param(global_param_struct *gp);
void print_license(void);
void print_option(option_struct *option);
void print_out_data(double **out_data, metadata_struct *metadata);
void print_out_metadata(metadata_struct *metadata, size_t nvars);
void print_output_streams(stream_struct *outf);
void print_param_set(param_set_struct *param_set);
void print_parameters(parameters_struct *param);
void print_save_data(save_data_struct *save);
void print_snow_data(snow_data_struct *snow);
void print_soil_con(soil_con_struct *scon, size_t nlayers, size_t nnodes,
                    size_t nfrost, size_t nbands, size_t nzwt);
void print_stream(stream_struct *stream, metadata_struct *metadata);
void print_veg_con(veg_con_struct *vcon, size_t nroots,
                   char carbon, size_t ncanopy);
void print_veg_lib(veg_lib_struct *vlib, char carbon);
void print_veg_var(veg_var_struct *vvar, size_t ncanopy, size_t nswband);
void print_version(char *);
void print_usage(char *);
double q_to_vp(double q, double p);
bool raise_alarm(alarm_struct *alarm, dmy_struct *dmy_current);
void reset_alarm(alarm_struct *alarm, dmy_struct *dmy_current);
void reset_stream(stream_struct *stream, dmy_struct *dmy_current);
void set_output_var(stream_struct *stream, char *varname, size_t varnum,
                    char *format, unsigned short int type, double mult,
                    unsigned short int aggtype);
unsigned int get_default_outvar_aggtype(unsigned int varid);
void set_alarm(dmy_struct *dmy_current, unsigned int freq, void *value,
               alarm_struct *alarm);
void set_output_defaults(stream_struct **output_streams,
                         dmy_struct     *dmy_current,
                         unsigned short  default_file_format);
void set_output_met_data_info();
void setup_stream(stream_struct *stream, size_t nvars, size_t ngridcells);
void soil_moisture_from_water_table(soil_con_struct *soil_con, size_t nlayers);
void sprint_dmy(char *str, dmy_struct *dmy);
void str_from_calendar(unsigned short int calendar, char *calendar_str);
void str_from_time_units(unsigned short int time_units, char *unit_str);
unsigned short int str_to_agg_type(char aggstr[]);
void str_to_ascii_format(char *format);
bool str_to_bool(char str[]);
unsigned short int str_to_calendar(char *cal_chars);
unsigned short int str_to_freq_flag(char freq[]);
double str_to_out_mult(char multstr[]);
unsigned short int str_to_out_type(char typestr[]);
unsigned short int str_to_timeunits(char units_chars[]);
void strpdmy(const char *s, const char *format, dmy_struct *dmy);
double time_delta(dmy_struct *dmy_current, unsigned short int freq, int n);
void timer_continue(timer_struct *t);
void timer_init(timer_struct *t);
void timer_start(timer_struct *t);
void timer_stop(timer_struct *t);
int update_step_vars(all_vars_struct *, veg_con_struct *, veg_hist_struct *);
int invalid_date(unsigned short int calendar, dmy_struct *dmy);
void validate_parameters(void);
void validate_streams(stream_struct **stream);
void zero_output_list(double **);

#endif
