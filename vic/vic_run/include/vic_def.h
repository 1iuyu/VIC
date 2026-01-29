/******************************************************************************
 * @section DESCRIPTION
 *
 * Definition header file
 *****************************************************************************/

#ifndef VIC_DEF_H
#define VIC_DEF_H

#define _BSD_SOURCE
#define __USE_XOPEN
#define _GNU_SOURCE

#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdbool.h>
#include <stddef.h>
#include <unistd.h>
#include <time.h>
#include <pwd.h>
#include <sys/time.h>
#include <sys/types.h>

#include <vic_physical_constants.h>
#include <vic_log.h>

/***** Model Constants *****/
#define MAXSTRING    2048
#define MISSING      -99999.   /**< missing value */
#define MISSING_USI  99999.    /**< missing value for unsigned ints */
#define MISSING_S    "MISSING"    /**< missing value for strings */
#define NODATA_VH    -1        /**< missing value for veg_hist inputs */
#define NODATA_VEG   -1        /**< flag for veg types not in grid cell */
#define ERROR        -999      /**< Error Flag returned by subroutines */

/***** Define maximum array sizes for model source code *****/
#define MAX_LAYERS      6      /**< maximum number of soil moisture layers */
#define MAX_NODES       17     /**< maximum number of soil thermal nodes */
#define MAX_FRONTS      3      /**< maximum number of freezing and thawing front depths to store */
#define MAX_SNOWS       3      /**< maximum number of snowpack layers */
#define MAX_SWBANDS     2      /**< maximum number of solar radiation wave bands */
#define MAX_ZWTVMOIST   11     /**< maximum number of points in water table vs moisture curve for each soil layer; should include points at lower and upper boundaries of the layer */

/***** Define minimum values for model parameters *****/
#define MINSOILDEPTH    0.001  /**< Minimum layer depth with which model can work (m) */
#define MIN_FCANOPY    0.0001  /**< Minimum allowable canopy fraction */
#define MIN_SNOW_WETFRAC 0.01  /**< Minimum fraction of snow depth to be considered wet */

/***** Define minimum and maximum values for model timesteps *****/
#define MIN_SUBDAILY_STEPS_PER_DAY  4
#define MAX_SUBDAILY_STEPS_PER_DAY  1440

#ifndef WET
#define WET 0
#define DRY 1
#endif

#ifndef SNOW
#define RAIN 0
#define SNOW 1
#endif

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

// 通用的any宏
#define any(arr, size, condition) \
    ({ \
        bool result = false; \
        for (size_t i = 0; i < (size); i++) { \
            if ((arr)[i] condition) { \
                result = true; \
                break; \
            } \
        } \
        result; \
    })

extern size_t NR;       /**< array index for force struct that indicates
                             the model step avarage or sum */
extern size_t NF;       /**< array index loop counter limit for force
                             struct that indicates the SNOW_STEP values */
extern char          vic_run_ref_str[MAXSTRING];

/******************************************************************************
 * @brief   Snow Density parametrizations
 *****************************************************************************/
enum
{
    DENS_BRAS,
    DENS_SNTHRM
};

/******************************************************************************
 * @brief   Soil water transpiration factor options
 *****************************************************************************/
enum
{
    NOAH,
    CLM
};

/******************************************************************************
 * @brief   Soil water transpiration factor options
 *****************************************************************************/
enum
{
    BATS,
    SNICAR
};
/******************************************************************************
 * @brief   Baseflow parametrizations
 *****************************************************************************/
enum
{
    ARNO,
    NIJSSEN2001
};

/******************************************************************************
 * @brief   Canopy resistance parametrizations
 *****************************************************************************/
enum
{
    RC_JARVIS,
    RC_PHOTO
};

/******************************************************************************
 * @brief   Photosynthesis parametrizations
 *****************************************************************************/
enum
{
    PS_FARQUHAR,
    PS_MONTEITH
};

/******************************************************************************
 * @brief   Photosynthetic pathways
 *****************************************************************************/
enum
{
    PHOTO_C3,
    PHOTO_C4
};

/***** Data Structures *****/

/******************************************************************************
 * @brief   This structure stores model options.
 *****************************************************************************/
typedef struct {
    // simulation modes
    bool BLOWING;        /**< TRUE = calculate sublimation from blowing snow */
    bool BLOWING_VAR_THRESHOLD;
    bool BLOWING_CALC_PROB;
    bool BLOWING_SIMPLE;
    bool BLOWING_FETCH;
    bool BLOWING_SPATIAL_WIND;
    bool CARBON;         /**< TRUE = simulate carbon cycling processes;
                            FALSE = no carbon cycling (default) */
    bool CONTINUEONERROR; /**< TRUE = VIC will continue to run after a cell has an error */
    bool CORRPREC;       /**< TRUE = correct precipitation for gage undercatch */
    bool EQUAL_AREA;     /**< TRUE = RESOLUTION stores grid cell area in km^2;
                            FALSE = RESOLUTION stores grid cell side length in degrees */
    bool EXP_TRANS;
    bool FROZEN_SOIL;    /**< TRUE = Use frozen soils code */

    size_t Ncanopy;      /**< Number of canopy layers in the model. */
    size_t Nlayer;       /**< Number of layers in model */
    size_t Nnode;        /**< Number of soil thermal nodes in the model */
    size_t Nswband;
    size_t Nroot;
    bool NOFLUX;         /**< TRUE = Use no flux lower bondary when computing
                            soil thermal fluxes */
    size_t NVEGTYPES;    /**< number of vegetation types in veg_param file */
    unsigned short int RC_MODE;        /**< RC_JARVIS = compute canopy resistance via Jarvis formulation (default)
                                          RC_PHOTO = compute canopy resistance based on photosynthetic activity */
    unsigned short int SNOW_DENSITY;   /**< DENS_BRAS: Use algorithm of Bras, 1990; DENS_SNTHRM: Use algorithm of SNTHRM89 adapted for 1-layer pack */
    size_t SNOW_BAND;    /**< Number of elevation bands over which to solve the
                            snow model */
    bool TFALLBACK;      /**< TRUE = when any temperature iterations fail to converge,
                                   use temperature from previous time step; the number
                                   of instances when this occurs will be logged and
                                   reported at the end of the cell's simulation
                            FALSE = when iterations fail to converge, report an error
                                    and abort simulation for current grid cell
                            Default = TRUE */
    // input options
    unsigned short int BASEFLOW;     /**< ARNO: read Ds, Dm, Ws, c; NIJSSEN2001: read d1, d2, d3, d4 */
    unsigned short int SOIL_TRANSP;  /**< Noah; CLM */
    unsigned short int GRID_DECIMAL; /**< Number of decimal places in grid file extensions */
    bool VEGLIB_FCAN;    /**< TRUE = veg library file contains monthly fcanopy values */
    bool VEGLIB_PHOTO;   /**< TRUE = veg library contains photosynthesis parameters */
    bool VEGPARAM_ALB;   /**< TRUE = veg param file contains monthly albedo values */
    bool VEGPARAM_FCAN;  /**< TRUE = veg param file contains monthly fcanopy values */
    bool VEGPARAM_LAI;   /**< TRUE = veg param file contains monthly LAI values */
    unsigned short int ALB_SRC;        /**< FROM_VEGLIB = use albedo values from veg library file
                                          FROM_VEGPARAM = use albedo values from the veg param file */
    unsigned short int FCAN_SRC;       /**< FROM_VEGLIB = use fcanopy values from veg library file
                                          FROM_VEGPARAM = use fcanopy values from the veg param file */
    unsigned short int LAI_SRC;        /**< FROM_VEGLIB = use LAI values from veg library file
                                          FROM_VEGPARAM = use LAI values from the veg param file */
    bool ORGANIC_FRACT;  /**< TRUE = organic matter fraction of each layer is read from the soil parameter file; otherwise set to 0.0. */
    bool BULK_DENSITY_COMB; /**< TRUE = soil bulk density (combined mineral and organic matter) read from soil parameter file; otherwise set to 0.0 */
    bool MAX_SNOW_ALBEDO; /**< TRUE = maximum snow albedo taken from the parameter file if veg type is not bare soil; otherwise use param option */

    // state options
    unsigned short int STATE_FORMAT;  /**< TRUE = model state file is binary (default) */
    unsigned short int SNOW_AGING;    /**< BATS and SNICAR */
    bool INIT_STATE;     /**< TRUE = initialize model state from file */
    bool SAVE_STATE;     /**< TRUE = save state file */
    bool STATENAME_CESM; /**< TRUE = use CESM statefile naming conventions */
    // glacier options
    int  GLACIER_ID;        /* An index indicating which veg class in the vegetation library contains glacier information. */

    // output options
    size_t Noutstreams;  /**< Number of output stream */
    bool   PRT_HEADER;     /* TRUE = insert header at beginning of output file; FALSE = no header */
    bool   ROUT_PARAM;
} option_struct;

/******************************************************************************
 * @brief   This structure stores all model run global parameters.
 *****************************************************************************/
typedef struct {
    double wind_h;                 /**< height of wind measurements (m) */
    double resolution;             /**< Model resolution (degrees) */
    double dt;                     /**< Time step in seconds */
    double snow_dt;                /**< Snow model time step in seconds */
    double runoff_dt;              /**< Runoff time step in seconds */
    double atmos_dt;               /**< Atmos time step in seconds */
    size_t model_steps_per_day;    /**< Number of model timesteps per day */
    size_t snow_steps_per_day;     /**< Number of snow timesteps per day */
    size_t runoff_steps_per_day;   /**< Number of runoff timesteps per day */
    size_t atmos_steps_per_day;    /**< Number of atmos timesteps per day */
    unsigned short int endday;     /**< Last day of model simulation */
    unsigned short int endmonth;   /**< Last month of model simulation */
    unsigned short int endyear;    /**< Last year of model simulation */
    unsigned short int forceday[2];  /**< day forcing files starts */
    unsigned int forcesec[2];          /**< seconds since midnight when forcing
                                          files starts */
    unsigned short int forcemonth[2];  /**< month forcing files starts */
    unsigned short int forceoffset[2];  /**< counter to keep track of offset in reading
                                           forcing files; updated after every read */
    unsigned int forceskip[2];   /**< number of model time steps to skip at
                                      the start of the forcing file */
    unsigned short int forceyear[2];  /**< year forcing files start */
    size_t nrecs;                /**< Number of time steps simulated */
    unsigned short int startday;  /**< Starting day of the simulation */
    unsigned short int startmonth;  /**< Starting month of the simulation */
    unsigned int startsec;          /**< Seconds since midnight when simulation
                                       will start */
    unsigned short int startyear;  /**< Starting year of the simulation */
    unsigned short int stateday;   /**< Day of the simulation at which to save
                                      model state */
    unsigned short int statemonth;  /**< Month of the simulation at which to save
                                       model state */
    unsigned int statesec;          /**< Seconds since midnight at which to save state */
    unsigned short int stateyear;  /**< Year of the simulation at which to save
                                      model state */
    unsigned short int calendar;  /**< Date/time calendar */
    unsigned short int time_units;  /**< Units for numeric times */
    double time_origin_num;        /**< Numeric date origin */
    char time_origin_str[MAXSTRING];  /**< string date origin */
} global_param_struct;

/******************************************************************************
 * @brief    This structure holds the model parameters.
 *****************************************************************************/
typedef struct {
    // Lapse Rate
    double LAPSE_RATE;  /**< temperature lapse rate (C/m) */

    // Precipitation Guage Height
    double GAUGE_HEIGHT;   /**< precipitation gauge height (m) */

    // REF_HEIGHT
    double REF_HEIGHT;
    double REF_HEIGHT_WIND;

    double TOPOINDEX;

    // Huge Resistance Term
    double HUGE_RESIST;  /**< Extermely large resistance term (s/m) */

    // Surface Albedo Parameters
    double ALBEDO_BARE_SOIL;  /**< Broadband albedo of bare soil */

    // Surface Emissivities
    double EMISS_GRND;  /**< Emissivity of bare soil */
    double EMISS_VEG;  /**< Emissivity of vegetation */
    double EMISS_ICE;  /**< Emissivity of bare ice */
    double EMISS_SNOW;  /**< Emissivity of snow */
    double EMISS_H2O;  /**< Emissivity of open water surface */

    // Soil Constraints
    double SOIL_RARC;  /**< Architectural resistance (s/m) of soil when computing soil evaporation via Penman-Monteith eqn */
    double SOIL_RESID_MOIST;  /**< Default residual moisture content (fraction of porosity) of soil column */
    double SOIL_SLAB_MOIST_FRACT;  /**< Moisture content (fraction of porosity) in the soil/rock below the bottom soil layer; this assumes that the soil below the bottom layer has the same texture as the bottom layer. */
    double SOIL_WINDH;  /**< Default wind measurement height over soil (m) */
    double SOIL_RESIST_EXP;
    double SOIL_FROST;

    // Vegetation Parameters
    double VEG_LAI_SNOW_MULTIPLIER;  /**< multiplier to calculate the amount of available snow interception as a function of LAI (m) */
    double VEG_LAI_WATER_FACTOR;  /**< Coefficient multiplied by the LAI to determine the amount of water that can be stored in the canopy */
    double VEG_MIN_INTERCEPTION_STORAGE;  /**< the amount of snow on the canopy that can only be melted off. (m) */
    double VEG_RATIO_DH_HEIGHT;  /**< Ratio of displacement height (m) to vegetation height (m) */
    double VEG_RATIO_RL_HEIGHT;  /**< Ratio of roughness length (m) to vegetation height (m) */

    // Canopy Parameters
    double CANOPY_CLOSURE;  /**< Threshold vapor pressure deficit for stomatal closure (Pa) */
    double CANOPY_RSMAX;  /**< Maximum allowable resistance (s/m) */
    double CANOPY_VPDMINFACTOR;  /**< Minimum allowable vapor pressure deficit factor */

    // Saturation Vapor Pressure Parameters
    double SVP_A0, SVP_A1, SVP_A2, SVP_A3, SVP_A4, SVP_A5, SVP_A6;  /**< constant for saturated vapor pressure curve (kPa) */
    double SVP_B0, SVP_B1, SVP_B2, SVP_B3, SVP_B4, SVP_B5, SVP_B6;  /**< constant for saturated vapor pressure curve (kPa) */
    double SVP_C0, SVP_C1, SVP_C2, SVP_C3, SVP_C4, SVP_C5, SVP_C6;  /**< constant for saturated vapor pressure curve (kPa) */
    double SVP_D0, SVP_D1, SVP_D2, SVP_D3, SVP_D4, SVP_D5, SVP_D6;
    double SVP_FRZ;
    double SVP_RDAIR;
    double SVP1, SVP2, SVP3;
    // Photosynthesis Parameters
    double PHOTO_OMEGA;  /**< single leaf scattering albedo */
    double PHOTO_LAIMAX;  /**< Maximum LAI in nitrogen scaling */
    double PHOTO_LAILIMIT;  /**< Minimum LAI in nitrogen scaling and maximum LAI in PAR computation */
    double PHOTO_LAIMIN;  /**< Minimum LAI in PAR computation */
    double PHOTO_EPAR;  /**< Energy content of PAR [J/mol photons] = (4.6 mol/MJ PAR)^-1 */
    double PHOTO_FCMAX;  /**< Maximum fractional veg cover; (1-FcMax) = min amount of ground visible */
    double PHOTO_FCMIN;  /**< Minimum fractional veg cover; (1-FcMin) = max amount of ground visible */
    double PHOTO_ZENITHMIN;  /**< Check for solar zenith angle > 89 deg */
    double PHOTO_ZENITHMINPAR;  /**< Cosine of the minimum solar zenith angle for photosynthesis to take place */
    double PHOTO_ALBSOIPARMIN;  /**< Minimum soil reflectivity in PAR range */
    double PHOTO_MINMAXETRANS;  /**< Minimum of maximum electron transport rate [10e-12 mol/(m^2 s)] */
    double PHOTO_MINSTOMCOND;  /**< Minimum stomatal conductance [mol H2O/m2s] */
    double PHOTO_FCI1C3;  /**< C3 Plants factor that relate leaf internal CO2 concentration to ambient CO2 concentration */
    double PHOTO_FCI1C4;  /**< C4 Plants factor that relate leaf internal CO2 concentration to ambient CO2 concentration */
    double PHOTO_OX;  /**< OXYGEN CONCENTRATION [MOL(O2) / MOL(AIR)] */
    double PHOTO_KC;  /**< MICHAELIS-MENTEN CONSTANT FOR CO2 AT 25C [MOL(CO2) / MOL(AIR)] */
    double PHOTO_KO;  /**< MICHAELIS-MENTEN CONSTANT FOR O2 AT 25C [MOL(O2) / MOL(AIR)] */
    double PHOTO_EC;  /**< ACTIVATION ENERGY FOR KC [J / MOL] */
    double PHOTO_EO;  /**< ACTIVATION ENERGY FOR KO [J / MOL] */
    double PHOTO_EV;  /**< ACTIVATION ENERGY FOR VCMAX [J / MOL] */
    double PHOTO_ER;  /**< ACTIVATION ENERGY FOR DARK RESPIRATION [J / MOL] */
    double PHOTO_ALC3;  /**< EFFICIENCY OF OF PHOTON CAPTURE */
    double PHOTO_FRDC3;  /**< RATIO OF DARK RESPIRATION TO "PVM" AT 25C for C3 */
    double PHOTO_EK;  /**< = Q10=2 (Collatz et al. 1992) */
    double PHOTO_ALC4;  /**< EFFECTIVE QUANTUM EFFICIENCY */
    double PHOTO_FRDC4;  /**< RATIO OF DARK RESPIRATION TO "PVM" AT 25C for C4 */
    double PHOTO_THETA;  /**< CURVATURE PARAMETER */
    double PHOTO_FRLEAF;  /**< Ratio of canopy leaf respiration to whole plant maintenance respiration */
    double PHOTO_FRGROWTH;  /**< Ratio of plant growth respiration to NPP */

    // Soil Respiration Parameters
    double SRESP_E0_LT;  /**< Lloyd-Taylor E0 parameter [K] */
    double SRESP_T0_LT;  /**< Lloyd-Taylor T0 parameter [K] */
    double SRESP_WMINFM;  /**< minimum soil moisture (fraction) at which soil respiration can occur */
    double SRESP_WMAXFM;  /**< maximum soil moisture (fraction) at which soil respiration can occur */
    double SRESP_WOPTFM;  /**< soil moisture (fraction) at which maximum soil respiration occurs */
    double SRESP_RHSAT;  /**< ratio of soil respiration rate under saturated conditions (w=wmaxFM) to that under optimal conditions (w=woptFM) */
    double SRESP_RFACTOR;  /**< scaling factor to account for other (non-moisture) sources of inhibition of respiration */
    double SRESP_TAULITTER;  /**< Litter pool turnover time [y] */
    double SRESP_TAUINTER;  /**< Intermediate pool turnover time [y] */
    double SRESP_TAUSLOW;  /**< Slow pool turnover time [y] */
    double SRESP_FAIR;  /**< Fraction of respired carbon from litter pool that is lost to atmosphere */
    double SRESP_FINTER;  /**< Fraction of [respired carbon from litter pool that goes to soil] that goes to intermediate pool */
    double SRESP_EXP;
    double SRESP_PSIWILT; /** metric potential for wilting point (m) */
    // Snow Parameters
    double SNOW_MAX_SURFACE_SWE;  /**< maximum depth of the surface layer in water equivalent (m) */
    double SNOW_MAX_LIQUID_FRAC;
    double SNOW_LIQUID_WATER_CAPACITY;  /**< water holding capacity of snow as a fraction of snow-water-equivalent */
    double SNOW_NEW_SNOW_DENSITY;  /**< density of new fallen snow */
    double SNOW_NEW_SNOW_DENS_MAX; /**< new snow density max for Hedstrom and Pomeroy 1998 equation [Warren et al. 1999, Bormann et al. 2013, Maidment Figure 7.2.3] */
    double SNOW_DEPTH_THRES;  /**< Snow depth threshold below which we do not consider the ground flux out of the snowpack in calculating change in cold content (m) */
    double SNOW_DENS_DMLIMIT;  /**< Density limit used in calculation of destructive metamorphism (kg/m^3) */
    double SNOW_DENS_DMLIMIT_FACTOR;  /**< Density limit factor used in calculation of destructive metamorphism (kg/m^3) */
    double SNOW_DENS_MAX_CHANGE;  /**< maximum change in snowfall depth (fraction of swe) */
    double SNOW_DENS_ETA0;  /**< viscosity of snow at T = 0C and density = 0 used in calculation of true viscosity (Ns/m2) */
    double SNOW_DENS_C1;  /**< Constant in snow density computation */
    double SNOW_DENS_C2;  /**< Constant in snow density computation */
    double SNOW_DENS_C3;  /**< Constant in snow density computation */
    double SNOW_DENS_C3_CONST;  /**< Constant in snow density computation */
    double SNOW_DENS_C4;  /**< Constant in snow density computation */
    double SNOW_DENS_C4WET;  /**< Constant in snow density computation */
    double SNOW_DENS_C5;  /**< constant used in snow viscosity calculation, taken from SNTHRM.89 (/C) */
    double SNOW_DENS_C6;  /**< constant used in snow viscosity calculation, taken from SNTHRM.89 (kg/m3) */
    double SNOW_DENS_F;  /**< internal compaction rate coefficient */
    double SNOW_DENS_EXP;  /**< exponent in snow density compaction equation [Bras pg. 257 ]*/
    double SNOW_DENS_DENOM;  /**< denomenator in snow density compaction equation [Bras pg. 257] */
    double SNOW_NEW_SNT_C1; /**< Constant in Sntherm new snow density computation. */
    double SNOW_NEW_SNT_C2; /**< Constant in Sntherm new snow density computation. */
    double SNOW_NEW_SNT_C3; /**< Constant in Sntherm new snow density computation. */
    double SNOW_NEW_BRAS_DENOM;  /**< Constant in Bras new snow density computation. */
    double SNOW_MIN_SWQ_EB_THRES;  /**< Minimum SWQ for which the snowpack energy balance is computed independent of the soil surface temperature */
    double SNOW_A1;  /**< Attenuation coefficient for shortwave in a snowpack. Value and equation taken from Patterson and Hamblin, 1988 */
    double SNOW_A2;  /**< Attenuation coefficient for shortwave in a snowpack. Value and equation taken from Patterson and Hamblin, 1988 */
    double SNOW_L1;  /**< Attenuation coefficient for shortwave in a snowpack. Value and equation taken from Patterson and Hamblin, 1988 (1/m) */
    double SNOW_L2;  /**< Attenuation coefficient for shortwave in a snowpack. Value and equation taken from Patterson and Hamblin, 1988 (1/m) */
    double SNOW_TRACESNOW;  /**< Defines the minimum amount of new snow (mm) which will reset the snowpack albedo to new snow */
    double SNOW_NEW_SNOW_ALB;  /**< Snow albedo curve parameters. */
    double SNOW_ALB_ACCUM_A;  /**< Snow albedo curve parameters. */
    double SNOW_ALB_ACCUM_B;  /**< Snow albedo curve parameters. */
    double SNOW_ALB_THAW_A;  /**< Snow albedo curve parameters. */
    double SNOW_ALB_THAW_B;  /**< Snow albedo curve parameters. */
    double SNOW_CONDUCT;  /**< conductivity of snow (W/mK) */
    double SNOW_PGRAD;
    double SNOW_COMPACT_A;
    double SNOW_COMPACT_B;
    double SNOW_COMPACT_C;
    double SNOW_COMPACT_P;
    double SNOW_COMPACT_DM;
    double SNOW_COMPACT_ETA;
    double SNOW_RADIUS_MIN;
    double SNOW_RADIUS_MAX;
    double SNOW_NEW_RADIUS;
    double SNOW_AGE_FACT;
    double SNOW_AGE_VAPF;
    double SNOW_AGE_FRZF;
    double SNOW_AGE_SOTF;
    double SNOW_NEW_SNOW_COVER;
    double SNOW_COSZEN_B;
    double SNOW_AGE_DIR_VIS;
    double SNOW_AGE_DIR_NIR;
    double SNOW_AGE_DFS_VIS;
    double SNOW_AGE_DFS_NIR;
    double SNOW_NEW_SNOW_VIS;
    double SNOW_NEW_SNOW_NIR;
    double SNOW_RELEASE_FAC;
    double SNOW_RASNOW;     /**< aerodynamic resistance over snow surface (s/m) */
    double SNOW_OMEGAS[MAX_SWBANDS];
    double GLAC_ALBEDO[MAX_SWBANDS];

    // Blowing Snow Parameters
    double BLOWING_KA;  /**< thermal conductivity of air (W/mK) */
    double BLOWING_CSALT;  /**< saltation constant m/s */
    double BLOWING_UTHRESH;  /**< threshold shear velocity m/s */
    double BLOWING_KIN_VIS;  /**< Kinemativ viscosity of air (m2/s) */
    int BLOWING_MAX_ITER;     /**< Max. iterations for numerical integration */
    int BLOWING_K;
    double BLOWING_SETTLING;  /**< Particle settling velocity m/s */
    int BLOWING_NUMINCS;     /**< Number of prob intervals to solve for wind. */

    // Treeline temperature
    double TREELINE_TEMPERATURE;  /**< Number of prob intervals to solve for wind. */

    // Iteration Bracket Widths
    double SNOW_DT;  /**< Used to bracket snow surface temperatures while computing the snow surface energy balance (C) */
    double SURF_DT;  /**< Used to bracket soil surface temperatures while computing energy balance (C) */
    double SOIL_DT;  /**< Used to bracket soil temperatures while solving the soil thermal flux (C) */
    double CANOPY_DT;  /**< Used to bracket canopy air temperatures while computing energy balance (C) */
    double CANOPY_VP;  /**< Used to bracket canopy vapor pressures while computing moisture balance (Pa) */

    // Solar radiation fraction
    double RAD_DIR_F;
    double RAD_VIS_F;

    // Convergence Tolerances
    double TOL_A;
    double TOL_B;
    double MAX_LIMIT;
    double TOL_WETBULB;
    double MAX_ITER_WETBULB;
    double MAX_ITER_MOST;
    double MAX_ITER_OVER;

    // Frozen Soil Parameters
    int FROZEN_MAXITER;

    // Canopy Iterations
    int MAX_ITER_GRND_CANOPY;
    int MAX_ITER_INFLI_RUNOFF;

    // Newton-Raphson Solver Parameters
    int NEWT_RAPH_MAXTRIAL;
    double NEWT_RAPH_TOLX;
    double NEWT_RAPH_TOLF;
    double NEWT_RAPH_R_MAX;
    double NEWT_RAPH_R_MIN;
    double NEWT_RAPH_RELAX1;
    double NEWT_RAPH_RELAX2;
    double NEWT_RAPH_RELAX3;
    double NEWT_RAPH_EPS2;

    // Root-Brent parameters
    int ROOT_BRENT_MAXTRIES;
    int ROOT_BRENT_MAXITER;
    double ROOT_BRENT_TSTEP;
    double ROOT_BRENT_T;
} parameters_struct;

/******************************************************************************
 * @brief   This structure stores the soil parameters for a grid cell.
 *****************************************************************************/
typedef struct {
    bool FS_ACTIVE;                   /**< if TRUE frozen soil algorithm is
                                         active in current grid cell */
    double Ksat[MAX_LAYERS];          /**< saturated hydraulic  conductivity
                                         (mm/day) */
    double Ksat_node[MAX_NODES];
    double Wfc[MAX_LAYERS];           /**< critical moisture level for soil
                                         layer, evaporation is no longer
                                         affected moisture stress in the
                                         soil (mm) */
    double Wfc_node[MAX_NODES];
    double Wpwp[MAX_LAYERS];          /**< soil moisture content at permanent
                                         wilting point (mm) */
    double Wpwp_node[MAX_NODES];
    double Wsat[MAX_LAYERS];
    double Wsat_node[MAX_NODES];
    double AlbedoSat[MAX_SWBANDS];
    double AlbedoDry[MAX_SWBANDS];
    double AlbedoPar;                 /**< soil albedo in PAR range (400-700nm) */
    double b_infilt;                  /**< infiltration parameter */
    double b_dynamic;                 /**< Dynamic VIC heterogeniety parameter for infiltration */
    double bexp[MAX_LAYERS];
    double bexp_node[MAX_NODES];
    double bubble[MAX_LAYERS];        /**< bubbling pressure, HBH 5.15 (cm) */
    double bubble_node[MAX_NODES];    /**< bubbling pressure (cm) */
    double bulk_density[MAX_LAYERS];  /**< soil bulk density (kg/m^3) */
    double bulk_dens_min[MAX_LAYERS]; /**< bulk density of mineral soil (kg/m^3) */
    double bulk_dens_min_node[MAX_NODES];
    double bulk_dens_org[MAX_LAYERS]; /**< bulk density of organic soil (kg/m^3) */
    double cap_drive;                 /**< mean capilary drive (m) for dynamic VIC runoff */
    double Dsat[MAX_LAYERS];          /**< soil moisture diffusion parameter (mm/mm) */
    double Dsat_node[MAX_NODES];
    double depth[MAX_LAYERS];         /**< thickness of each soil moisture layer (m) */
    double dp;                        /**< soil thermal damping depth (m) */
    double dz_soil[MAX_NODES];        /**< thermal node thickness (m) */
    double Zsum_soil[MAX_NODES];      /**< thermal node depth (m) */
    double zc_soil[MAX_NODES];        /**< depth of thermal nodes below soil surface (m) */
    double expt[MAX_LAYERS];          /**< layer-specific exponent n (=3+2/lambda) in Campbell's eqn for hydraulic conductivity, HBH 5.6 */
    double expt_node[MAX_NODES];      /**< node-specific exponent n (=3+2/lambda) in Campbell's eqn for hydraulic conductivity, HBH 5.6 */
    double max_infil;                 /**< maximum infiltration rate */
    double max_moist[MAX_LAYERS];     /**< maximum moisture content (mm) per layer */
    double max_snow_distrib_slope;    /**< Maximum slope of snow depth distribution [m].  This should equal 2*depth_min, where depth_min = minimum snow pack depth below which coverage < 1.  Comment, ported from user_def.h, with questionable units: SiB uses 0.076; Rosemount data imply 0.155cm depth ~ 0.028mm swq. */         
    double porosity[MAX_LAYERS];      /**< porosity (fraction) */
    double porosity_node[MAX_NODES];  /**< porosity (fraction) per node */
    double psi_sat[MAX_LAYERS];
    double psi_sat_node[MAX_NODES];
    double quartz[MAX_LAYERS];        /**< quartz content of soil (fraction of mineral soil volume) */
    double quartz_node[MAX_NODES];
    double organic[MAX_LAYERS];       /**< organic content of soil (fraction of total soil volume) */
    double organic_node[MAX_NODES];
    double resid_moist[MAX_LAYERS];   /**< residual moisture content of soil layer (mm) */
    double resid_moist_node[MAX_NODES];
    double rough;                     /**< soil surface roughness (m) */
    double snow_rough;                /**< snow surface roughness (m) */
    double soil_density[MAX_LAYERS];  /**< soil particle density (kg/m^3) */
    double soil_dens_min[MAX_LAYERS]; /**< particle density of mineral soil (kg/m^3) */
    double soil_dens_min_node[MAX_NODES];
    double soil_dens_org[MAX_LAYERS]; /**< particle density of organic soil (kg/m^3) */
    double *BandElev;                 /**< Elevation of each snow elevation band */
    double *AreaFract;                /**< Fraction of grid cell included in each snow elevation band */
    double *Pfactor;                  /**< Change in Precipitation due to elevation (fract) in each snow elevation band */
    double *Tfactor;                  /**< Change in temperature due to elevation (C) in each snow elevation band */

    double elevation;                 /**< grid cell elevation (m) */
    double lat;                       /**< grid cell central latitude */
    double lng;                       /**< grid cell central longitude */
    double cell_area;                 /**< Area of grid cell (m^2) */
    double time_zone_lng;             /**< central meridian of the time zone */
    unsigned int gridcel;             /**< grid cell number */
    double zwtvmoist_zwt[MAX_LAYERS + 2][MAX_ZWTVMOIST]; /**< zwt values in the zwt-v-moist curve for each layer */
    double zwtvmoist_moist[MAX_LAYERS + 2][MAX_ZWTVMOIST]; /**< moist values in the zwt-v-moist curve for each layer */
    double slope;
    double aspect;
    double ehoriz;
    double whoriz;

} soil_con_struct;

/******************************************************************************
 * @brief   This structure stores information about the vegetation coverage of
 *          the current grid cell.
 *****************************************************************************/
typedef struct {
    double albedo[MONTHS_PER_YEAR];   /**< climatological vegetation albedo
                                         (fraction) */
    double displacement[MONTHS_PER_YEAR]; /**< climatological vegetation
                                             displacement height (m) */
    double fcanopy[MONTHS_PER_YEAR];    /**< fractional area covered by plant
                                                            canopy (fraction) */
    double LAI[MONTHS_PER_YEAR];        /**< leaf area index */
    double SAI[MONTHS_PER_YEAR];        /**< stem area index */
    double roughness[MONTHS_PER_YEAR]; /**< climatological vegetation
                                          roughness length (m) */
    double *CanopLayerBnd;  /**< Upper boundary of each canopy layer,
                               expressed as fraction of total LAI */
    double Cv;              /**< fraction of vegetation coverage */
    double root[MAX_NODES]; /**< percent of roots in each soil layer
                                (fraction) */
    int veg_class;          /**< vegetation class id number */
    size_t vegetat_type_num; /**< number of vegetation types in the grid
                                cell */
    double a;               /**< Empirical parameter a in eqa(2) */
    double b;               /**< Empirical parameter b in eqa(2) */
    double d;               /**< Maximum root depth (m) */
    int BandIndex;
    bool IS_GLAC;
    size_t Nroot;
} veg_con_struct;

/******************************************************************************
 * @brief   This structure stores parameters for individual vegetation types.
 *****************************************************************************/
typedef struct {
    double albedo[MONTHS_PER_YEAR];     /**< vegetation albedo (added for full
                                                           energy) (fraction) */
    double displacement[MONTHS_PER_YEAR]; /**< climatological vegetation
                                             displacement height (m) */
    double fcanopy[MONTHS_PER_YEAR];    /**< fractional area covered by plant
                                                            canopy (fraction) */
    double LAI[MONTHS_PER_YEAR];        /**< leaf area index */
    double SAI[MONTHS_PER_YEAR];        /**< stem area index */
    double roughness[MONTHS_PER_YEAR]; /**< climatological vegetation
                                          roughness length (m) */

    double reflleaf[MAX_SWBANDS];
    double reflstem[MAX_SWBANDS];
    double transleaf[MAX_SWBANDS];
    double transstem[MAX_SWBANDS];
    double alpha_canopy;
    double Canopy_Upper;                /**< top of canopy (m) */
    double Canopy_Lower;                /**< bottom of canopy (m) */
    double Canopy_Radius;
    double COI;
    double c_biomass;
    double d_leaf;
    size_t NVegLibTypes;   /**< number of vegetation classes defined in
                              library */
    double rad_atten;      /**< radiation attenuation due to canopy,
                              default = 0.5 (N/A) */
    double rmax;           /**< Maximal stomatal resistance (s/m) */
    double rmin;           /**< minimum stomatal resistance (s/m) */
    double trunk_ratio;    /**< ratio of trunk height to tree height,
                              default = 0.2 (fraction) */
    double wind_atten;     /**< wind attenuation through canopy,
                              default = 0.5 (N/A) */
    double RGL;            /**< Value of solar radiation below which there
                              will be no transpiration (ranges from
                              ~30 W/m^2 for trees to ~100 W/m^2 for crops) */
    double m_vpd;
    double T_opt_trans;

    unsigned short int veg_class; /**< vegetation class reference number */
    // Carbon terms
    char Ctype;            /**< Photosynthetic pathway; 0 = C3; 1 = C4 */
    double CO2Specificity; /**< CO2 specificity at 25 deg C (mol(CO2)/m2s)
                              (C4 plants) */
    double LightUseEff;    /**< Light-use efficiency (mol(CO2)/mol(photons)) */
    double MaxCarboxRate;  /**< maximum carboxlyation rate at 25 deg C
                              (mol(CO2)/m2s) */
    double MaxETransport;  /**< maximum electron transport rate at 25 deg C
                              (mol(CO2)/m2s) (C3 plants) */
    double NPPfactor_sat;  /**< photosynthesis multiplier (fraction of
                              maximum) when top soil layer is saturated */
    bool NscaleFlag;       /**< TRUE = nitrogen-scaling factors are
                              applicable to this veg class */
    double Wnpp_inhib;     /**< moisture level (fraction of maximum moisture)
                              above which photosynthesis experiencing
                              saturation inhibition, i.e. too wet for optimal
                              photosynthesis; only applies to top soil layer */
} veg_lib_struct;

/******************************************************************************
 * @brief   This structure stores vegetation parameter forcing data for each
 * model time step for a single veg tile.  Each array stores the values for the
 * SNOW_STEPs during the current model step and the value for the entire model
 * step.  The latter is referred to by array[NR].  Looping over the SNOW_STEPs
 * is done by for (i = 0; i < NF; i++)
 *****************************************************************************/
typedef struct {
    double *albedo;       /**< vegetation albedo (fraction) */
    double *displacement; /**< vegetation displacement height (m) */
    double *fcanopy;      /**< fractional area covered by plant canopy
                             (fraction) */
    double *LAI;          /**< leaf area index (m2/m2) */
    double *SAI;
    double *roughness;    /**< vegetation roughness length (m) */
} veg_hist_struct;

/******************************************************************************
 * @brief   This structure stores the forcing data for each model
 * time step for a single grid cell.  Each array stores the values for the
 * SNOW_STEPs during the current model step and the value for the entire model
 * step.  The latter is referred to by array[NR].  Looping over the SNOW_STEPs
 * is done by for (i = 0; i < NF; i++)
 *****************************************************************************/
typedef struct {
    double *air_temp;   /**< air temperature (C) */
    double *Catm;       /**< atmospheric CO2 mixing ratio (mol CO2/ mol air) */
    double *channel_in; /**< incoming channel inflow for time step (mm) */
    double *coszen;     /**< cosine of the solar zenith angle */
    double *density;    /**< atmospheric density (kg/m^3) */
    double *fdir;       /**< fraction of incoming shortwave that is direct (fraction) */
    double *longwave;   /**< incoming longwave radiation (W/m^2) (net incoming
                                           longwave for water balance model) */
    double *Qair;
    double *par;        /**< incoming photosynthetically active radiation () */
    double *prec;       /**< average precipitation in grid cell (mm) */
    double *snowf;      /**< snowfall partitioned from precipitation (mm) */
    double *rainf;      /**< rainfall partitioned from precipitation (mm) */
    double *pressure;   /**< atmospheric pressure (kPa) */
    double *shortwave;  /**< incoming shortwave radiation (W/m^2) */
    double *vp;         /**< atmospheric vapor pressure (kPa) */
    double *rel_humid;  /**< atmospheric vapor pressure deficit (kPa) */
    double *wind;       /**< wind speed (m/s) */
} force_data_struct;

/******************************************************************************
 * @brief   This structure stores information about the time and date of the
 *          current time step.
 *****************************************************************************/
typedef struct {
    unsigned short int day;         /**< current day */
    unsigned short int day_in_year; /**< julian day in year */
    unsigned short int month;       /**< current month */
    int year;                       /**< current year */
    unsigned int dayseconds;        /**< seconds since midnight */
} dmy_struct;                       /**< array of length nrec created */


/******************************************************************************
 * @brief   This structure stores soil variables for the complete soil column
 *          for each grid cell.
 *****************************************************************************/
typedef struct {
    // State variables
    double Ra_over[3];                 /**< The (stability-corrected) aerodynamic
                                          resistance (s/m) that was actually used
                                          in flux calculations.
                                          [0] = surface (bare soil, non-overstory veg, or snow pack)
                                          [1] = overstory */
    double Ra_sub[3];
    double Ra_grnd[3];
    double Ra_evap;                    /**< ground surface resistance [s/m] to evaporation */
    double Ra_leaf;
    double asat;                       /**< saturated area fraction */
    double CLitter;                    /**< carbon storage in litter pool [gC/m2] */
    double CInter;                     /**< carbon storage in intermediate pool [gC/m2] */
    double CSlow;                      /**< carbon storage in slow pool [gC/m2] */
    double rootmoist;                  /**< total of layer.moist over all layers
                                          in the root zone (mm) */
    double wetness;                    /**< average of
                                          (layer.moist - Wpwp)/(porosity*depth - Wpwp)
                                          over all layers (fraction) */
    double zwt;                        /**< average water table position [cm] - using lowest unsaturated layer */
    double zwt_lumped;                 /**< average water table position [cm] - lumping all layers' moisture together */
    double VPcanopy;
//    double pack_transp;
    double ice[MAX_NODES];             /**< ice content of the frozen sublayer (mm) */
    double liq[MAX_NODES];             /**< liq content of the frozen sublayer (mm) */
    double moist[MAX_NODES];           /**< moisture content of the unfrozen sublayer
                                            (mm) */
    double soil_T[MAX_NODES];
    double SnowFrac[MAX_SNOWS];
    double zc_node[MAX_NODES + MAX_SNOWS];     /**< depth of thermal nodes (m) */
    double dz_node[MAX_NODES + MAX_SNOWS];     /**< the thickness of each layer (m) */
    double Zsum_node[MAX_NODES + MAX_SNOWS];   /**< depth of bottom of each thermal node (m) */
    int PhaseChange[MAX_NODES + MAX_SNOWS];
    // Fluxes
    double baseflow;                   /**< baseflow from current cell (mm/TS) */
    double runoff;                     /**< runoff from current cell (mm/TS) */
    double soil_inflow;                /**< moisture that reaches the top of
                                            the soil column (mm) */
    double evap;
    double vapor_grnd;
    double conden_grnd;
    double snowfrost;
    double snow_sublim;
    // Canopy terms
    double transp;
    double canopyevap;
    double canopydew;
    double canopyfrost;
    double canopy_sublim;
    double canopy_vapor;
    // Soil terms
    double esoil;                      /**< soil evaporation from soil layer (mm) */
    double dewsoil;                    /**< evapotranspiration from soil layer (mm) */
    double total_transp;               /**< accumulated soil water transpiration factor (0 to 1) */
    double soil_excess;
    double soil_frost[MAX_NODES];      /**< frost content of the frozen sublayer */
    double soil_transp[MAX_NODES];     /**< transpiration from soil layer (mm) */
    double transp_fact[MAX_NODES];     /**< soil layer transpiration facter */
    double diffusivity[MAX_NODES];     /**< soil water diffusivity [m2/s] */
    double conductivity[MAX_NODES];    /**< soil hydraulic conductivity [m/s] */
    // Carbon terms
    double RhLitter;                   /**< soil respiration from litter pool [gC/m2] */
    double RhLitter2Atm;               /**< soil respiration from litter pool [gC/m2] that goes to atmosphere */
    double RhInter;                    /**< soil respiration from intermediate pool [gC/m2] */
    double RhSlow;                     /**< soil respiration from slow pool [gC/m2] */
    double RhTot;                      /**< total soil respiration over all pools [gC/m2] (=RhLitter2Atm+RhInter+RhSlow) */
} cell_data_struct;

/******************************************************************************
 * @brief   This structure stores energy balance components, and variables used
 *          to solve the thermal fluxes through the soil column.
 *****************************************************************************/
typedef struct {
    // State variables
    bool FrozenGrnd;                   /**< TRUE = frozen soil present */
    bool FrozenOver;                   /**< TRUE = frozen canopy present */

    double kappa_node[MAX_NODES + MAX_SNOWS]; /**< thermal conductivity of the soil thermal nodes (W/m/K) */
    double Cs_node[MAX_NODES + MAX_SNOWS];   /**< heat capacity of the soil thermal nodes (J/m^3/K) */    
    double T[MAX_NODES + MAX_SNOWS];         /**< thermal node temperatures (C) */
    double kappa_int[MAX_NODES + MAX_SNOWS]; /**< thermal conductivity used for interface between nodes (W/m/K) */
    double fact[MAX_NODES + MAX_SNOWS];        /**< thermal diffusivity of each soil thermal node (m^2/s) */
    double fn[MAX_NODES + MAX_SNOWS];          /**< Fourier number for each soil thermal node */
    double snow_flux;            /**< heat flux into or out of the snowpack (Wm-2) */
    double soil_flux;            /**< heat flux into or out of the soil column (Wm-2) */
    double flux_slope;          /**< slope of the soil thermal flux with respect to soil surface temperature (Wm-2/C) */
    double Nthaw;
    double Nfrost;
    double fdepth;
    double tdepth;
    double Tcanopy;              /**< temperature of the canopy */
    double Tair;                 /**< temperature of the canopy air */
    double Tsurf;                /**< temperature of the understory */
    double Tgrnd;                /**< temperature of the bare ground */
    double Tupper;               /**< temperature of upper soil/snow layer */
    double Tlower;               /**< temperature of lower soil/snow layer */
    double PsyCh_grnd;
    double PsyCh_canopy;

    // Fluxes
    double APAR_sunlit;
    double APAR_shade;  
    double advection;            /**< advective flux (Wm-2) */
    double AdvectSub;
    double AdvectGrnd;
    double AdvectOver;
    double grnd_flux;            /**< ground heat flux (Wm-2) */
    double GroundSub;
    double GroundGrnd;
    double latent;               /**< net latent heat flux (Wm-2) */
    double LatentSub;
    double LatentGrnd;
    double sensible;             /**< net sensible heat flux (Wm-2) */
    double SensibleSub;
    double SensibleGrnd;
    double SensibleOver;
    double LatentCanopy;
    double LatentTransp;
    double LatentVapOver;
    double LatentVapGrnd;

    double ReflShortSurf;
    double ReflShortGrnd;
    double ReflShortSub;
    double EmissLongSub;
    double EmissLongGrnd;
    double EmissLongSurf;
    double ReflectVeg[MAX_SWBANDS];
    double TransmitVeg[MAX_SWBANDS];
    double fusion_fact[MAX_NODES + MAX_SNOWS];
    double fusion_flux[MAX_NODES + MAX_SNOWS];

    double longwave;             /**< net longwave flux (Wm-2) */
    double NetLongSurf;
    double NetLongGrnd;          /**< net longwave radiation to the atmosphere (W/m^2) */
    double NetLongOver;          /**< net longwave radiation from the overstory (W/m^2) */
    double NetLongSub;         /**< net longwave radiation from the understory (W/m^2) */
    double shortwave;
    double NetShortSurf;         /**< net shortwave to the atmosphere */
    double NetShortGrnd;         /**< net shortwave penetrating snowpack */
    double NetShortSub;        /**< net shortwave radiation from the understory (W/m^2) */
    double NetShortSoil;      /**< net shortwave radiation to the soil (W/m^2) */
    double NetShortSnow;      /**< net shortwave radiation to the snow (W/m^2) */

    double AlbedoSoilDir[MAX_SWBANDS];
    double AlbedoSoilDfs[MAX_SWBANDS];
    double AlbedoSnowDir[MAX_SWBANDS];
    double AlbedoSnowDfs[MAX_SWBANDS];
    double AlbedoGrndDir[MAX_SWBANDS];
    double AlbedoGrndDfs[MAX_SWBANDS];
    double AbsSubDir[MAX_SWBANDS];
    double AbsSubDfs[MAX_SWBANDS];
    double ShortDir2Dir[MAX_SWBANDS];
    double ShortDfs2Dir[MAX_SWBANDS];
    double ShortDir2Dfs[MAX_SWBANDS];
    double ShortDfs2Dfs[MAX_SWBANDS];
    double AlbedoSurfDir[MAX_SWBANDS];
    double AlbedoSurfDfs[MAX_SWBANDS];
    double ReflGrndDir[MAX_SWBANDS];
    double ReflGrndDfs[MAX_SWBANDS];
    double ReflSubDir[MAX_SWBANDS];
    double ReflSubDfs[MAX_SWBANDS];

} energy_bal_struct;

/******************************************************************************
 * @brief   This structure stores vegetation variables for each vegetation type
 *          in a grid cell.
 *****************************************************************************/
typedef struct {
    // State variables
    double albedo;              /**< current vegetation albedo (fraction) */
    double displacement;        /**< current vegetation displacement height
                                   (m) */
    double fcanopy;             /**< current fractional area of plant canopy
                                   (fraction) */
    double LAI;                 /**< current leaf area index (m2/m2) */
    double SAI;                 /**< current stem area index (m2/m2) */
    double NetLAI;
    double NetSAI;
    double roughness;           /**< current vegetation roughness length
                                   (m) */
    double Wdew;                /**< maximum intercepted water per unit lai+sai (mm) */
    double MaxSnowInt;          /**< maximum canopy capacity for snow interception [mm] */
    double MaxRainInt;          /**< maximum canopy capacity for rain interception [mm] */
    double wetFrac;
    double leaf_sun;
    double leaf_shade;

    // Fluxes
    double RainThroughFall;     /**< rain that reaches the ground through
                                    the canopy (mm/TS) */
    double SnowThroughFall;     /**< snow that reaches the ground through
                                    the canopy (mm/TS) */
    double RainDrip;
    double SnowDrip;
    double int_rain;            /**< rain intercepted on canopy (mm) */
    double int_snow;            /**< snow intercepted on canopy (mm) */
    double RS_sunlit;           /**< sunlit leaf stomatal resistance [s/m] */
    double RS_shade;            /**< shaded leaf stomatal resistance [s/m] */
    // Carbon terms - states
    double AnnualNPP;           /**< running total annual NPP [gC/m2] */
    double AnnualNPPPrev;       /**< total annual NPP from previous year
                                   [gC/m2] */
    double Ci;                  /**< whole-canopy leaf-internal CO2 mixing
                                   ratio (mol CO2/mol air) */
    double *CiLayer;            /**< array of per-layer leaf-internal CO2
                                   mixing ratio (mol CO2/mol air) */
    double NPPfactor;           /**< whole-canopy photosynthesis multiplier
                                   to account for inhibition separate from
                                   stomatal resistance */
    double *NscaleFactor;       /**< array of per-layer nitrogen scaling
                                   factors */
    double rc;                  /**< whole-canopy stomatal resistance (s/m) */
    double *rsLayer;            /**< array of per-layer stomatal resistance
                                   (s/m) */

    // Carbon terms - fluxes
    double aPAR;                /**< whole-canopy absorbed PAR
                                   (mol(photons)/m2 leaf area s) */
    double *aPARLayer;          /**< array of per-layer absorbed PAR
                                   (mol(photons)/m2 leaf area s) */
    double GPP;                 /**< whole-canopy gross assimilation
                                   (photosynthesis) (umol(CO2)/m2s) */
    double Litterfall;          /**< flux of carbon from living biomass to
                                   litter pool [gC/m2] */
    double NPP;                 /**< net primary productivity (= GPP - Raut)
                                   (umol(CO2)/m2s) */
    double Raut;                /**< total plant respiration (= Rmaint
                                 + Rgrowth) (umol(CO2)/m2s) */
    double Rdark;               /**< whole-canopy 'dark' respiration
                                   (umol(CO2)/m2s) */
    double Rgrowth;             /**< growth respiration ( = (GPP-Rmaint)
                                 * FRGrowth/(1+FRGrowth) ) (umol(CO2)/m2s) */
    double Rmaint;              /**< plant maintenance respiration
                                   (= Rdark/FRLeaf) (umol(CO2)/m2s) */
    double Rphoto;              /**< whole-canopy photorespiration
                                   (umol(CO2)/m2s) */
} veg_var_struct;

/******************************************************************************
 * @brief   This structure stores snow pack variables needed to run the snow
 *          model.
 *****************************************************************************/
typedef struct {
    // State variables
    size_t Nsnow;                   /**< Number of snow layers in the model */
    double albedo;                  /**< snow surface albedo (fraction) */
    double coverage;                /**< fraction of snow band that is covered with snow */
    double density;                 /**< snow density (kg/m^3) */
    double dz_snow[MAX_SNOWS];        /**< each snow pack depth (m) */
    double delta_depth;             /**< snow depth increasing rate [m/s] due to snowfall */
    double Zsum_snow[MAX_SNOWS];
    double zc_snow[MAX_SNOWS];      /**< depth of snow thermal nodes (m) */
    double snow_thresholds[MAX_SNOWS]; /**< snow depth thresholds for layer remobilization (m) */
    double glac_excess;
    double snow_depth;              /**< snow depth (m) */
    double snow_canopy;             /**< amount of snow and water on canopy (mm) */
    double SnowAge;
    double old_swq;
    double pack_T[MAX_SNOWS];
    double pack_ice[MAX_SNOWS];
    double pack_liq[MAX_SNOWS];     /**< liquid water content of the snow pack (m) */
    double pack_outflow[MAX_SNOWS];
    double theta_ice[MAX_SNOWS];
    double theta_liq[MAX_SNOWS];
    double porosity[MAX_SNOWS];
    double new_snow_density;        /**< bulk density of snowfall [kg/m3] */
    double pack_melt;
    double pack_transp;
    double pack_comb;
    double swq;             /**< snow water equivalent of the entire pack (mm) */
} snow_data_struct;

/******************************************************************************
 * @brief   This structure stores all variables needed to solve, or save
 *          solututions for all versions of this model.
 *****************************************************************************/
typedef struct {
    cell_data_struct  *cell;     /**< Stores soil layer variables */
    energy_bal_struct *energy;   /**< Stores energy balance variables */
    snow_data_struct  *snow;     /**< Stores snow variables */
    veg_var_struct    *veg_var;  /**< Stores vegetation variables */
} all_vars_struct;

/******************************************************************************
 * @brief   This structure stores all variables needed to solve, or save
 *          solututions for all versions of this model.
 *****************************************************************************/
typedef struct {
    // State variables
    double S_hill;
    double out_hill;
    double q_surf;
    double delta_S_hill;

} rout_data_struct;

#endif
