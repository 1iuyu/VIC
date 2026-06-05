/******************************************************************************
 * @section DESCRIPTION
 *
 * Definition header file
 *****************************************************************************/

#ifndef VIC_DEF_H
#define VIC_DEF_H

#define _DEFAULT_SOURCE
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
#include <lapacke.h>
#include <pwd.h>
#include <sys/time.h>
#include <sys/types.h>

#include "vic_physical_constants.h"
#include "vic_log.h"

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
#define MAX_NODES       50     /**< maximum number of nodes = MAX_SNOWS + MAX_SOILS + 1 */
#define MAX_FRONTS      3      /**< maximum number of freezing and thawing front depths to store */
#define MAX_SNOWS       3      /**< maximum number of snowpack layers */
#define MAX_SOILS       45     /**< maximum number of soil nodes */
#define MAX_SWBANDS     2      /**< maximum number of solar radiation wave bands */
#define MAX_CANOPYS     10     /**< maximum number of canopy layers for radiative transfer */
#define MAX_HRUS        50     /**< maximum number of hydrological response units */

/***** Define minimum values for model parameters *****/
#define MINSOILDEPTH    0.001  /**< Minimum layer depth with which model can work (m) */
#define MIN_FCANOPY    0.0001  /**< Minimum allowable canopy fraction */
#define MIN_SNOW_WETFRAC 0.01  /**< Minimum fraction of snow depth to be considered wet */
#define MIN_SOILMOIST    0.01  /**< Minimum soil moisture content */
#define MIN_VEG_LAI 0.05       /**< */
#define MIN_TOL_LAI 0.001      /**< Minimum allowable leaf area index for transpiration calculations */

/***** Define minimum and maximum values for model timesteps *****/
#define MIN_SUBDAILY_STEPS_PER_DAY  4
#define MAX_SUBDAILY_STEPS_PER_DAY  1440

#ifndef SNOW
#define RAIN 0
#define SNOW 1
#endif

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

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
 * @brief   Aerodynamic Resistance options
 *****************************************************************************/
enum
{
    AR_ZENG,
    AR_MEIER
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
 * @brief   Photosynthetic pathways
 *****************************************************************************/
enum
{
    PHOTO_C3,
    PHOTO_C4
};

/******************************************************************************
 * @brief   Soil water retention curve parametrizations
 *          Defines which part is called.
 *****************************************************************************/
enum 
{
    MATRIC_FLAG,   // θ → ψ
    MOIST_FLAG,    // ψ → θ
    DERIV_FLAG,    // dψ/dθ
    CONDUCT_FLAG   // θ → k
};

/******************************************************************************
 * @brief   Calculation of canopy type options
 *****************************************************************************/
enum 
{
    SUNLIT,
    SHADE
};

/******************************************************************************
 * @brief   Calculation of saturation vapor pressure options
 *****************************************************************************/
typedef enum {

    ESAT = 1 << 0,  // 1
    QSAT = 1 << 1,  // 2
    ESDT = 1 << 2,  // 4
    QSDT = 1 << 3,  // 8

} svp_flag;

/***** Data Structures *****/

/******************************************************************************
 * @brief   This structure stores model options.
 *****************************************************************************/
typedef struct {
    // simulation modes
    bool ACTIVE_LAYER;
    bool CARBON;         /**< TRUE = simulate carbon cycling processes;
                            FALSE = no carbon cycling (default) */
    bool CONTINUEONERROR; /**< TRUE = VIC will continue to run after a cell has an error */
    bool CORRPREC;       /**< TRUE = correct precipitation for gage undercatch */
    bool BIOMASST;          /**< TRUE = Use biomass code to calculate canopy heat capacity and absorbed radiation */
    size_t Nlayer;       /**< Number of layers in model */
    size_t Nswband;      /**< Number of waveband in model */
    size_t Nfrost;       /**< Number of frost layer in model */
    size_t MAX_HRU;      /**< maximum number of hydrological response units */
    bool NOFLUX;         /**< TRUE = Use no flux lower bondary when computing
                            soil thermal fluxes */
    size_t NVEGTYPES;    /**< number of vegetation types in veg_param file */
    unsigned short int SNOW_DENSITY;   /**< DENS_BRAS: Use algorithm of Bras, 1990; DENS_SNTHRM: Use algorithm of SNTHRM89 adapted for 1-layer pack */
    size_t SNOW_BAND;    /**< Number of elevation bands over which to solve the
                            snow model */
    // glacier options
    int  GLACIER_ID;        /* An index indicating which veg class in the vegetation library contains glacier information. */
    bool TFALLBACK;      /**< TRUE = when any temperature iterations fail to converge,
                                   use temperature from previous time step; the number
                                   of instances when this occurs will be logged and
                                   reported at the end of the cell's simulation
                            FALSE = when iterations fail to converge, report an error
                                    and abort simulation for current grid cell
                            Default = TRUE */
    // input options
    unsigned short int AERO_RESIST;    /**< AR_ZENG = use Zeng et al. (2005) method to calculate canopy aerodynamic resistance
          AR_MEIER = use Meier et al. (2017) method to calculate canopy aerodynamic resistance */
    unsigned short int CANOPY_INTERCEP;
    unsigned short int SNOW_AGING;    /**< BATS and SNICAR */
    unsigned short int GRID_DECIMAL; /**< Number of decimal places in grid file extensions */
    unsigned short int FCAN_SRC;       /**< FROM_VEGLIB = use fcanopy values from veg library file
                                          FROM_VEGPARAM = use fcanopy values from the veg param file */
    unsigned short int LAI_SRC;        /**< FROM_VEGLIB = use LAI values from veg library file
                                          FROM_VEGPARAM = use LAI values from the veg param file */
    unsigned short int SAI_SRC;        /**< FROM_VEGLIB = use SAI values from veg library file
                                          FROM_VEGPARAM = use SAI values from the veg param file */
    bool PARAM_FROM_SOIL; /**< TRUE = bulk density and soil density (particle density) read from soil parameter file; otherwise set to 0.0 */
    bool ROUT;            /**< TRUE = */
    // state options
    unsigned short int STATE_FORMAT;  /**< TRUE = model state file is binary (default) */
    bool INIT_STATE;     /**< TRUE = initialize model state from file */
    bool SAVE_STATE;     /**< TRUE = save state file */

    // output options
    size_t Noutstreams;  /**< Number of output stream */
    bool   PRT_HEADER;     /* TRUE = insert header at beginning of output file; FALSE = no header */

} option_struct;

/******************************************************************************
 * @brief   This structure stores all model run global parameters.
 *****************************************************************************/
typedef struct {
    double resolution;             /**< Model resolution (degrees) */
    double step_dt;                /**< Time step in seconds */
    size_t model_steps_per_day;    /**< Number of model timesteps per day */
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
    // Vegetation Parameters
    double VEG_LAI_SNOW_MULTIPLIER;  /**< multiplier to calculate the amount of available snow interception as a function of LAI (m) */
    double VEG_LAI_WATER_FACTOR;  /**< Coefficient multiplied by the LAI to determine the amount of water that can be stored in the canopy */
    double VEG_RATIO_DH_A;  /**< Ratio of displacement height (m) to vegetation height (m). Canopy_Upper <= 1.0 */
    double VEG_RATIO_DH_B;  /**< Ratio of displacement height (m) to vegetation height (m). Canopy_Upper > 1.0 */
    double VEG_RATIO_RL_A;  /**< Ratio of roughness length (m) to vegetation height (m). Canopy_Upper <= 1.0 */
    double VEG_RATIO_RL_B;  /**< Ratio of roughness length (m) to vegetation height (m). Canopy_Upper > 1.0 */
    // Canopy Parameters
    double CANOPY_CLOSURE;  /**< Threshold vapor pressure deficit for stomatal closure (Pa) */
    double CANOPY_RSMAX;  /**< Maximum allowable resistance (s/m) */
    // Saturation Vapor Pressure Parameters
    double SVP_A0, SVP_A1, SVP_A2, SVP_A3, SVP_A4, SVP_A5, SVP_A6, SVP_A7, SVP_A8;  /**< constant for saturated vapor pressure curve (kPa) */
    double SVP_B0, SVP_B1, SVP_B2, SVP_B3, SVP_B4, SVP_B5, SVP_B6, SVP_B7, SVP_B8;  /**< constant for saturated vapor pressure curve (kPa) */
    double SVP_C0, SVP_C1, SVP_C2, SVP_C3, SVP_C4, SVP_C5, SVP_C6, SVP_C7, SVP_C8;  /**< constant for saturated vapor pressure curve (kPa) */
    double SVP_D0, SVP_D1, SVP_D2, SVP_D3, SVP_D4, SVP_D5, SVP_D6, SVP_D7, SVP_D8;  /**< constant for saturated vapor pressure curve (kPa) */
    double SVP_FRZ;
    double SVP_RDAIR;
    // Photosynthesis Parameters
    double PHOTO_MINCONDUCT;  /**< Minimum stomatal conductance [mol H2O/m2s] */
    double PHOTO_RSMAX;  /**< Maximum stomatal resistance (s/m) */
    double PHOTO_LRESC3;  /**< C3 Plants factor for scalar constant of leaf respiration with Vcmax */
    double PHOTO_LRESC4;  /**< C4 Plants factor for scalar constant of leaf respiration with Vcmax */
    double PHOTO_MAXCS;
    double PHOTO_OX;  /**< OXYGEN CONCENTRATION [MOL(O2) / MOL(AIR)] */
    double PHOTO_CX;
    double PHOTO_KC;  /**< MICHAELIS-MENTEN CONSTANT FOR CO2 AT 25C [MOL(CO2) / MOL(AIR)] */
    double PHOTO_KO;  /**< MICHAELIS-MENTEN CONSTANT FOR O2 AT 25C [MOL(O2) / MOL(AIR)] */
    double PHOTO_CP;  /**< CO2 compensation point at 25degC at present day O2 [mol/mol] */
    double PHOTO_EC;  /**< ACTIVATION ENERGY FOR KC [J / MOL] */
    double PHOTO_EO;  /**< ACTIVATION ENERGY FOR KO [J / MOL] */
    double PHOTO_EP;  /**< ACTIVATION ENERGY FOR CP [J / MOL] */
    double PHOTO_EL;  /**< ACTIVATION ENERGY FOR LMR [J / MOL] */
    double PHOTO_EV;  /**< ACTIVATION ENERGY FOR VCMAX [J / MOL] */
    double PHOTO_EJ;  /**< ACTIVATION ENERGY FOR JMAX [J / MOL] */
    double PHOTO_ET;  /**< ACTIVATION ENERGY FOR TPU [J / MOL] */
    double PHOTO_DV;
    double PHOTO_DJ;
    double PHOTO_DT;
    double PHOTO_FTPU;
    double PHOTO_FKP;
    double PHOTO_FNPS;
    double PHOTO_LMRHD;
    double PHOTO_LMRSE;
    double PHOTO_CROOT;
    double PHOTO_ER;  /**< ACTIVATION ENERGY FOR DARK RESPIRATION [J / MOL] */
    double PHOTO_FNR; /**< Mass ratio of total Rubisco molecular mass to N in Rubisco */
    double PHOTO_SACT;
    // Surface roughness constants
    double ROUGH3; 
    double ROUGH_BETA;
    double ROUGH_NU;
    double SNOW_ROUGH;
    double SOIL_ROUGH;
    double GLAC_ROUGH;
    double SOIL_RROOT;
    double SOIL_RHOROOT;
    // Snow Parameters
    double SNOW_MAX_SURFACE_SWE;  /**< maximum depth of the surface layer in water equivalent (m) */
    double SNOW_MAX_LIQUID_FRAC;
    double SNOW_LIQUID_WATER_CAPACITY;  /**< water holding capacity of snow as a fraction of snow-water-equivalent */
    double SNOW_NEW_SNOW_DENSITY;  /**< density of new fallen snow */
    double SNOW_NEW_SNOW_DENS_MAX; /**< new snow density max for Hedstrom and Pomeroy 1998 equation [Warren et al. 1999, Bormann et al. 2013, Maidment Figure 7.2.3] */
    double SNOW_NEW_SNT_C1; /**< Constant in Sntherm new snow density computation. */
    double SNOW_NEW_SNT_C2; /**< Constant in Sntherm new snow density computation. */
    double SNOW_NEW_SNT_C3; /**< Constant in Sntherm new snow density computation. */
    double SNOW_NEW_BRAS_DENOM;  /**< Constant in Bras new snow density computation. */
    double SNOW_CONDUCT;  /**< conductivity of snow (W/mK) */
    double SNOW_PGRAD;
    double SNOW_COMPACT_A;
    double SNOW_COMPACT_B;
    double SNOW_COMPACT_C;
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
    double SNOW_BETADS;
    double SNOW_BETAIS;
    double SNOW_OMEGAS[MAX_SWBANDS];
    double GLAC_ALBEDO[MAX_SWBANDS];
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
    // Crank Nicholson factor between 0 and 1
    double CN_FACTOR;
} parameters_struct;

/******************************************************************************
 * @brief   This structure stores the soil parameters for a grid cell.
 *****************************************************************************/
typedef struct {

    size_t Nbedrock;                  /**< Number of thermal nodes in the model */
    double Ksat_node[MAX_SOILS];      /**< saturated hydraulic conductivity (m/s) */          
    double Wpwp_node[MAX_SOILS];      /**< soil moisture content at permanent
                                           wilting point (m3/m3) */
    double Wsat_node[MAX_SOILS];      /** soil moisture content at saturation (m3/m3) */  
    double AlbedoSat[MAX_SWBANDS];    /** soil albedo at saturation */
    double AlbedoDry[MAX_SWBANDS];      
    double AlbedoPar;                 /**< soil albedo in PAR range (400-700nm) */
    double alpha_node[MAX_SOILS];     /**< retention shape parameter in van Genuchten equation [1/m] */
    double b_infilt;                  /**< infiltration parameter */
    double b_dynamic;                 /**< Dynamic VIC heterogeniety parameter for infiltration */
    double expt_node[MAX_SOILS];      /**< layer-specific exponent n in Campbell or van Genuchten eqn */
    double bubble_node[MAX_SOILS];    /**< bubbling pressure (cm) */
    double bulk_dens_min[MAX_LAYERS]; /**< bulk density of mineral soil (kg/m^3) */
    double bulk_dens_node[MAX_SOILS]; /**< soil bulk density (kg/m^3) */
    double bulk_dens_org[MAX_LAYERS]; /**< bulk density of organic soil (kg/m^3) */
    double capil_drive;               /**< mean capilary drive (m) for dynamic VIC runoff */
    double clay_node[MAX_SOILS];      /**< clay content of soil (fraction of mineral soil volume) */
    double depth[MAX_LAYERS];         /**< thickness of each soil moisture layer (m) */
    double dz_soil[MAX_SOILS];        /**< thermal node thickness (m) */
    double Zsum_soil[MAX_SOILS];      /**< thermal node depth (m) */
    double zc_soil[MAX_SOILS];        /**< depth of thermal nodes below soil surface (m) */
    double gravel_node[MAX_SOILS];    /**< gravel content of soil (fraction of mineral soil weight) */
    double organic_node[MAX_SOILS];   /**< organic content of soil (fraction of total soil volume) */
    double soil_dens_min[MAX_SOILS];  /**< particle density of mineral soil (kg/m^3) */
    double soil_dens_node[MAX_SOILS]; /**< soil particle density (kg/m^3) */
    double soil_dens_org[MAX_SOILS];  /**< particle density of organic soil [kg/m^3] */
    double sand_node[MAX_SOILS];      /**< sand content of soil (fraction of mineral soil volume) */
    double silt_node[MAX_SOILS];      /**< silt content of soil (fraction of mineral soil volume) */
    double lpar_node[MAX_SOILS];      /**< unsaturated hydraulic conductivity exponent in van Genuchten eqn. */
    double mpar_node[MAX_SOILS];      /**< unsaturated hydraulic conductivity exponent in van Genuchten eqn. */
    double *BandElev;                 /**< Elevation of each snow elevation band */
    double *AreaFract;                /**< Fraction of grid cell included in each snow elevation band */
    double *Pfactor;                  /**< Change in Precipitation due to elevation (fract) in each snow elevation band */
    double *Tfactor;                  /**< Change in temperature due to elevation (K) in each snow elevation band */

    double elevation;                 /**< grid cell elevation [m] */
    double lat;                       /**< grid cell central latitude */
    double lng;                       /**< grid cell central longitude */
    double cell_area;                 /**< Area of grid cell (m^2) */
    double time_zone_lng;             /**< central meridian of the time zone */
    unsigned int gridcel;             /**< grid cell number */
    double off_gmt;
    double slope;
    double z_bedrock;                 /**< Depth to bedrock [m] */
} soil_con_struct;

/******************************************************************************
 * @brief   This structure stores information about the vegetation coverage of
 *          the current grid cell.
 *****************************************************************************/
typedef struct {
    double fcanopy[MONTHS_PER_YEAR];    /**< fractional area covered by plant
                                                            canopy (fraction) */
    double LAI[MONTHS_PER_YEAR];        /**< leaf area index */
    double SAI[MONTHS_PER_YEAR];        /**< stem area index */
    double Cv;              /**< fraction of vegetation coverage */
    int veg_class;          /**< vegetation class id number */
    double root[MAX_SOILS];
    size_t vegetat_type_num; /**< number of vegetation types in the grid
                                cell */
    int BandIndex;
    bool IS_GLAC;
    size_t Nroot;
} veg_con_struct;

/******************************************************************************
 * @brief   This structure stores parameters for individual vegetation types.
 *****************************************************************************/
typedef struct {
    double fcanopy[MONTHS_PER_YEAR];    /**< fractional area covered by plant
                                                            canopy (fraction) */
    double LAI[MONTHS_PER_YEAR];        /**< leaf area index */
    double SAI[MONTHS_PER_YEAR];        /**< stem area index */
    double reflleaf[MAX_SWBANDS];
    double reflstem[MAX_SWBANDS];
    double transleaf[MAX_SWBANDS];
    double transstem[MAX_SWBANDS];
    double Canopy_Upper;                /**< top of canopy (m) */
    double Canopy_Lower;                /**< bottom of canopy (m) */
    double Canopy_Radius;
    double COI;
    double c_biomass;
    double d_leaf;
    double root_a;               /**< Empirical parameter a in eqa(2) */
    double root_b;               /**< Empirical parameter b in eqa(2) */
    double root_d;               /**< Maximum root depth (m) */
    unsigned short int veg_class; /**< vegetation class reference number */
    size_t NVegLibTypes;   /**< number of vegetation classes defined in
                              library */
    double trunk_dia;    /**< ratio of trunk height to tree height,
                              default = 0.2 (fraction) */
    double liq_bioms;
    double slatop;
    double stem_num;
    double Z0sub_LAImax;
    double Z0sub_Cs;
    double Z0sub_Cr;
    double Z0sub_c;
    double Z0sub_cw;
    double smpsc;
    double smpso;
    // Carbon terms
    char Ctype;                   /**< Photosynthetic pathway; 0 = C3; 1 = C4 */
    double froot_leaf;            /**< ratio of fine root mass to leaf mass */
    double matric50;              /**< matric potential at which stomatal conductance is reduced by 50% (m) */
    double kcano_max;             /**< plant segment max conductance. m h2o (transpired)/m h2o (water potential gradient)/sec [1/s] */
    double kroot_max;             /**< root segment max conductance. m h2o (transpired)/m h2o (water potential gradient)/sec [1/s] */
    double theta_cj;
    double leaf_CN;               /**< Leaf C:N [gC/gN] */
    double SLA_top;               /**< specific leaf area at top of canopy [m2/gC] */
    double fN_rub;                /**< fraction of leaf N in Rubisco enzyme (gN Rubisco/gN leaf) */ 
    double medlynslope;           /**< slope of Medlyn conductance-photosynthesis relationship */
    double medlynint;             /**< intercept of Medlyn conductance-photosynthesis relationship */
} veg_lib_struct;

/******************************************************************************
 * @brief   This structure stores vegetation parameter forcing data for each
 * model time step for a single veg tile.  Each array stores the values for the
 * SNOW_STEPs during the current model step and the value for the entire model
 * step.  The latter is referred to by array[NR].  Looping over the SNOW_STEPs
 * is done by for (i = 0; i < NF; i++)
 *****************************************************************************/
typedef struct {
    double *fcanopy;      /**< fractional area covered by plant canopy
                             (fraction) */
    double *LAI;          /**< leaf area index [m2/m2] */
    double *SAI;          /**< stem area index [m2/m2] */
} veg_hist_struct;

/******************************************************************************
 * @brief   This structure stores the forcing data for each model
 * time step for a single grid cell.  Each array stores the values for the
 * SNOW_STEPs during the current model step and the value for the entire model
 * step.  The latter is referred to by array[NR].  Looping over the SNOW_STEPs
 * is done by for (i = 0; i < NF; i++)
 *****************************************************************************/
typedef struct {
    double *air_temp;   /**< air temperature (K) */
    double *prec;       /**< average precipitation in grid cell (mm/s) */
    double *snowf;      /**< snowfall partitioned from precipitation (mm/s) */
    double *rainf;      /**< rainfall partitioned from precipitation (mm/s) */
    double *wind;       /**< wind speed (m/s) */
    double *Qair;       /**< specific humidity (kg/kg) */
    double *pressure;   /**< atmospheric pressure (kPa) */
    double *shortwave;  /**< incoming shortwave radiation (W/m^2) */
    double *longwave;   /**< incoming longwave radiation (W/m^2) (net incoming
                                           longwave for water balance model) */
    double *Catm;       /**< atmospheric CO2 mixing ratio (mol CO2/ mol air) */
    double *channel_in; /**< upstream input river channel inflow for time step (m/s) */
    double *coszen;     /**< cosine of the solar zenith angle */
    double *daylen;     /**< day length in seconds */
    double *density;    /**< atmospheric density (kg/m^3) */
    double *fdir;       /**< fraction of incoming shortwave that is direct (fraction) */
    double *par;        /**< incoming photosynthetically active radiation (μmol/m^2/s) */
    double *vp;         /**< atmospheric vapor pressure (kPa) */
    double *rel_humid;  /**< atmospheric relative humidity (%) */
    double *theta_pot;  /**< atmospheric potential temperature (K) */
    double *theta_v;    /**< atmospheric virtual potential temperature (K) */
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
    bool   IS_VEG;
    bool   IS_GLAC;
    size_t Nsoil;                       /**< Number of soil nodes in the model */
    size_t Nroot;                       /**< Number of root nodes in the model */
    size_t Nnode;                       /**< Number of thermal nodes in the model */
    size_t Ncanopy;                     /**< Number of canopy layer in model */
    double Ra_over[3];                  /**< The aerodynamic resistance (s/m) that was actually used
                                          in flux calculations. [0] = wind speed at reference height, 
                                          [1] = sensible heat flux, [2] = latent heat flux */
    double Ra_sub[3];                   /**< canopy ground surface resistance */
    double Ra_grnd[3];
    double Ra_root[MAX_SOILS];
    double Ra_evap;                    /**< ground surface resistance [s/m] to evaporation */
    double Ra_leaf;                    /**< canopy leaf resistance [s/m] to transpiration */
    double Ra_stem;
    double Z0m_grnd[3];
    double Z0m_sub[3];
    double displacement[2];            /**< displacement height (m) */
    double ref_height[3];              /**< height of meteorological forcing data (m) */
    double rh_grnd;                    /**< relative humidity on the ground (fraction) */
    double asat;                       /**< saturated area fraction */
    double rootmoist;                  /**< total of layer.moist over all layers in the root zone (mm) */
    double zwt;                        /**< average water table position [cm] - using lowest unsaturated layer */
    double max_daylen;                 /**< maximum daylength for this grid cell (s) */
    //double Qair_over;                  /**< specific humidity of the air at the canopy layer (kg/kg) */
    double Qair_grnd;                  /**< specific humidity of the air at the ground surface (kg/kg) */
    double Qair_soil;                  /**< specific humidity of the air at the soil surface (kg/kg) */
    double Qair_snow;                  /**< specific humidity of the air at the snow surface (kg/kg) */
    double ice[MAX_SOILS];             /**< ice content of the soil sublayer [m3/m3] */
    double liq[MAX_SOILS];             /**< liq content of the soil sublayer [m3/m3] */
    double last_ice[MAX_SOILS];
    double last_liq[MAX_SOILS];
    double moist[MAX_SOILS];           /**< moisture content of the unfrozen sublayer [m3/m3] */
    double porosity[MAX_SOILS];
    double soil_T[MAX_SOILS];
    double zc_node[MAX_NODES];     /**< depth of thermal nodes (m) */
    double dz_node[MAX_NODES];     /**< the thickness of each layer (m) */
    double Zsum_node[MAX_NODES];   /**< depth of bottom of each thermal node (m) */
    // Fluxes
    double baseflow;                   /**< baseflow from current cell (mm/TS) */
    double runoff;                     /**< runoff from current cell (mm/TS) */
    double soil_inflow;                /**< moisture that reaches the top of
                                            the soil column (mm) */
    double recharge;                   /**< aquifer recharge rate (mm/s) */
    double storage_aqf;
    double evap;
    double snowfrost;
    double snow_sublim;
    double lateral_flow[MAX_SOILS];
    double deriv_vapor[MAX_SOILS];
    // Canopy terms
    double transp;
    double canopyevap;
    double canopydew;
    double canopyfrost;
    double canopy_sublim;
    double canopy_vapor;
    // Soil terms
    size_t Nthaw;
    size_t Nfrost;
    double fdepth;
    double tdepth;
    double esoil;                      /**< soil evaporation from soil layer (mm) */
    double esoil_sub;
    double esoil_grnd;
    double dewsoil;                    /**< evapotranspiration from soil layer (mm) */
    double transp_fact;                   /**< soil water transpiration factor (0 to 1) */
    double soil_excess;
    double root[MAX_SOILS];
    double hksr_int[MAX_SOILS];        /**< soil-root interface conductance (mm/s) */
    double Netroot[MAX_SOILS];
    double liquid_flux[MAX_SOILS];
    double vapor_flux[MAX_SOILS];
    double soil_imped[MAX_SOILS];      /**< frost content of the frozen sublayer */
    double transp_sink[MAX_SOILS];     /**< transpiration sink term [m/s] */
    double conductivity[MAX_SOILS];    /**< soil hydraulic conductivity [m/s] */
    double conduct_int[MAX_SOILS];     /**< soil hydraulic conductivity for interface between layers */
    double matric[MAX_SOILS];          /**< soil matric potential [mm] */
    double last_matric[MAX_SOILS];
} cell_data_struct;

/******************************************************************************
 * @brief   This structure stores energy balance components, and variables used
 *          to solve the thermal fluxes through the soil column.
 *****************************************************************************/
typedef struct {
    // State variables
    bool FrozenGrnd;                   /**< TRUE = frozen soil present */
    bool FrozenOver;                   /**< TRUE = frozen canopy present */
    double kappa_node[MAX_NODES];      /**< thermal conductivity of the soil thermal nodes (W/m/K) */
    double Cs_node[MAX_NODES];         /**< heat capacity of the soil thermal nodes (J/m^3/K) */   
    double last_Cs[MAX_NODES]; 
    double T[MAX_NODES];               /**< thermal node temperatures (k) */
    double last_T[MAX_NODES];
    double kappa_int[MAX_NODES];       /**< thermal conductivity used for interface between nodes (W/m/K) */
    double fact[MAX_NODES];            /**< thermal diffusivity of each soil thermal node (m^2/s) */
    double Tcanopy;              /**< temperature of the canopy */
    double Tsurf;                /**< temperature of the understory */
    double Tgrnd;
    double Tfoliage;
    double Tstem;                /**< temperature of the stem */
    double deriv_grnd;           /**< terms in the energy balance that are linear with respect to the surface temperature (W/m^2/K) */
    double deriv_sub;            /**< terms in the energy balance that are linear with respect to the surface temperature (W/m^2/K) */
    double deriv_terms;          /**< sum of all terms in the energy balance that are linear with respect to the surface temperature (W/m^2/K) */
    double deriv_egrnd;
    double deriv_esub;
    double deriv_evap;
    double qsdT;                /**< temperature derivative of "Qair_grnd" */
    double delt_T;
    double delt_Q;              
    double error;                /**< energy balance error (W/m^2) */
    // Fluxes
    double advection;            /**< advective flux (Wm-2) */
    double AdvectSub;
    double AdvectGrnd;
    double AdvectOver;
    // 地表热通量
    double grnd_flux;            /**< ground heat flux (Wm-2) */
    // 感热通量
    double sensible;
    double SensibleGrnd;
    double SensibleSub;
    double SensibleStem;
    double SensibleLeaf;
    // 潜热通量
    double latent;
    double LatentGrnd;
    double LatentSub;
    double LatentLeaf;
    double LatentVapOver;
    double LatentVapGrnd;
    // 辐射通量
    double ReflShortSurf;
    double ReflShortGrnd;
    double ReflShortSub;
    double EmissLongSub;
    double EmissLongGrnd;
    double EmissLongSurf;
    // 长波辐射项
    double longwave;             /**< net longwave flux (Wm-2) */
    double NetLongSurf;
    double NetLongGrnd;          /**< net longwave radiation to the atmosphere (W/m^2) */
    double NetLongOver;          /**< net longwave radiation from the overstory (W/m^2) */
    double NetLongSub;         /**< net longwave radiation from the understory (W/m^2) */
    double NetLongOut;       /**< net longwave radiation from the ground and understory to the atmosphere (W/m^2) */
    // 短波辐射项
    double shortwave;
    double NetShortOver;
    double NetShortSurf;         /**< net shortwave to the atmosphere */
    double NetShortGrnd;         /**< net shortwave penetrating snowpack */
    double NetShortSub;        /**< net shortwave radiation from the understory (W/m^2) */
    double NetShortSoil;      /**< net shortwave radiation to the soil (W/m^2) */
    double NetShortSnow;      /**< net shortwave radiation to the snow (W/m^2) */
    // 辐射项
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
    double ReflectVeg[MAX_SWBANDS];
    double TransmitVeg[MAX_SWBANDS];
} energy_bal_struct;

/******************************************************************************
 * @brief   This structure stores vegetation variables for each vegetation type
 *          in a grid cell.
 *****************************************************************************/
typedef struct {
    // State variables
    double fcanopy;             /**< current fractional area of plant canopy (fraction) */
    double LAI;                 /**< current leaf area index (m2/m2) */
    double SAI;                 /**< current stem area index (m2/m2) */
    double NetLAI;              /**< net leaf area index (m2/m2) */
    double NetSAI;              /**< net stem area index (m2/m2) */
    double Wdew;                /**< maximum intercepted water per unit lai+sai (mm) */
    double MaxSnowInt;          /**< maximum canopy capacity for snow interception [mm] */
    double MaxRainInt;          /**< maximum canopy capacity for rain interception [mm] */
    double wetFrac;
    double dryFrac;
    double leaf_sun;
    double leaf_sha;
    double f_sun;
    double f_shade;
    double LAI_z[MAX_CANOPYS];  /**< leaf area index above the center of each canopy layer (m2/m2) */
    double SAI_z[MAX_CANOPYS];  /**< stem area index above the center of each canopy layer (m2/m2) */
    double mat_VEG[4];          /**< vegetation water matric potential (mm) [sun, shade, xylem, root] */
    // Fluxes
    double RainThroughFall;     /**< rain that reaches the ground through the canopy (mm/s) */
    double SnowThroughFall;     /**< snow that reaches the ground through the canopy (mm/s) */
    double SnowUnload;
    double RainDrip;
    double SnowDrip;
    double int_rain;            /**< rain intercepted on canopy (mm) */
    double int_snow;            /**< snow intercepted on canopy (mm) */
    double canopy_swq;          /**< snow water equivalent of the canopy (mm) */
    // PHS terms
    double aPAR_sun;            /**< par absorbed per unit lai for canopy layer (w/m**2) */
    double aPAR_sha;            /**< par absorbed per unit lai for canopy layer (w/m**2) */
    double ac_sun;              /**< Rubisco-limited gross photosynthesis (umol CO2/m**2/s) */
    double ac_sha;              /**< Rubisco-limited gross photosynthesis (umol CO2/m**2/s) */
    double ag_sun;              /**< co-limited gross leaf photosynthesis (umol CO2/m**2/s) */
    double ag_sha;              /**< co-limited gross leaf photosynthesis (umol CO2/m**2/s) */
    double aj_sun;              /**< RuBP-limited gross photosynthesis (umol CO2/m**2/s) */
    double aj_sha;              /**< RuBP-limited gross photosynthesis (umol CO2/m**2/s) */
    double ap_sun;              /**< product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s) */
    double ap_sha;              /**< product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s) */
    double an_sun;              /**< net sunlit leaf photosynthesis (umol CO2/m**2/s) */
    double an_sha;              /**< net shaded leaf photosynthesis (umol CO2/m**2/s) */
    double RS_sunlit;           /**< sunlit leaf stomatal resistance [s/m] */
    double RS_shade;            /**< shaded leaf stomatal resistance [s/m] */
    double ksun_vcmax;          /**< leaf to canopy scaling coefficient, sunlit leaf vcmax */
    double ksha_vcmax;          /**< leaf to canopy scaling coefficient, shaded leaf vcmax */
    double NetPhotosha;         /**< net shaded leaf photosynthesis (umol CO2/m**2/s) */
    double NetPhotosun;         /**< net sunlit leaf photosynthesis (umol CO2/m**2/s) */
    double PhotoError[2];          /**< photo synthetic error */
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
    double density[MAX_SNOWS];                 /**< snow density (kg/m^3) */
    double dz_snow[MAX_SNOWS];        /**< each snow pack depth (m) */
    double delta_depth;             /**< snow depth increasing rate [m/s] due to snowfall */
    double Zsum_snow[MAX_SNOWS];
    double zc_snow[MAX_SNOWS];      /**< depth of snow thermal nodes (m) */
    double snow_thresholds[MAX_SNOWS]; /**< snow depth thresholds for layer remobilization (m) */
    double glac_excess;
    double snow_depth;              /**< snow depth (m) */
    double snowage;                 /**< snow age (s) */
    double pack_T[MAX_SNOWS];       /**< temperature of each snow pack (k) */
    double pack_ice[MAX_SNOWS];     /**< ice content of the snow pack (mm) */
    double pack_liq[MAX_SNOWS];     /**< liquid water content of the snow pack (mm) */
    double pack_outflow[MAX_SNOWS]; /**< outflow of liq water from each snow pack (m/s) */
    double theta_ice[MAX_SNOWS];    
    double theta_liq[MAX_SNOWS];
    double porosity[MAX_SNOWS];     /**< porosity of each snow pack (fraction) */
    double snow_frac[MAX_SNOWS];    /**< fraction of each snow pack that is snow (fraction) */
    double last_snowfrac[MAX_SNOWS];
    double snow_outflow;            /**< outflow of liquid water from the snowpack bottom (m/s) */
    double new_snow_density;        /**< bulk density of snowfall [kg/m3] */
    double pack_melt;
    double pack_frze;
    double pack_transp;             /**< transpiration from each snow pack (m/s) */
    double pack_comb;               /**< combined heat and moisture of each snow pack (J/m^3) */
    double swq;                     /**< snow water equivalent of the entire pack (mm) */
    double last_swq;                 /**< snow water equivalent of the entire pack from previous time step (mm) */
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

#endif
