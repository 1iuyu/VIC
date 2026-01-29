/******************************************************************************
 * @section DESCRIPTION
 *
 * Initialize parameters structure.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    Initialize parameters structure.
 *****************************************************************************/
void
initialize_parameters()
{
    extern parameters_struct param;
    // Initialize temporary parameters

    // Lapse Rate
    param.LAPSE_RATE = -0.0065;

    // Precipitation Guage Height
    param.GAUGE_HEIGHT = 1.0;

    // Huge Resistance Term
    param.HUGE_RESIST = 1e20;

    // Surface Albedo Parameters
    param.ALBEDO_BARE_SOIL = 0.2;

    param.REF_HEIGHT = 2.0;
    param.REF_HEIGHT_WIND = 10.0;
    param.TOPOINDEX = 10.5;
    // Surface Emissivities
    param.EMISS_GRND = 0.97;
    param.EMISS_VEG = 0.97;
    param.EMISS_ICE = 0.98;
    param.EMISS_SNOW = 0.95;
    param.EMISS_H2O = 0.98;

    // Soil Constraints
    param.SOIL_RARC = 100.0;
    param.SOIL_RESID_MOIST = 0.0;
    param.SOIL_SLAB_MOIST_FRACT = 1.0;
    param.SOIL_WINDH = 10.0;
    param.SOIL_RESIST_EXP = 5.0;
    param.SOIL_FROST = 4.0;

    // Vegetation Parameters
    param.VEG_LAI_SNOW_MULTIPLIER = 0.0005;
    param.VEG_LAI_WATER_FACTOR = 0.1;
    param.VEG_MIN_INTERCEPTION_STORAGE = 0.005;
    param.VEG_RATIO_DH_HEIGHT = 0.67;
    param.VEG_RATIO_RL_HEIGHT = 0.123;

    // Canopy Parameters
    param.CANOPY_CLOSURE = 4000.0;
    param.CANOPY_RSMAX = 5000.0;
    param.CANOPY_VPDMINFACTOR = 0.1;

    // Saturation Vapor Pressure Parameters
    param.SVP_A0 = 6.107799961;
    param.SVP_A1 = 4.436518521e-01;
    param.SVP_A2 = 1.428945805e-02;
    param.SVP_A3 = 2.650648471e-04;
    param.SVP_A4 = 3.031240396e-06;
    param.SVP_A5 = 2.034080948e-08;
    param.SVP_A6 = 6.136820929e-11;
    param.SVP_B0 = 6.109177956;
    param.SVP_B1 = 5.034698970e-01;
    param.SVP_B2 = 1.886013408e-02;
    param.SVP_B3 = 4.176223716e-04;
    param.SVP_B4 = 5.824720280e-06;
    param.SVP_B5 = 4.838803174e-08;
    param.SVP_B6 = 1.838826904e-10;
    param.SVP_C0 = 4.438099984e-01;
    param.SVP_C1 = 2.857002636e-02;
    param.SVP_C2 = 7.938054040e-04;
    param.SVP_C3 = 1.215215065e-05;
    param.SVP_C4 = 1.036561403e-07;
    param.SVP_C5 = 3.532421810e-10;
    param.SVP_C6 = -7.090244804e-13;
    param.SVP_D0 = 5.030305237e-01;
    param.SVP_D1 = 3.773255020e-02;
    param.SVP_D2 = 1.267995369e-03;
    param.SVP_D3 = 2.477563108e-05;
    param.SVP_D4 = 3.005693132e-07;
    param.SVP_D5 = 2.158542548e-09;
    param.SVP_D6 = 7.131097725e-12;

    param.SVP_FRZ = 0.611;
    param.SVP_RDAIR = 0.622;
    param.SVP1 = 17.67;
    param.SVP2 = 29.65;
    param.SVP3 = param.SVP1 * (CONST_TKFRZ - param.SVP2);

    // Photosynthesis Parameters
    param.PHOTO_OMEGA = 0.12;
    param.PHOTO_LAIMAX = 8.0;
    param.PHOTO_LAILIMIT = 3.0;
    param.PHOTO_LAIMIN = 1.0e-9;
    param.PHOTO_EPAR = 2.2e5;
    param.PHOTO_FCMAX = 0.9;
    param.PHOTO_FCMIN = 1.0e-3;
    param.PHOTO_ZENITHMIN = 0.0174524;
    param.PHOTO_ZENITHMINPAR = 1.0e-3;
    param.PHOTO_ALBSOIPARMIN = 0.0;
    param.PHOTO_MINMAXETRANS = 1.0e-12;
    param.PHOTO_MINSTOMCOND = 0.0;
    param.PHOTO_FCI1C3 = 0.87;
    param.PHOTO_FCI1C4 = 0.67;
    param.PHOTO_OX = 0.21;
    param.PHOTO_KC = 460.0e-6;
    param.PHOTO_KO = 330.0e-3;
    param.PHOTO_EC = 59356.0;
    param.PHOTO_EO = 35948.0;
    param.PHOTO_EV = 58520.0;
    param.PHOTO_ER = 45000.0;
    param.PHOTO_ALC3 = 0.28;
    param.PHOTO_FRDC3 = 0.011;
    param.PHOTO_EK = 50967.0;
    param.PHOTO_ALC4 = 0.04;
    param.PHOTO_FRDC4 = 0.042;
    param.PHOTO_THETA = 0.83;
    param.PHOTO_FRLEAF = 0.4;
    param.PHOTO_FRGROWTH = 0.25;

    // Soil Respiration Parameters
    param.SRESP_E0_LT = 308.56;
    param.SRESP_T0_LT = 227.13;
    param.SRESP_WMINFM = 0.0;
    param.SRESP_WMAXFM = 1.0;
    param.SRESP_WOPTFM = 0.5;
    param.SRESP_RHSAT = 0.15;
    param.SRESP_RFACTOR = 0.5;
    param.SRESP_TAULITTER = 2.86;
    param.SRESP_TAUINTER = 33.3;
    param.SRESP_TAUSLOW = 1000.0;
    param.SRESP_FAIR = 0.7;
    param.SRESP_FINTER = 0.985;
    param.SRESP_PSIWILT = -150; 
    param.SRESP_EXP = 5.0;
    // Snow Parameters
    param.SNOW_MAX_SURFACE_SWE = 5000.0;
    param.SNOW_LIQUID_WATER_CAPACITY = 0.035;   // 0.03
    param.SNOW_MAX_LIQUID_FRAC = 0.4;
    param.SNOW_RELEASE_FAC = 5.0e-5;      // snowpack water release timescale factor (1/s)
    param.SNOW_NEW_SNOW_DENSITY = 50.0;
    param.SNOW_NEW_SNOW_DENS_MAX = 120.0;
    param.SNOW_DEPTH_THRES = 1.e-8;
    param.SNOW_DENS_DMLIMIT = 100.0;
    param.SNOW_DENS_DMLIMIT_FACTOR = 1.15;
    param.SNOW_DENS_MAX_CHANGE = 0.9;
    param.SNOW_DENS_ETA0 = 3.6e6;
    param.SNOW_DENS_C1 = 0.04;
    param.SNOW_DENS_C2 = 2.778e-6;
    param.SNOW_DENS_C3 = 1.0;
    param.SNOW_DENS_C3_CONST = -0.046;
    param.SNOW_DENS_C4 = 1.0;
    param.SNOW_DENS_C4WET = 2.0;
    param.SNOW_DENS_C5 = 0.08;
    param.SNOW_DENS_C6 = 0.021;
    param.SNOW_DENS_F = 0.6;
    param.SNOW_DENS_EXP = 0.35;
    param.SNOW_DENS_DENOM = 10.;
    param.SNOW_NEW_SNT_C1 = 67.92;
    param.SNOW_NEW_SNT_C2 = 51.25;
    param.SNOW_NEW_SNT_C3 = 2.59;
    param.SNOW_NEW_BRAS_DENOM = 100.;
    param.SNOW_MIN_SWQ_EB_THRES = 0.0010;
    param.SNOW_PGRAD = 0.0005;
    param.SNOW_RASNOW = 50.0;
    param.SNOW_A1 = 0.7;
    param.SNOW_A2 = 0.3;
    param.SNOW_L1 = 6.0;
    param.SNOW_L2 = 20.0;
    param.SNOW_TRACESNOW = 0.03;
    param.SNOW_CONDUCT = 2.9302e-6;
    param.SNOW_NEW_SNOW_ALB = 0.85;
    param.SNOW_ALB_ACCUM_A = 0.94;
    param.SNOW_ALB_ACCUM_B = 0.58;
    param.SNOW_ALB_THAW_A = 0.82;
    param.SNOW_ALB_THAW_B = 0.46;
    param.SNOW_COMPACT_A = 2.5e-6;
    param.SNOW_COMPACT_B = 0.04;
    param.SNOW_COMPACT_C = 2.0;
    param.SNOW_COMPACT_P = 21.0e-3;
    param.SNOW_COMPACT_DM = 100.0;
    param.SNOW_COMPACT_ETA = 1.33e+6;
    param.SNOW_RADIUS_MIN = 54.526;
    param.SNOW_RADIUS_MAX = 1500.;
    param.SNOW_NEW_RADIUS = 204.526;
    // BATS snow aging parameters
    param.SNOW_AGE_FACT = 1.0e6;
    param.SNOW_AGE_VAPF = 5000.;
    param.SNOW_AGE_FRZF = 10.;
    param.SNOW_AGE_SOTF = 0.3;
    param.SNOW_NEW_SNOW_COVER = 1.0;
    // BATS snow albedo parameters
    param.SNOW_COSZEN_B = 2.0;
    param.SNOW_AGE_DIR_VIS = 0.4;
    param.SNOW_AGE_DIR_NIR = 0.4;
    param.SNOW_AGE_DFS_VIS = 0.2;
    param.SNOW_AGE_DFS_NIR = 0.5;
    param.SNOW_NEW_SNOW_VIS = 0.95;
    param.SNOW_NEW_SNOW_NIR = 0.65;

    // Blowing Snow Parameters
    param.BLOWING_KA = 0.0245187;
    param.BLOWING_CSALT = 0.68;
    param.BLOWING_UTHRESH = 0.25;
    param.BLOWING_KIN_VIS = 1.3e-5;
    param.BLOWING_MAX_ITER = 100;
    param.BLOWING_K = 5;
    param.BLOWING_SETTLING = 0.3;
    param.BLOWING_NUMINCS = 10;

    // Treeline temperature
    param.TREELINE_TEMPERATURE = 10.0;

    // Iteration bracket widths
    param.SNOW_DT = 5.0;
    param.SURF_DT = 1.0;
    param.SOIL_DT = 0.25;
    param.CANOPY_DT = 1.0;
    param.CANOPY_VP = 25.0;

    // Solar radiation fraction
    param.RAD_DIR_F = 0.7;
    param.RAD_VIS_F = 0.5;

    // Convergence Tolerances
    param.TOL_A = 1.0e-6;
    param.TOL_B = 1.0e-8;
    param.MAX_LIMIT = 1.0e6;

    // Wetbulb Iterations
    param.TOL_WETBULB = 0.001;
    param.MAX_ITER_WETBULB = 10;

    // MOST Iterations
    param.MAX_ITER_MOST = 5;
    param.MAX_ITER_OVER = 50;

    // Frozen Soil Parameters
    param.FROZEN_MAXITER = 1000;

    // Canopy Iterations
    // initialized to 10, set to 0 if
    // options.CLOSE_ENERGY is false
    // this allows for flexibility in
    // changing the maximum number of
    // iterations
    param.MAX_ITER_GRND_CANOPY = 10;
    param.MAX_ITER_INFLI_RUNOFF = 3;

    // Newton-Raphson solver parameters
    param.NEWT_RAPH_MAXTRIAL = 150;
    param.NEWT_RAPH_TOLX = 1.0e-4;
    param.NEWT_RAPH_TOLF = 1.0e-1;
    param.NEWT_RAPH_R_MAX = 2.0;
    param.NEWT_RAPH_R_MIN = -5.0;
    param.NEWT_RAPH_RELAX1 = 0.9;
    param.NEWT_RAPH_RELAX2 = 0.7;
    param.NEWT_RAPH_RELAX3 = 0.2;
    param.NEWT_RAPH_EPS2 = 1.0e-4;

    // Root-Brent parameters
    param.ROOT_BRENT_MAXTRIES = 5;
    param.ROOT_BRENT_MAXITER = 1000;
    param.ROOT_BRENT_TSTEP = 10;
    param.ROOT_BRENT_T = 1.0e-7;
}
