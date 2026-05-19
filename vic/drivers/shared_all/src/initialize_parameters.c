/******************************************************************************
 * @section DESCRIPTION
 *
 * Initialize parameters structure.
 *****************************************************************************/

#include "vic_driver_shared_all.h"

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
    param.SVP_A0 = 6.11213476;
    param.SVP_A1 = 0.444007856;
    param.SVP_A2 = 0.143064234e-01;
    param.SVP_A3 = 0.264461437e-03;
    param.SVP_A4 = 0.305903558e-05;
    param.SVP_A5 = 0.196237241e-07;
    param.SVP_A6 = 0.892344772e-10;
    param.SVP_A7 = -0.373208410e-12;
    param.SVP_A8 = 0.209339997e-15;
    param.SVP_B0 = 6.11123516;
    param.SVP_B1 = 0.503109514;
    param.SVP_B2 = 0.188369801e-01;
    param.SVP_B3 = 0.420547422e-03;
    param.SVP_B4 = 0.614396778e-05;
    param.SVP_B5 = 0.602780717e-07;
    param.SVP_B6 = 0.387940929e-09;
    param.SVP_B7 = 0.149436277e-11;
    param.SVP_B8 = 0.262655803e-14;
    param.SVP_C0 = 0.444017302;
    param.SVP_C1 = 0.286064092e-01;
    param.SVP_C2 = 0.794683137e-03;
    param.SVP_C3 = 0.121211669e-04;
    param.SVP_C4 = 0.103354611e-06;
    param.SVP_C5 = 0.404125005e-09;
    param.SVP_C6 = -0.788037859e-12;
    param.SVP_C7 = -0.114596802e-13;
    param.SVP_C8 = 0.381294516e-16;
    param.SVP_D0 = 0.503277922;
    param.SVP_D1 = 0.377289173e-01;
    param.SVP_D2 = 0.126801703e-02;
    param.SVP_D3 = 0.249468427e-04;
    param.SVP_D4 = 0.313703411e-06;
    param.SVP_D5 = 0.257180651e-08;
    param.SVP_D6 = 0.133268878e-10;
    param.SVP_D7 = 0.394116744e-13;
    param.SVP_D8 = 0.498070196e-16;

    param.SVP_FRZ = 0.611;
    param.SVP_RDAIR = 0.622;
//    param.SVP1 = 17.67;
//    param.SVP2 = 29.65;
//    param.SVP3 = param.SVP1 * (CONST_TKFRZ - param.SVP2);

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
    param.SRESP_PSIWILT = -150; // [m]
    param.SRESP_EXP = 5.0;

    // Surface roughness constants
    param.ROUGH3 = 70.0;
    param.ROUGH_BETA = 7.2;
    param.ROUGH_NU = 1.5e-5;
    param.SNOW_ROUGH = 0.00085;
    param.SOIL_ROUGH = 0.000775;
    param.GLAC_ROUGH = 0.00230000005;

    // Snow Parameters
    param.SNOW_MAX_SURFACE_SWE = 5000.0;
    param.SNOW_LIQUID_WATER_CAPACITY = 0.03;   // 0.03
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
    param.SNOW_BETADS = 0.5;
    param.SNOW_BETAIS = 0.5;

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

    // two-stream parameter omega for snow
    param.SNOW_OMEGAS[0] = 0.8;
    param.SNOW_OMEGAS[1] = 0.4;

    // glacier parameters
    param.GLAC_ALBEDO[0] = 0.80;
    param.GLAC_ALBEDO[1] = 0.55;

    // MOST Iterations
    param.MAX_ITER_MOST = 5;
    param.MAX_ITER_OVER = 50;

    param.CN_FACTOR = 0.5;

}
