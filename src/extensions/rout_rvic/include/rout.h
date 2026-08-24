/******************************************************************************
 * @section DESCRIPTION
 *
 * Header file for rvic routing routines
 *****************************************************************************/
#ifndef ROUT_RVIC_H
#define ROUT_RVIC_H

#define ROUT_EXT "rout_rvic"

#include "vic_def.h"
#include "vic_driver_shared_image.h"

#define DLevelH2R 5 /* number of sub-timesteps for hillslope routing */
#define DLevelR2C 3 /* number of sub-timesteps for channel routing */
#define TOL_VALUE 1.0e-14 /* tolerance value for convergence of routing */

/******************************************************************************
 * @brief   Routing Structs
 *****************************************************************************/
typedef struct {
    int downstream;                          /*1d array - downstream outlet index for each source grid cell */
    int indegree;                            /*1d array - number of upstream source grid cells for each source grid cell */
    int queue;                               /*1d array - queue of source grid cells to be processed */
    int routing_order;                       /*1d array - order of source grid cells to be processed */
    int source_torow;                        /*1d array - source y location*/
    int source_tocol;                        /*1d array - source x location*/
} rout_param_struct;

/******************************************************************************
 * @brief   main channel Struct
 *****************************************************************************/
typedef struct {
    double erin;                      /*scalar - channel inflow (mm/s)*/
    double erout;                     /*scalar - channel outflow (mm/s)*/
    double vr;                      /*scalar - channel velocity (m/s)*/
    double wr;                       /*scalar - channel storage (mm)*/
    double mr;
    double yr;
    double pr;
    double nr;                      /*scalar - Manning's n for the channel*/
    double rr;
    double erlateral;                       /*scalar - lateral inflow (mm/s)*/
    double rslpsqrt;                  /*scalar - channel slope (-)*/
    double dwr;                  /*scalar - change in channel storage (mm)*/
    double rlen;                          /*scalar - channel length (m)*/
    double rwidth;                          /*scalar - channel width (m)*/
    double rwidth0;                          /*scalar - channel width at bankfull (m)*/
    double rdepth;                          /*scalar - channel depth (m)*/
} main_channel_struct;

/******************************************************************************
 * @brief   subnetwork channel Struct
 *****************************************************************************/
typedef struct {
    double etin;                      /*scalar - channel inflow (mm/s)*/
    double etout;                     /*scalar - channel outflow (mm/s)*/
    double vt;                      /*scalar - channel velocity (m/s)*/
    double wt;                       /*scalar - channel storage (mm)*/
    double nt;                      /*scalar - Manning's n for the channel*/
    double mt;                       /*scalar - channel cross-sectional area (m^2)*/
    double tslpsqrt;                  /*scalar - channel slope (-)*/
    double dwt;                  /*scalar - change in channel storage (mm)*/
    double tlen;                          /*scalar - channel length (m)*/
    double twidth;                          /*scalar - channel width (m)*/
    double yt;                          /*scalar - channel depth (m)*/
    double pt;                          /*scalar - channel wetted perimeter (m)*/
    double rt;                          /*scalar - channel hydraulic radius (m)*/
} sub_channel_struct;

/******************************************************************************
 * @brief   hillslope Struct
 *****************************************************************************/
typedef struct {
    double ehout;
    double yh;
    double nh;                         /*scalar - Manning's n for the hillslope*/
    double wh;
    double hslpsqrt;
    double dwh;
    double hlen;                          /*scalar - channel length (m)*/
} hillslope_struct;

/******************************************************************************
 * @brief   main routing Struct
 *****************************************************************************/
typedef struct {
    rout_param_struct rout_param;
    main_channel_struct main_channel;
    hillslope_struct hillslope;
    sub_channel_struct sub_channel;
    double discharge;
    double runoff;
    double baseflow;
    double total_length;                             /*scalar - glacier runoff (mm/s)*/  
    double upstream;
    double drainage_density;                             /* drainage density within the cell, [1/m] */
    double acc_area;
    double total_storage_prev;
    size_t river_steps;
    size_t sub_steps;
} rout_struct;

/******************************************************************************
 * @brief   Function prototypes for the rout_rvic extension
 *****************************************************************************/
void rout_alloc(void);                 // allocate memory
void rout_init(void);                  // initialize model parameters from parameter files
void rout_run(void);                   // run routing over the domain

/******************************************************************************
 * @brief   MPI Function prototypes for the rout_rvic extension
 *****************************************************************************/
void scatter_var_double(double *, double *);
void gather_var_double(double *, double *);

/******************************************************************************
 * @brief   Convolution function adapted from the RVIC scheme
 *****************************************************************************/
void get_global_param_rout(FILE *gp);
int Euler_Routing(size_t, double);
void SearchCatchment(int *direction);
double GRMR(double, double);
double GRPT(double, double);
double GRHT(double, double);
double GRHR(double, double, double, double);
double GRPR(double, double, double, double);
double GRRR(double, double);
double CRVRMAN_nosqrt(double, double, double);
double CREHT_nosqrt(double, double, double, double);

#endif
