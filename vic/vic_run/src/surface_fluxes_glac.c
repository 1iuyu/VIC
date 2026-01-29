/******************************************************************************
* @section DESCRIPTION
*
* This routine computes all surface fluxes, and solves the snow accumulation
* and ablation algorithm. Solutions are for the current snow band and
* vegetation type.
******************************************************************************/

#include <vic_run.h>

/******************************************************************************
* @brief  This is a modified version of surface_fluxes.c specific for glaciers.
******************************************************************************/
int
surface_fluxes_glac(double               air_temp,
                    double               snowfall,
                    double               rainfall,
                    double               ref_height,
                    double              *roughness,
                    double              *displacement,
                    force_data_struct   *force,
                    energy_bal_struct   *energy,
                    global_param_struct *gp,
                    cell_data_struct    *cell,
                    snow_data_struct    *snow,
                    soil_con_struct     *soil_con)
{

    int                      ErrorFlag;
    int                      N_steps;
    size_t                   hidx; // index of initial element of atmos array
    size_t                   step_inc; // number of atmos array elements to skip per surface fluxes step
    size_t                   endhidx; // index of final element of atmos array
    double                   step_dt; // time length of surface fluxes step (in seconds)

    /*********************************
       Set-up sub-time step controls
    *********************************/
    hidx = 0;
    step_inc = 1;
    endhidx = hidx + NF;
    step_dt = gp->snow_dt;

    N_steps = 0;
    /*********************************
      Compute glacier surface fluxes
    *********************************/
    do {

        /** Solve energy balence processes **/
        ErrorFlag = calc_energy_bal_glac(hidx, step_dt,
                                         air_temp, rainfall,
                                         snowfall, ref_height,
                                         roughness,
                                         displacement,
                                         force, energy, cell,
                                         snow, soil_con);

        if (ErrorFlag == ERROR) {
            return (ERROR);
        }
        /** Solve water balence processes **/
        ErrorFlag = calc_water_bal_glac(hidx, air_temp, step_dt,
                                        snowfall, rainfall,
                                        energy->latent,
                                        energy->LatentVapGrnd,
                                        force, cell, snow, soil_con);
        if (ErrorFlag == ERROR) {

            return (ERROR);
        }

        /* increment time step */
        N_steps++;
        hidx += step_inc;

    }
    while (hidx < endhidx);

    return(0);
}
