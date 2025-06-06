/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine computes the state variables (energy balance, water balance,
 * and snow components) that are derived from the variables that are stored in
 * state files.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    Compute the state variables (energy balance, water balance,
 *           and snow components) that are derived from the variables that
 *           are stored in state files.
 *****************************************************************************/
void
compute_derived_state_vars(all_vars_struct *all_vars,
                           soil_con_struct *soil_con,
                           veg_con_struct  *veg_con)
{
    extern global_param_struct global_param;
    extern option_struct       options;

    char                       FIRST_VEG;
    size_t                     Nveg;
    size_t                     veg;
    size_t                     lidx;
    size_t                     band;
    size_t                     tmpMshape[] = {
        veg_con[0].vegetat_type_num, options.Nlayer
    };
    size_t                     tmpTshape[] = {
        options.Nlayer, options.Nnode,
        options.Nfrost + 1
    };
    size_t                     tmpZshape[] = {
        options.Nlayer, options.Nnode
    };
    int                        ErrorFlag;
    double                     Cv;
    double                     dt_thresh;
    double                     tmp_runoff;
    double                   **tmpM;
    double                  ***tmpT;
    double                   **tmpZ;

    cell_data_struct          *cell;
    energy_bal_struct         *energy;
    snow_data_struct          *snow;

    cell = all_vars->cell;
    energy = all_vars->energy;
    snow = all_vars->snow;
    Nveg = veg_con[0].vegetat_type_num;

    // allocate memory for tmp* arrays
    malloc_2d_double(tmpMshape, &tmpM);
    if (!options.QUICK_FLUX) {
        malloc_3d_double(tmpTshape, &tmpT);
        malloc_2d_double(tmpZshape, &tmpZ);
    }

    /******************************************
       Compute derived soil layer vars
    ******************************************/
    for (veg = 0; veg <= Nveg; veg++) {
        // Initialize soil for existing vegetation types
        Cv = veg_con[veg].Cv;

        if (Cv > 0) {
            band = veg_con[veg].BandIndex;
            // Initialize soil for existing snow elevation bands
            if (soil_con->AreaFract[band] > 0.) {
                // set up temporary moist array
                for (lidx = 0; lidx < options.Nlayer; lidx++) {
                        tmpM[veg][lidx] =
                            cell[veg].layer[lidx].moist;
                }

                // compute saturated area and water table
                compute_runoff_and_asat(soil_con, tmpM[veg], 0,
                                        &(cell[veg].asat),
                                        &tmp_runoff);
                wrap_compute_zwt(soil_con, &(cell[veg]));
            }
        }
    }

    /******************************************
       Compute derived soil snow state vars
    ******************************************/
    for (veg = 0; veg <= Nveg; veg++) {
        band = veg_con[veg].BandIndex;
        if (snow[veg].density > 0.) {
            snow[veg].depth = CONST_RHOFW * snow[veg].swq /
                                    snow[veg].density;
        }
    }

    /******************************************
       Compute soil thermal node properties
    ******************************************/
    FIRST_VEG = true;
    for (veg = 0; veg <= Nveg; veg++) {
        // Initialize soil for existing vegetation types
        Cv = veg_con[veg].Cv;

        if (Cv > 0) {
            band = veg_con[veg].BandIndex;
            // Initialize soil for existing snow elevation bands
            if (soil_con->AreaFract[band] > 0.) {
                /** Set soil properties for all soil nodes **/
                if (FIRST_VEG) {
                    FIRST_VEG = false;
                    set_node_parameters(soil_con->Zsum_node,
                                        soil_con->porosity_node,
                                        soil_con->expt_node,
                                        soil_con->bubble_node,
                                        soil_con->alpha, soil_con->beta,
                                        soil_con->gamma, soil_con->depth,
                                        soil_con->porosity, soil_con->expt,
                                        soil_con->bubble,
                                        options.Nnode, options.Nlayer);
                }

                // set soil moisture properties for all soil thermal nodes
                if (options.FULL_ENERGY || options.FROZEN_SOIL) {
                    ErrorFlag =
                            distribute_node_moisture_properties(
                                energy[veg].moist,
                                energy[veg].ice,
                                energy[veg].kappa_node,
                                energy[veg].Cs_node,
                                soil_con->Zsum_node,
                                energy[veg].T,
                                soil_con->porosity_node,
                                soil_con->expt_node,
                                soil_con->bubble_node,
                                tmpM[veg],
                                soil_con->depth,
                                soil_con->soil_dens_min,
                                soil_con->bulk_dens_min,
                                soil_con->quartz,
                                soil_con->soil_density,
                                soil_con->bulk_density,
                                soil_con->organic,
                                options.Nnode, options.Nlayer,
                                soil_con->FS_ACTIVE);
                    if (ErrorFlag == ERROR) {
                        log_err("Error setting physical properties for "
                                "soil thermal nodes");
                    }
                }

                // Check node spacing v time step
                // (note this is only approximate since heat capacity and
                // conductivity can change considerably during the
                // simulation depending on soil moisture and ice content)
                if ((options.FROZEN_SOIL &&
                        !options.QUICK_FLUX) && !options.IMPLICIT) {
                        // in seconds
                        dt_thresh = 0.5 * energy[veg].Cs_node[1] /
                                    energy[veg].kappa_node[1] *
                                    pow((soil_con->dz_node[1]),
                                        2);
                        if (global_param.dt > dt_thresh) {
                            log_err("You are currently running FROZEN SOIL "
                                    "with an explicit method (IMPLICIT is "
                                    "set to FALSE).  For the explicit method "
                                    "to be stable, time step %f seconds is too "
                                    "large for the given thermal node spacing "
                                    "%f m, soil heat capacity %f J/m3/K, and "
                                    "soil thermal conductivity %f J/m/s/K.  "
                                    "Either set IMPLICIT to TRUE in your "
                                    "global parameter file (this is the "
                                    "recommended action), or decrease time "
                                    "step length to <= %f seconds, or decrease "
                                    "the number of soil thermal nodes.",
                                    global_param.dt,
                                    soil_con->dz_node[1],
                                    energy[veg].Cs_node[1],
                                    energy[veg].kappa_node[1], dt_thresh);
                        }
                }

                /* calculate soil layer temperatures  */
                if (options.QUICK_FLUX) {
                        ErrorFlag =
                            estimate_layer_temperature_quick_flux(
                                cell[veg].layer,
                                soil_con->depth, soil_con->dp,
                                energy[veg].T[0],
                                energy[veg].T[1],
                                soil_con->avg_temp);
                    if (ErrorFlag == ERROR) {
                        log_err("Error calculating layer temperature "
                                "using QUICK_FLUX option");
                    }
                }
                else {
                    estimate_frost_temperature_and_depth(
                            tmpT,
                            tmpZ,
                            soil_con->Zsum_node,
                            energy[veg].T,
                            soil_con->depth,
                            soil_con->frost_fract,
                            soil_con->frost_slope,
                            options.Nnode,
                            options.Nlayer);
                    ErrorFlag = estimate_layer_temperature(
                            cell[veg].layer,
                            tmpT,
                            tmpZ,
                            soil_con->Zsum_node,
                            soil_con->depth,
                            options.Nnode,
                            options.Nlayer);
                    if (ErrorFlag == ERROR) {
                            log_err("Error calculating layer temperature");
                    }
                }

                /* Find freezing and thawing front depths */
                if (!options.QUICK_FLUX && soil_con->FS_ACTIVE) {
                        find_0_degree_fronts(&energy[veg],
                                             soil_con->Zsum_node,
                                             energy[veg].T,
                                             options.Nnode);
                }
            }
        }
    }
    // free memory for tmp* arrays
    free_2d_double(tmpMshape, tmpM);
    if (!options.QUICK_FLUX) {
        free_3d_double(tmpTshape, tmpT);
        free_2d_double(tmpZshape, tmpZ);
    }
}
