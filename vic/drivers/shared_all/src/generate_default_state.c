/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes the model state (energy balance, water balance, and
 * snow components) to default values.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    Initialize the model state (energy balance, water balance, and
 *           snow components) to default values.
 *****************************************************************************/
void
generate_default_state(all_vars_struct *all_vars,
                       soil_con_struct *soil_con,
                       veg_con_struct  *veg_con,
                       dmy_struct      *dmy_current)
{
    extern option_struct     options;
    extern parameters_struct param;

    size_t                   Nveg;
    size_t                   veg;
    size_t                   band;
    size_t                   lidx;
    size_t                   k;
    size_t                   tmpTshape[] = {
        options.Nlayer, options.Nnode,
        options.Nfrost + 1
    };
    size_t                   tmpZshape[] = {
        options.Nlayer, options.Nnode
    };
    double                   Cv;
    double                   tmp;
    double                   AreaFactor;
    double                   TreeAdjustFactor = 1.;
    double                   lakefactor = 1.;
    double                   albedo_sum;
    double                ***tmpT;
    double                 **tmpZ;
    int                      ErrorFlag;

    cell_data_struct        *cell;
    energy_bal_struct       *energy;

    cell = all_vars->cell;
    energy = all_vars->energy;
    Nveg = veg_con[0].vegetat_type_num;

    // allocate memory for tmpT and tmpZ
    malloc_3d_double(tmpTshape, &tmpT);
    malloc_2d_double(tmpZshape, &tmpZ);

    /************************************************************************
       Initialize soil moistures
       TBD: currently setting moist to init_moist from parameter file, but
            in future we should initialize to max_moist as default and
            eliminate init_moist (require user to use a state file if they
            want control over initial moist)
    ************************************************************************/

    for (veg = 0; veg <= Nveg; veg++) {
        Cv = veg_con[veg].Cv;
        if (Cv > 0) {
            band = veg_con[veg].BandIndex;
            if (soil_con->AreaFract[band] > 0.) {
                /* Initialize soil moistures */
                for (lidx = 0; lidx < options.Nlayer; lidx++) {
                    cell[veg].layer[lidx].moist =
                            soil_con->init_moist[lidx];
                    if (cell[veg].layer[lidx].moist >
                        soil_con->max_moist[lidx]) {
                        cell[veg].layer[lidx].moist =
                                soil_con->max_moist[lidx];
                    }
                }
            }
        }
    }

    /************************************************************************
       Initialize soil temperatures
    ************************************************************************/

    for (veg = 0; veg <= Nveg; veg++) {
        Cv = veg_con[veg].Cv;
        if (Cv > 0) {
            band = veg_con[veg].BandIndex;
            if (soil_con->AreaFract[band] > 0.) {
                /* Initialize soil node temperatures */
                for (k = 0; k < options.Nnode; k++) {
                    if (options.FULL_ENERGY || options.FROZEN_SOIL) {
                        energy[veg].T[k] = soil_con->avg_temp;
                    }
                    else {
                        energy[veg].T[k] = 0;
                    }
                }
                /* Initial estimate of LongUnderOut for use by snow_intercept() */
                tmp = energy[veg].T[0] + CONST_TKFRZ;
                energy[veg].LongUnderOut = calc_outgoing_longwave(tmp,
                                                                  param.EMISS_SNOW);
                energy[veg].Tfoliage = energy[veg].T[0] +
                                                 soil_con->Tfactor[band];
            }
        }
    }


    /************************************************************************
       Initialize gridcell-averaged albedo
    ************************************************************************/
    // vegetation class-weighted albedo over gridcell
    albedo_sum = 0;
    for (veg = 0; veg <= Nveg; veg++) {
        Cv = veg_con[veg].Cv;
        if (Cv > 0) {
            // TO-DO: account for treeline and lake factors
            AreaFactor = (Cv * TreeAdjustFactor * lakefactor);

            if (veg != Nveg) {
                // cold start, so using climatological albedo for all veg classes
                // except for bare soil
                albedo_sum += (AreaFactor *
                               veg_con[veg].albedo[dmy_current->month - 1]);
            }
            else {
                // bare soil class, use bare soil albedo
                albedo_sum += AreaFactor * param.ALBEDO_BARE_SOIL;
            }
        }
    }
    all_vars->gridcell_avg.avg_albedo = albedo_sum;

    /************************************************************************
       Initialize soil layer ice content
    ************************************************************************/

    for (veg = 0; veg <= Nveg; veg++) {
        Cv = veg_con[veg].Cv;
        if (Cv > 0) {
            band = veg_con[veg].BandIndex;
            if (soil_con->AreaFract[band] > 0.) {
                if (options.QUICK_FLUX) {
                    // TBD: calculation of layer ice content for quick flux
                    // depends on layer temperatures; so this initial
                    // estimation here is not ideal, since the layer
                    // temperature calculation is later in compute_derived_
                    // state_vars.
                    ErrorFlag =
                        estimate_layer_ice_content_quick_flux(
                                cell[veg].layer,
                                soil_con->depth,
                                soil_con->max_moist,
                                soil_con->expt, soil_con->bubble,
                                soil_con->frost_fract, soil_con->frost_slope,
                                soil_con->FS_ACTIVE);
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
                    ErrorFlag = estimate_layer_ice_content(
                            cell[veg].layer,
                            tmpT,
                            tmpZ,
                            soil_con->Zsum_node,
                            soil_con->depth,
                            soil_con->max_moist,
                            soil_con->expt,
                            soil_con->bubble,
                            options.Nnode,
                            options.Nlayer,
                            soil_con->FS_ACTIVE);
                    if (ErrorFlag == ERROR) {
                        log_err("Error calculating layer ice content");
                    }
                }
            }
        }
    }
    // free memory for tmpT and tmpZ
    free_3d_double(tmpTshape, tmpT);
    free_2d_double(tmpZshape, tmpZ);
}
