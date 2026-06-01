/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine calculate the surface albedo and two-stream fluxes.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief  
 * Compute the surface albedo and two-stream fluxes.
 *****************************************************************************/
int
surface_albedo(double             step_dt,
               double             coszen,
               double             snowfall,
               energy_bal_struct *energy,
               cell_data_struct  *cell,
               snow_data_struct  *snow,
			   soil_con_struct   *soil_con,
               veg_var_struct    *veg_var,
               veg_lib_struct    *veg_lib)
{
    extern option_struct options;
    extern parameters_struct param;
    size_t i, Nswband;
    double f_sun = 0.0;
    double f_shade = 0.0;
    double tau_beam = 0.0;
    double leaf_frac = 0.0;
    double stem_frac = 0.0;
    double proj_area = 0.0;
    double f_snowage = 0.0;
    double coverage = snow->coverage;
    double NetLAI = veg_var->NetLAI;
    double NetSAI = veg_var->NetSAI;
    double wetFrac = veg_var->wetFrac;
    double *LAI_z = veg_var->LAI_z;
    double *SAI_z = veg_var->SAI_z;

    Nswband = options.Nswband;
    // initialize albedo and two-stream fluxes
    for (i = 0; i < Nswband; i++) {
        energy->AlbedoSnowDir[i] = 0.0;
        energy->AlbedoSnowDfs[i] = 0.0;
        energy->AlbedoSoilDir[i] = 0.0;
        energy->AlbedoSoilDfs[i] = 0.0;
        energy->AlbedoGrndDir[i] = 0.0;
        energy->AlbedoGrndDfs[i] = 0.0;
        energy->AbsSubDir[i] = 0.0;
        energy->AbsSubDfs[i] = 0.0;
        energy->ShortDir2Dir[i] = 0.0;
        energy->ShortDfs2Dir[i] = 0.0;
        energy->ShortDir2Dfs[i] = 0.0;
        energy->ShortDfs2Dfs[i] = 0.0;
        energy->AlbedoSurfDir[i] = 0.0;
        energy->AlbedoSurfDfs[i] = 0.0;
        energy->ReflGrndDir[i] = 0.0;
        energy->ReflGrndDfs[i] = 0.0;
        energy->ReflSubDir[i] = 0.0;
        energy->ReflSubDfs[i] = 0.0;
        energy->ReflectVeg[i] = 0.0;
        energy->TransmitVeg[i] = 0.0;
    }
    // compute snow age factor
    f_snowage = snow_aging(step_dt,
                           energy->Tsurf,
                           snowfall, snow);

    /** compute understory albedo and net shortwave radiation **/
    if (coszen > 0.) {

        leaf_frac = NetLAI / max(NetLAI + NetSAI, param.TOL_A);
        stem_frac = NetSAI / max(NetLAI + NetSAI, param.TOL_A);

        for (i = 0; i < Nswband; i++) {
            energy->ReflectVeg[i] = max(veg_lib->reflleaf[i] * leaf_frac + 
                                        veg_lib->reflstem[i] * 
                                            stem_frac, param.TOL_A);
            energy->TransmitVeg[i] = max(veg_lib->transleaf[i] * leaf_frac +
                                        veg_lib->transstem[i] * 
                                            stem_frac, param.TOL_A);
        }
        // age snow albedo if no new snowfall
        // solar radiation process is only done if there is light
        snow_albedo(coszen, f_snowage, energy);

        /* Compute ground albedo based on soil and snow albedo */
        GroundAlbedo(cell->moist[0],
                     coverage,
                     energy, soil_con);
        
        // Diagnose number of canopy layers for radiative transfer
        size_t nrad = 0;
        double frac_veg = 0.25;
        double veg_sum = 0.0;
        for (i = 0; i < MAX_CANOPYS; i++) {
            if (options.Ncanopy == 1) {
                nrad = 1;
                LAI_z[i] = NetLAI;
                SAI_z[i] = NetSAI;
            }
            else if (options.Ncanopy > 1) {
                if (NetLAI + NetSAI <= 0.0) {
                    nrad = 0;
                    LAI_z[i] = 0.0;
                    SAI_z[i] = 0.0;
                }
                else {
                    veg_sum += frac_veg;
                    if (NetLAI + NetSAI - veg_sum > param.TOL_A) {
                        nrad = i + 1;
                        LAI_z[i] = NetLAI * frac_veg / max(NetLAI + NetSAI, param.TOL_A);
                        SAI_z[i] = NetSAI * frac_veg / max(NetLAI + NetSAI, param.TOL_A);
                    }
                    else {
                        nrad = i + 1;
                        LAI_z[i] = max(NetLAI - LAI_z[i - 1], 0.0);
                        SAI_z[i] = max(NetSAI - SAI_z[i - 1], 0.0);
                        break;
                    }
                }
            }
        }
        // if the canopy is too sparse, set it to 4 layers with equal LAI and SAI
        if (nrad < 4) {
            options.Ncanopy = 4;
            for (i = 0; i < 4; i++) {
                LAI_z[i] = NetLAI / 4.0;
                SAI_z[i] = NetSAI / 4.0;
            }
        }
        double sum_LAI = 0.0;
        double sum_SAI = 0.0;
        for (i = 0; i < options.Ncanopy; i++) {
            sum_LAI += LAI_z[i];
            sum_SAI += SAI_z[i];
        }
        if (fabs(sum_LAI - NetLAI) > param.TOL_A || fabs(sum_SAI - NetSAI) > param.TOL_A) {
            log_err("Error in diagnosing canopy layers: sum of LAI_z or SAI_z does not equal NetLAI or NetSAI");
        }

        /* Compute canopy radiative transfer 
           using two-stream approximation */
        for (i = 0; i < Nswband; i++) {
            // direct
            canopy_two_stream(i, SUNLIT,
                              coszen, &proj_area,
                              energy, cell,
                              veg_var, veg_lib);
            // diffuse
            canopy_two_stream(i, SHADE,
                              coszen, &proj_area, 
                              energy, cell,
                              veg_var, veg_lib);
        }
        tau_beam = proj_area / coszen *
                        sqrt(1.0 - energy->ReflectVeg[0] -
                                            energy->TransmitVeg[0]);
        f_sun = (1.0 - exp(-tau_beam *
                            (NetLAI + NetSAI))) /
                                        max(tau_beam * 
                                                (NetLAI + NetSAI), 
                                                        param.TOL_A);
        tau_beam = f_sun;
        if (tau_beam < 0.01) {
            leaf_frac = 0.0;
        }
        else {
            leaf_frac = tau_beam;
        }
        f_sun = leaf_frac;
    } // end of coszen > 0.

    /* shaded canopy fraction */
    f_shade = 1.0 - f_sun;
    // update veg_var
    veg_var->f_sun = f_sun;
    veg_var->f_shade = f_shade;
    veg_var->leaf_sun = NetLAI * f_sun;
    veg_var->leaf_sha = NetLAI * f_shade;

    return (0);
}