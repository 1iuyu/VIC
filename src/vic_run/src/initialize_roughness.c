/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes the roughness and forcing heights for each tiles.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    This routine initializes roughness lengths and 
 *           forcing heights [0]momentum, [1] sensible heat, [2] latent heat.
 *****************************************************************************/
void
initialize_roughness(double            Canopy_Upper,
                     double            coverage,
                     cell_data_struct *cell,
                     veg_var_struct   *veg_var)
{
    extern option_struct     options;
    extern parameters_struct param;
    double *Z0m_grnd = cell->Z0m_grnd;
    double *Z0m_sub = cell->Z0m_sub;
    double *displacement = cell->displacement;
    double beta = 1.0;
    double V, ratio_dis, ratio_Z0;
    double Net_VEG = veg_var->NetLAI + veg_var->NetSAI;

    if (options.AERO_RESIST == AR_ZENG) {
        if (coverage > 0.0) {
            Z0m_grnd[0] = param.SNOW_ROUGH;
        }
        else {
            Z0m_grnd[0] = param.SOIL_ROUGH;
        }
        Z0m_grnd[1] = Z0m_grnd[0];  // bare soil or glacier
        Z0m_grnd[2] = Z0m_grnd[0];  // veg cover
        if (cell->IS_VEG || cell->IS_URBAN) {
            if (Canopy_Upper > 0.0 && Canopy_Upper <= 1.0) {
                ratio_Z0 = param.VEG_RATIO_RL_A;
                ratio_dis = param.VEG_RATIO_DH_A;
            }
            else {
                ratio_Z0 = param.VEG_RATIO_RL_B;
                ratio_dis = param.VEG_RATIO_DH_B;
            }
            V = (1 - exp(-beta * min(Net_VEG, 2.0))) / (1 - exp(-beta * 2.0));
            Z0m_sub[0] = exp(V * log(Canopy_Upper * ratio_Z0) + (1 - V) * Z0m_grnd[0]);
            Z0m_sub[1] = Z0m_sub[0];
            Z0m_sub[2] = Z0m_sub[0];
            displacement[0] = 0.0;
            displacement[1] = calc_veg_displacement(V, ratio_dis, Canopy_Upper);
        }
        else {
            Z0m_sub[0] = 0.0;
            Z0m_sub[0] = 0.0;
            Z0m_sub[0] = 0.0;
            displacement[0] = 0.0;
            displacement[1] = 0.0;
        }
    }
    else if (options.AERO_RESIST == AR_MEIER) {
        if (coverage > 0.0) {
            Z0m_grnd[0] = param.SNOW_ROUGH;
        }
        else {
            Z0m_grnd[0] = param.SOIL_ROUGH;
        }
        Z0m_grnd[1] = Z0m_grnd[0];  // bare soil or glacier
        Z0m_grnd[2] = Z0m_grnd[0];  // veg cover
        if (cell->IS_VEG || cell->IS_URBAN) {
            if (Canopy_Upper > 0.0 && Canopy_Upper <= 1.0) {
                ratio_Z0 = param.VEG_RATIO_RL_A;
                ratio_dis = param.VEG_RATIO_DH_A;
            }
            else {
                ratio_Z0 = param.VEG_RATIO_RL_B;
                ratio_dis = param.VEG_RATIO_DH_B;
            }
            V = (1 - exp(-beta * min(Net_VEG, 2.0))) / (1 - exp(-beta * 2.0));
            Z0m_sub[0] = exp(V * log(Canopy_Upper * ratio_Z0) + (1 - V) * Z0m_grnd[0]);
            Z0m_sub[1] = Z0m_sub[0];
            Z0m_sub[2] = Z0m_sub[0];
            displacement[0] = 0.0; // 裸土
            displacement[1] = calc_veg_displacement(V, ratio_dis, Canopy_Upper); // 植被
        }
        else {
            Z0m_sub[0] = 0.0;
            Z0m_sub[0] = 0.0;
            Z0m_sub[0] = 0.0;
            displacement[0] = 0.0;
            displacement[1] = 0.0;
        }
    }
    else {
        log_err("Unknown AERO_RESIST option");
    }
}

