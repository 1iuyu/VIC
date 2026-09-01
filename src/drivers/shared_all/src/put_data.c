/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine converts data units, and stores finalized values in an array
 * for later output to the output files.
 *****************************************************************************/

#include "vic_driver_shared_all.h"

/******************************************************************************
 * @brief    This routine converts data units, and stores finalized values in
 *           an array for later output to the output files.
 *****************************************************************************/
void
put_data(size_t             nveg,
         all_vars_struct   *all_vars,
         force_data_struct *force,
         veg_con_struct    *veg_con,
         double          ***out_data,
         timer_struct      *timer,
         stream_struct     *streams)
{
    extern global_param_struct global_param;
    extern option_struct       options;

    size_t veg;
    size_t i, j;
    bool HasVeg = false;
    bool HasGlac = false;
    bool output_var[N_OUTVAR_TYPES] = {false};

    cell_data_struct  *cell;
    energy_bal_struct *energy;
    snow_data_struct  *snow;
    veg_var_struct    *veg_var;

    cell = all_vars->cell;
    energy = all_vars->energy;
    snow = all_vars->snow;
    veg_var = all_vars->veg_var;

    double dt_sec = global_param.step_dt;

    for (j = 0; j < streams->nvars; j++) {
        size_t varid = streams->varid[j];
        output_var[varid] = true;
    }

    // Initialize output data to zero
    zero_output_list(nveg, out_data, streams);

    // Set output versions of input forcings
    if (output_var[OUT_AIR_TEMP]) {
        out_data[OUT_AIR_TEMP][0][0] = force->air_temp[NR];
    }
    if (output_var[OUT_DENSITY]) {
        out_data[OUT_DENSITY][0][0] = force->density[NR];
    }
    if (output_var[OUT_LWDOWN]) {
        out_data[OUT_LWDOWN][0][0] = force->longwave[NR];
    }
    if (output_var[OUT_PREC]) {
        out_data[OUT_PREC][0][0] = force->prec[NR] * dt_sec;  // mm over grid cell
    }
    if (output_var[OUT_PRESSURE]) {
        out_data[OUT_PRESSURE][0][0] = force->pressure[NR] / PA_PER_KPA;
    }
    if (output_var[OUT_QAIR]) {
        out_data[OUT_QAIR][0][0] = force->Qair[NR];
    }
    if (output_var[OUT_RAINF]) {
        out_data[OUT_RAINF][0][0] = force->rainf[NR] * dt_sec;   // mm over grid cell
    }
    if (output_var[OUT_REL_HUMID]) {
        out_data[OUT_REL_HUMID][0][0] = force->rel_humid[NR];
    }
    if (output_var[OUT_SWDOWN]) {
        out_data[OUT_SWDOWN][0][0] = force->shortwave[NR];
    }
    if (output_var[OUT_SNOWF]) {
        out_data[OUT_SNOWF][0][0] = force->snowf[NR] * dt_sec;   // mm over grid cell
    }
    if (output_var[OUT_VP]) {
        out_data[OUT_VP][0][0] = force->vp[NR] / PA_PER_KPA;
    }
    if (output_var[OUT_WIND]) {
        out_data[OUT_WIND][0][0] = force->wind[NR];
    }
    if (output_var[OUT_COSZEN]) {
        out_data[OUT_COSZEN][0][0] = force->coszen[NR];
    }
    if (options.CARBON) {
        if (output_var[OUT_CATM]) {
            out_data[OUT_CATM][0][0] = force->Catm[NR] / PPM_to_MIXRATIO;
        }
        if (output_var[OUT_FDIR]) {
            out_data[OUT_FDIR][0][0] = force->fdir[NR];
        }
        if (output_var[OUT_PAR]) {
            out_data[OUT_PAR][0][0] = force->par[NR];
        }
    }
    else {
        if (output_var[OUT_CATM]) {
            out_data[OUT_CATM][0][0] = MISSING;
        }
        if (output_var[OUT_FDIR]) {
            out_data[OUT_FDIR][0][0] = MISSING;
        }
        if (output_var[OUT_PAR]) {
            out_data[OUT_PAR][0][0] = MISSING;
        }
    }

    /****************************************
       Store Output for all Vegetation Types
    ****************************************/
    for (veg = 0; veg <= veg_con[0].vegetat_type_num; veg++) {
        // Determine whether the land cover is vegetation
        if (cell[veg].IS_VEG) {
            HasVeg = true;
        }
        else if (cell[veg].IS_GLAC) {
            HasGlac = true;
        }
        if (veg_con[veg].Cv > 0.0) {
            // Ground (soil/snow) surface
            double tmp_evap = 0.0;
            tmp_evap += cell[veg].evap;
            if (output_var[OUT_NET_EVAP]) {
                out_data[OUT_NET_EVAP][veg][0] = cell[veg].evap * dt_sec;
            }
            // Canopy evaporation
            if (HasVeg) {
                tmp_evap += cell[veg].canopyevap;
                tmp_evap += cell[veg].transp;
                if (output_var[OUT_EVAP_CANOP]) {
                    out_data[OUT_EVAP_CANOP][veg][0] = cell[veg].canopyevap * dt_sec;
                }
                if (output_var[OUT_TRANSP_VEG]) {
                    out_data[OUT_TRANSP_VEG][veg][0] = cell[veg].transp * dt_sec;
                }
                if (output_var[OUT_DEW_CANOP]) {
                    out_data[OUT_DEW_CANOP][veg][0] = cell[veg].canopydew * dt_sec;
                }
                if (output_var[OUT_FROST_CANOP]) {
                    out_data[OUT_FROST_CANOP][veg][0] = cell[veg].canopyfrost * dt_sec;
                }
                if (output_var[OUT_SUB_CANOP]) {
                    out_data[OUT_SUB_CANOP][veg][0] = cell[veg].canopy_sublim * dt_sec;
                }
                if (output_var[OUT_VAPOR_CANOP]) {
                    out_data[OUT_VAPOR_CANOP][veg][0] = cell[veg].canopy_vapor * dt_sec;
                }
                if (output_var[OUT_WDEW]) {
                    out_data[OUT_WDEW][veg][0] = veg_var[veg].Wdew;
                }
                if (output_var[OUT_INT_SNOW]) {
                    out_data[OUT_INT_SNOW][veg][0] = veg_var[veg].int_snow;
                }
                if (output_var[OUT_INT_RAIN]) {
                    out_data[OUT_INT_RAIN][veg][0] = veg_var[veg].int_rain;
                }
                if (output_var[OUT_CANOPY_SWQ]) {
                    out_data[OUT_CANOPY_SWQ][veg][0] = veg_var[veg].canopy_swq;
                }
                if (output_var[OUT_SNOW_DRIP]) {
                    out_data[OUT_SNOW_DRIP][veg][0] = veg_var[veg].SnowDrip;
                }
                if (output_var[OUT_RAIN_DRIP]) {
                    out_data[OUT_RAIN_DRIP][veg][0] = veg_var[veg].RainDrip;
                }
                if (output_var[OUT_SNOWTHROUGHFALL]) {
                    out_data[OUT_SNOWTHROUGHFALL][veg][0] = veg_var[veg].SnowThroughFall;
                }
                if (output_var[OUT_RAINTHROUGHFALL]) {
                    out_data[OUT_RAINTHROUGHFALL][veg][0] += veg_var[veg].RainThroughFall;
                }
            }

            // net evaporation = evap + canopyevap + transp
            if (output_var[OUT_EVAP]) {
                out_data[OUT_EVAP][veg][0] = tmp_evap * dt_sec;  // mm over gridcell
            }

            /** record saturated area fraction **/
            if (output_var[OUT_ASAT]) {
                out_data[OUT_ASAT][veg][0] = cell[veg].asat;
            }

            /** record runoff **/
            if (output_var[OUT_RUNOFF]) {
                out_data[OUT_RUNOFF][veg][0] = cell[veg].runoff;
            }

            /** record baseflow **/
            if (output_var[OUT_BASEFLOW]) {
                out_data[OUT_BASEFLOW][veg][0] = cell[veg].baseflow;
            }

            /** record recharge to groundwater storage */
            if (output_var[OUT_RECHARGE]) {
                out_data[OUT_RECHARGE][veg][0] = cell[veg].recharge;
            }

            /** record soil surface dew rate [mm/s] */
            if (output_var[OUT_SOIL_DEW]) {
                out_data[OUT_SOIL_DEW][veg][0] = cell[veg].soil_dew * dt_sec * MM_PER_M;
            }

            /** record soil surface frost rate [mm/s] */
            if (output_var[OUT_SOIL_FROST]) {
                out_data[OUT_SOIL_FROST][veg][0] = cell[veg].soil_frost * dt_sec * MM_PER_M;
            }

            /** record soil evaporation from soil layer [mm/s] */
            if (output_var[OUT_SOIL_EVAP]) {
                out_data[OUT_SOIL_EVAP][veg][0] = cell[veg].soil_evap * dt_sec * MM_PER_M;
            }

            /** record soil surface sublimation rate [mm/s] */
            if (output_var[OUT_SOIL_SUBLIM]) {
                out_data[OUT_SOIL_SUBLIM][veg][0] += cell[veg].soil_sublim * dt_sec * MM_PER_M;
            }

            /** record soil_inflow[mm] **/
            if (output_var[OUT_INFLOW]) {
                out_data[OUT_INFLOW][veg][0] = cell[veg].soil_inflow * dt_sec * MM_PER_M;
            }

            /** record LAI **/
            if (output_var[OUT_LAI]) {
                out_data[OUT_LAI][veg][0] = veg_var[veg].LAI;
            }

            /** record fcanopy **/
            if (output_var[OUT_FCANOPY]) {
                out_data[OUT_FCANOPY][veg][0] = veg_var[veg].fcanopy;
            }

            /** record aerodynamic conductance and resistance **/
            if (output_var[OUT_RA_OVER]) {
                for (i = 0; i < 3; i++) {
                    out_data[OUT_RA_OVER][veg][i] = cell[veg].Ra_over[i];
                }
            }
            if (output_var[OUT_RA_GRND]) {
                for (i = 0; i < 3; i++) {
                    out_data[OUT_RA_GRND][veg][i] = cell[veg].Ra_grnd[i];
                }
            }
            if (output_var[OUT_RA_SUB]) {
                for (i = 0; i < 3; i++) {
                    out_data[OUT_RA_SUB][veg][i] = cell[veg].Ra_sub[i];
                }
            }
            if (output_var[OUT_RA_EVAP]) {
                out_data[OUT_RA_EVAP][veg][0] = cell[veg].Ra_evap;
            }
            if (HasVeg) {
                if (output_var[OUT_RA_LEAF]) {
                    out_data[OUT_RA_LEAF][veg][0] = cell[veg].Ra_leaf;
                }
                if (output_var[OUT_RA_STEM]) {
                    out_data[OUT_RA_STEM][veg][0] = cell[veg].Ra_stem;
                }
            }

            /** record nodes moistures **/
            if (output_var[OUT_SOIL_LIQ]) {
                for (i = 0; i < cell[veg].Nsoil; i++) {
                    out_data[OUT_SOIL_LIQ][veg][i] = cell[veg].liq[i];
                }
            }
            if (output_var[OUT_SOIL_ICE]) {
                for (i = 0; i < cell[veg].Nsoil; i++) {
                    out_data[OUT_SOIL_ICE][veg][i] = cell[veg].ice[i];
                }
            }
            if (output_var[OUT_SOIL_MOIST]) {
                for (i = 0; i < cell[veg].Nsoil; i++) {
                    out_data[OUT_SOIL_MOIST][veg][i] = cell[veg].moist[i];
                }
            }
            if (output_var[OUT_MATRIC]) {
                for (i = 0; i < cell[veg].Nsoil; i++) {
                    out_data[OUT_MATRIC][veg][i] = cell[veg].matric[i];
                }
            }
            if (output_var[OUT_ROOTMOIST]) {
                out_data[OUT_ROOTMOIST][veg][0] = cell[veg].rootmoist;
            }

            /** record water table position **/
            if (output_var[OUT_ZWT]) {
                out_data[OUT_ZWT][veg][0] = cell[veg].zwt;
            }

            /*****************************
             Record Snow Pack Variables
            *****************************/
            /** record snow water equivalence **/
            if (output_var[OUT_SWE]) {
                out_data[OUT_SWE][veg][0] = snow[veg].swq;
            }

            /** record snowpack depth **/
            if (output_var[OUT_SNOW_DEPTH]) {
                out_data[OUT_SNOW_DEPTH][veg][0] = snow[veg].snow_depth;
            }

            size_t Nsnow = snow[veg].Nsnow;
            /** record snowpack temperature water and ice content **/
            if (output_var[OUT_SNOW_PACK_TEMP]) {
                for (i = 0; i < Nsnow; i++) {
                    out_data[OUT_SNOW_PACK_TEMP][veg][i] = snow[veg].pack_T[i];
                }
            }

            if (output_var[OUT_SNOW_PACK_ICE]) {
                for (i = 0; i < Nsnow; i++) {
                    out_data[OUT_SNOW_PACK_ICE][veg][i] = snow[veg].pack_ice[i];
                }
            }

            if (output_var[OUT_SNOW_PACK_LIQ]) {
                for (i = 0; i < Nsnow; i++) {
                    out_data[OUT_SNOW_PACK_LIQ][veg][i] = snow[veg].pack_liq[i];
                }
            }

            if (output_var[OUT_PACK_OUTFLOW]) {
                for (i = 0; i < Nsnow; i++) {
                    out_data[OUT_PACK_OUTFLOW][veg][i] = snow[veg].pack_outflow[i];
                }
            }

            if (output_var[OUT_SNOW_ICEFRAC]) {
                for (i = 0; i < Nsnow; i++) {
                    out_data[OUT_SNOW_ICEFRAC][veg][i] = snow[veg].theta_ice[i];
                }
            }

            if (output_var[OUT_SNOW_LIQFRAC]) {
                for (i = 0; i < Nsnow; i++) {
                    out_data[OUT_SNOW_LIQFRAC][veg][i] = snow[veg].theta_liq[i];
                }
            }

            if (output_var[OUT_SNOW_POROSITY]) {
                for (i = 0; i < Nsnow; i++) {
                    out_data[OUT_SNOW_POROSITY][veg][i] = snow[veg].porosity[i];
                }
            }

            /** record snow density **/
            if (output_var[OUT_SNOW_DENSITY]) {
                for (i = 0; i < Nsnow; i++) {
                    out_data[OUT_SNOW_DENSITY][veg][i] = snow[veg].density[i];
                }
            }

            /** record snowpack melt **/
            if (output_var[OUT_SNOW_MELT]) {
                for (i = 0; i < Nsnow; i++) {
                    out_data[OUT_SNOW_MELT][veg][i] = snow[veg].pack_melt[i];
                }
            }

            /** record snow surface freezing rate [mm/s] */
            if (output_var[OUT_SNOW_FRZE]) {
                for (i = 0; i < Nsnow; i++) {
                    out_data[OUT_SNOW_FRZE][veg][i] = snow[veg].pack_frze[i];
                }
            }

            /** record snow surface evaporation **/
            if (output_var[OUT_SNOW_EVAP]) {
                out_data[OUT_SNOW_EVAP][veg][0] = snow[veg].snow_evap * dt_sec * MM_PER_M;
            }

            /** record snowpack combination **/
            if (output_var[OUT_SNOW_COMB]) {
                out_data[OUT_SNOW_COMB][veg][0] = snow[veg].pack_comb;
            }

            /** record snow surface frost [mm] */
            if (output_var[OUT_SNOW_FROST]) {
                out_data[OUT_SNOW_FROST][veg][0] = snow[veg].snow_frost * dt_sec * MM_PER_M;
            }

            /** record snow surface dew [mm] */
            if (output_var[OUT_SNOW_DEW]) {
                out_data[OUT_SNOW_DEW][veg][0] = snow[veg].snow_dew * dt_sec * MM_PER_M;
            }

            /** record snow surface sublimation [mm] */
            if (output_var[OUT_SNOW_SUBLIM]) {
                out_data[OUT_SNOW_SUBLIM][veg][0] = snow[veg].snow_sublim * dt_sec * MM_PER_M;
            }

            /** record glacier snow excess flow **/
            if (output_var[OUT_GLAC_EXCESS]) {
                out_data[OUT_GLAC_EXCESS][veg][0] = snow[veg].glac_excess;
            }

            /** record snow cover fraction **/
            if (output_var[OUT_SNOW_COVER]) {
                out_data[OUT_SNOW_COVER][veg][0] = snow[veg].coverage;
            }

            /** record snow depth increasing rate **/
            if (output_var[OUT_DELDEPTH]) {
                out_data[OUT_DELDEPTH][veg][0] = snow[veg].delta_depth;
            }

            /** record NEW snow density **/
            if (output_var[OUT_NEW_DENSITY]) {
                out_data[OUT_NEW_DENSITY][veg][0] = snow[veg].new_snow_density;
            }

            /** record SnowAge **/
            if (output_var[OUT_SNOW_AGE]) {
                out_data[OUT_SNOW_AGE][veg][0] = snow[veg].snowage;
            }

            /** outflow of liquid water from the snowpack bottom (mm) */
            if (output_var[OUT_SNOW_OUTFLOW]) {
                out_data[OUT_SNOW_OUTFLOW][veg][0] = snow[veg].snow_outflow * dt_sec;
            }

            // Glacier Water Balance Terms
            if (HasGlac) {
                if (output_var[OUT_GLAC_INFLOW]) {
                    out_data[OUT_GLAC_INFLOW][veg][0] = cell[veg].soil_inflow;
                }
                if (output_var[OUT_GLAC_OUTFLOW]) {
                    out_data[OUT_GLAC_OUTFLOW][veg][0] = cell[veg].runoff + cell[veg].baseflow;
                }
            }

            /**********************************
             Record Energy Balance Variables
            **********************************/           
            if (output_var[OUT_TRND_FBFLAG]) {
                out_data[OUT_TRND_FBFLAG][veg][0] = energy[veg].FrozenGrnd;
            }
            if (HasVeg) {
                if (output_var[OUT_TCAN_FBFLAG]) {
                    out_data[OUT_TCAN_FBFLAG][veg][0] = energy[veg].FrozenOver;
                }
            }

            /** record landcover temperature **/
            if (HasVeg) {
                // landcover is vegetation
                if (output_var[OUT_VEGT]) {
                    out_data[OUT_VEGT][veg][0] = energy[veg].Tfoliage;
                }
                if (output_var[OUT_VEGTAIR]) {
                    out_data[OUT_VEGTAIR][veg][0] = energy[veg].Tcanopy;
                }
            }

            // landcover is bare soil
            if (output_var[OUT_BARESOILT]) {
                out_data[OUT_BARESOILT][veg][0] = energy[veg].Tgrnd;
            }

            /** record mean surface temperature [C] **/
            if (output_var[OUT_SURF_TEMP]) {
                out_data[OUT_SURF_TEMP][veg][0] = energy[veg].Tsurf;
            }

            /** record NODES temperatures **/
            if (output_var[OUT_SOIL_TEMP]) {
                for (i = 0; i < cell[veg].Nsoil; i++) {
                    out_data[OUT_SOIL_TEMP][veg][i] = cell[veg].soil_T[i];
                }
            }

            /** record advective flux from prec **/
            if (output_var[OUT_ADVECTION]) {
                out_data[OUT_ADVECTION][veg][0] = energy[veg].advection;
            }

            if (output_var[OUT_ADVECTGRND]) {
                out_data[OUT_ADVECTGRND][veg][0] = energy[veg].AdvectGrnd;
            }

            if (HasVeg) {
                if (output_var[OUT_ADVECTSUB]) {
                    out_data[OUT_ADVECTSUB][veg][0] = energy[veg].AdvectSub;
                }
                if (output_var[OUT_ADVECTOVER]) {
                    out_data[OUT_ADVECTOVER][veg][0] = energy[veg].AdvectOver;
                }
            }

            /** record net shortwave radiation **/
            if (output_var[OUT_SWNET]) {
                out_data[OUT_SWNET][veg][0] = energy[veg].shortwave;
            }

            /** record longwave radiation flux **/
            if (output_var[OUT_LWNET]) {
                out_data[OUT_LWNET][veg][0] = energy[veg].longwave;
            }

            /** record latent heat flux **/
            if (output_var[OUT_LATENT]) {
                out_data[OUT_LATENT][veg][0] = energy[veg].latent;
            }

            /** record sensible heat flux **/
            if (output_var[OUT_SENSIBLE]) {
                out_data[OUT_SENSIBLE][veg][0] = energy[veg].sensible;
            }

            /** record ground heat flux **/
            if (output_var[OUT_GRND_FLUX]) {
                out_data[OUT_GRND_FLUX][veg][0] = energy[veg].grnd_flux;
            }

            // 冰川变量
            if (HasGlac) {
                if (output_var[OUT_H2OSFC_T]) {
                    out_data[OUT_H2OSFC_T][veg][0] = energy[veg].Tgrnd;
                }
            }
        }
    } // End loop over veg

    // vic_run run time
    if (output_var[OUT_TIME_VICRUN_WALL]) {
        out_data[OUT_TIME_VICRUN_WALL][0][0] = timer->delta_wall;
    }
    if (output_var[OUT_TIME_VICRUN_CPU]) {
        out_data[OUT_TIME_VICRUN_CPU][0][0] = timer->delta_cpu;
    }
}