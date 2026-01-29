/******************************************************************************
 * @section DESCRIPTION
 *
 * Iteratively solve the atmospheric energy balance equation to estimate the
 * canopy air temperature.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    Iteratively solve the atmospheric energy balance equation to
 *           estimate the canopy air temperature.
 *****************************************************************************/
double
calc_atmos_energy_bal(double             step_dt,
                      double             air_temp,
                      double             air_density,
                      double             esat_Tsub,
                      double             vp,
                      double             coef_lw_atmos,
                      double             coef_lw_canopy,
                      double             c_biomass,
                      double             esat_Tcanopy,
                      double             esat_slope,
                      double            *delta_over,
                      energy_bal_struct *energy,
                      cell_data_struct  *cell,
                      veg_var_struct    *veg_var)
{
    /* Initialization */
    double Tupper = energy->Tupper;
    double Tcanopy = energy->Tcanopy;
    double Ra_leaf = cell->Ra_leaf;
    double Ra_evap = cell->Ra_evap;
    double *Ra_over = cell->Ra_over;
    double *Ra_sub = cell->Ra_sub;

    double fcanopy = veg_var->fcanopy;
    double wetFrac = veg_var->wetFrac;
    double int_snow = veg_var->int_snow;
    double int_rain = veg_var->int_rain;
    double RS_shade = veg_var->RS_shade;
    double RS_sunlit = veg_var->RS_sunlit;
    double PsyCh_canopy = energy->PsyCh_canopy;
    double LatentVapOver = energy->LatentVapOver;
    double leaf_sun = min(6.0, veg_var->leaf_sun);
    double leaf_shade = min(6.0, veg_var->leaf_shade);
    double NetVAI = min(6.0, veg_var->NetLAI + veg_var->NetSAI);

    // compute incoming sensible heat
    double coef_SHover = 1.0 / Ra_over[1];
    double coef_SHleaf = 2.0 * NetVAI / Ra_leaf;
    double coef_ground = 1.0 / Ra_sub[1];
    double coef_flux = coef_SHover + coef_SHleaf + coef_ground;
    double SH_temp = (air_temp * coef_SHover +
                             Tupper * coef_ground) / coef_flux;
    double SH_frac = coef_SHleaf / coef_flux;
    double coef_sensible = (1.0 - SH_frac) * air_density *
                                    CONST_CPDAIR * coef_SHleaf;

    // compute total latent heat flux
    double coef_LHover = 1.0 / Ra_over[2];
    double coef_LHevap = wetFrac * NetVAI / Ra_leaf;
    double coef_LHtransp = (1.0 - wetFrac) * 
                                (leaf_sun / (Ra_leaf + RS_sunlit) +
                                    leaf_shade / (Ra_leaf + RS_shade));
    double coef_LHgrnd = 1.0 / (Ra_sub[2] + Ra_evap);
    coef_flux = coef_LHover + coef_LHevap + coef_LHtransp + coef_LHgrnd;
    double LH_vp = (vp * coef_LHover + esat_Tsub * coef_LHgrnd) / coef_flux;
    double coef_evapfrac = (coef_LHevap + coef_LHtransp) / coef_flux;
    double coef_latent = (1.0 - coef_evapfrac) * coef_LHevap * air_density * 
                                                CONST_CPDAIR / PsyCh_canopy;
    double coef_transp = (1.0 - coef_evapfrac) * coef_LHtransp * air_density * 
                                                CONST_CPDAIR / PsyCh_canopy;
    double Tair = SH_temp + SH_frac * Tcanopy;
    double VPcanopy = LH_vp + coef_evapfrac * esat_Tcanopy;

    double NetLongOver = fcanopy * (coef_lw_atmos + coef_lw_canopy * pow(Tcanopy, 4.0));
    double SensibleOver = fcanopy * air_density * CONST_CPDAIR * 
                                              coef_SHleaf * (Tcanopy - Tair);

    double LatentCanopy = fcanopy * air_density * CONST_CPDAIR * coef_LHevap * 
                                (esat_Tcanopy - VPcanopy) / PsyCh_canopy;
    double LatentTransp = fcanopy * air_density * CONST_CPDAIR * coef_LHtransp * 
                                (esat_Tcanopy - VPcanopy) / PsyCh_canopy;
    if (Tcanopy > CONST_TKFRZ) {
        LatentCanopy = min(int_rain * LatentVapOver / step_dt, LatentCanopy);
    }
    else {
        LatentCanopy = min(int_snow * LatentVapOver / step_dt, LatentCanopy);
    }

    /* canopy heat capacity */
    double C_canopy = c_biomass * NetVAI *
                    CONST_CPFWICE + int_rain * CONST_CPFWICE /
                        CONST_RHOFW + int_snow * CONST_CPICE / CONST_RHOICE;  // [J/m2/K]

    double RestTerm = energy->NetShortSub - NetLongOver -
                    SensibleOver - LatentCanopy -
                            LatentTransp + energy->AdvectOver;
    double tmp_coef = fcanopy * (4.0 * coef_lw_canopy * pow(Tcanopy, 3) + coef_sensible +
                        (coef_latent + coef_transp) * esat_slope + C_canopy / step_dt);
    (*delta_over) = RestTerm / tmp_coef;

    /* update fluxes with temperature change */
    NetLongOver += fcanopy * 4.0 * coef_lw_canopy * pow(Tcanopy, 3.0) * (*delta_over);
    SensibleOver += fcanopy * coef_sensible * (*delta_over);
    LatentCanopy += fcanopy * coef_latent * esat_slope * (*delta_over);
    LatentTransp += fcanopy * coef_transp * esat_slope * (*delta_over);
    double delta_Qc = fcanopy * C_canopy / step_dt * (*delta_over);
    /* update Tcanopy */
    Tcanopy += (*delta_over);

    energy->LatentCanopy = LatentCanopy;
    energy->LatentTransp = LatentTransp;
    energy->NetLongOver = NetLongOver;
    energy->SensibleOver = SensibleOver;
    energy->Tair = Tair;
    cell->VPcanopy = VPcanopy;

    return(Tcanopy);
}
