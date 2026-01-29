/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine computes the evaporation, traspiration and throughfall of the
 * vegetation types for multi-layered model.
 *
 * The value of x, the fraction of precipitation that exceeds the canopy
 * storage capacity, is returned by the subroutine.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    Calculation of evaporation from the canopy, including the
 *           possibility of potential evaporation exhausting ppt+canopy storage
 *****************************************************************************/
double
canopy_evap(double             step_dt,
            energy_bal_struct *energy,
            cell_data_struct  *cell,
            snow_data_struct  *snow,
            veg_var_struct    *veg_var)
{
    /** declare global variables **/
    extern parameters_struct param;

    /** declare local variables **/
    double          snow_canopy;
    double          transp;
    double          canopy_sublim;
    double          canopydew;
    double          canopy_vapor;
    double          canopyfrost;
    double          canopy_melt;
    double          canopy_freez;
    double          Evap;
    double          wetFrac;

    /* Initialization */
    double Tcanopy = energy->Tcanopy;
    double int_snow = veg_var->int_snow;
    double int_rain = veg_var->int_rain;
    double LatentTransp = energy->LatentTransp;
    double LatentCanopy = energy->LatentCanopy;

    /* Canopy Hydrology processes */
    if (energy->FrozenOver == true) {
        transp = max(LatentTransp / CONST_LATSUB, 0.0);
        canopy_vapor = 0.0;
        canopydew = 0.0;
        canopy_sublim = max(LatentCanopy / CONST_LATSUB, 0.0);
        canopyfrost = fabs(min(LatentCanopy / CONST_LATSUB, 0.0));
    }
    else {
        transp = max(LatentTransp / CONST_LATVAP, 0.0);
        canopy_vapor = max(LatentCanopy / CONST_LATVAP, 0.0);
        canopydew = fabs(min(LatentCanopy / CONST_LATVAP, 0.0));
        canopy_sublim = 0.0;
        canopyfrost = 0.0;        
    }
    if (canopy_vapor > int_rain / step_dt) {
        canopy_vapor = int_rain / step_dt;
    }
    int_rain = max(0.0, int_rain + (canopydew - canopy_vapor) * step_dt);
    if (int_rain < param.TOL_A) {
        int_rain = 0.0;
    }

    /* canopy ice */
    if (canopy_sublim > int_snow / step_dt) {
        canopy_sublim = int_snow / step_dt;
    }
    int_snow = max(0.0, int_snow + (canopyfrost - canopy_sublim) * step_dt);
    if (int_snow < param.TOL_A) {
        int_snow = 0.0;
    }

    /* wetted fraction of canopy */
    if (int_snow > 0.0 && int_snow > int_rain) {
        wetFrac = max(0.0, int_snow) / max(veg_var->MaxSnowInt, param.TOL_A);
    }
    else {
        wetFrac = max(0.0, int_rain) / max(veg_var->MaxRainInt, param.TOL_A);
    }
    wetFrac = pow(min(wetFrac, 1.0), 0.667);
    snow_canopy = int_snow + int_rain;

    /* phase change */
    if (int_snow > param.TOL_A && Tcanopy > CONST_TKFRZ) {
        canopy_melt = min(int_snow / step_dt, (Tcanopy - CONST_TKFRZ) * CONST_CPICE *
                                int_snow / CONST_RHOICE / (step_dt * CONST_LATICE));
        int_snow = max(0.0, int_snow - canopy_melt * step_dt);
        int_rain = max(0.0, snow_canopy - int_snow);
        Tcanopy = wetFrac * CONST_TKFRZ + (1.0 - wetFrac) * Tcanopy;
    }

    /* canopy water refreeezing */
    if (int_rain > param.TOL_A && Tcanopy < CONST_TKFRZ) {
        canopy_freez = min(int_rain / step_dt, (CONST_TKFRZ - Tcanopy) * CONST_CPFW * 
                                int_rain / CONST_RHOFW / (step_dt * CONST_LATICE));
        int_rain = max(0.0, int_rain - canopy_freez * step_dt);
        int_snow = max(0.0, snow_canopy - int_rain);
        Tcanopy = wetFrac * CONST_TKFRZ + (1.0 - wetFrac) * Tcanopy;
    }

    /* update total canopy water */
    snow->snow_canopy = int_snow + int_rain;
    veg_var->int_snow = int_snow;
    veg_var->int_rain = int_rain;
    cell->transp = transp;
    cell->canopydew = canopydew;
    cell->canopyfrost = canopyfrost;
    cell->canopy_sublim = canopy_sublim;
    cell->canopy_vapor = canopy_vapor;
    energy->Tcanopy = Tcanopy;
    Evap = canopy_vapor + canopy_sublim - canopydew - canopyfrost;
    
    return Evap;
}