/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine calculates evapotranspiration using the Penman-Monteith
 * approach.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    Calculate canopy stomatal resistance
 *****************************************************************************/
int
calc_rc(size_t             index,
        double             VPcanopy,
        double             pressure,
        double             total_transp,
        energy_bal_struct *energy,
        veg_var_struct    *veg_var,
        veg_lib_struct    *veg_lib)
{   
    extern option_struct     options;
    extern parameters_struct param;

    if (options.RC_MODE == RC_JARVIS) {
        double      Sfactor;    /* factor for canopy resistance based on photosynthesis */
        double      Tfactor;    /* factor for canopy resistance based on temperature */
        double      vpdfactor;  /* factor for canopy resistance based on vpd */
        double      RS_tmp;

        double ratio = 0.0;
        double ratio_slope = 0.0;
        double Tcanopy = energy->Tcanopy;
        double APAR_sunlit = energy->APAR_sunlit;
        double APAR_shade = energy->APAR_shade;
        double fcanopy = veg_var->fcanopy;
        double RGL = veg_lib->RGL;
        double rmax = veg_lib->rmax;
        double rmin = veg_lib->rmin;
        double m_vpd = veg_lib->m_vpd;
        double T_opt_trans = veg_lib->T_opt_trans;
        double tmp_APAR = 0.;
        if (index == 0) {
            tmp_APAR = APAR_sunlit / max(fcanopy, param.TOL_A);
        }
        if (index == 1) {
            tmp_APAR = APAR_shade / max(fcanopy, param.TOL_A);
        }

        double tmp_Qair = 0.622 * VPcanopy / (pressure - 0.378 * VPcanopy); // specific humidity
        double tmp_ratio = tmp_Qair / (1.0 - tmp_Qair);

        s_humid(Tcanopy, pressure, &ratio, &ratio_slope);

        /* solar radiation factor */
        double f_rad = 2.0 * tmp_APAR / RGL;
        Sfactor = (f_rad + rmin / rmax) / (1.0 + f_rad);
        Sfactor = max(Sfactor, 0.0001);

        /* air temperature factor */
        Tfactor = 1.0 - 0.0016 * pow(T_opt_trans - Tcanopy, 2.0);
        Tfactor = max(Tfactor, 0.0001);

        /* vpd factor */
        vpdfactor = 1.0 / (1.0 + m_vpd * max(0.0, ratio - tmp_ratio));
        vpdfactor = max(vpdfactor, 0.01);
        RS_tmp = rmin / (Sfactor * Tfactor * vpdfactor * total_transp);

        // Maximal stomatal resistance [s/m]  
//        if (RS_tmp > 50000.0) {
//            RS_tmp = 50000.0;
//        }
        /* assign updated values */
        if (index == 0) {
            veg_var->RS_sunlit = RS_tmp;
        }
        else {
            veg_var->RS_shade = RS_tmp;
        }
    }
    else if (options.RC_MODE == RC_PHOTO) {

/*      double CF = pressure / (8.314 * air_temp) * param.MAX_LIMIT;
        double RS_tmp = 1.0 / g_min * CF;
        double tmp_photoleaf = 0.;
        if (index == 0) {
            tmp_APAR = APAR_sunlit / max(fcanopy, param.TOL_A);
        }
        if (index == 1) {
            tmp_APAR = APAR_shade / max(fcanopy, param.TOL_A);
        }

        if (tmp_APAR > 0.) {
            f_N = min()
        }
        else {
            tmp_photoleaf = 0.;
        }
*/
    }


    return (0);
}
