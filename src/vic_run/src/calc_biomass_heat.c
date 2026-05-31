/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine calculate the biomass heat capacities.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    Calculate the biomass heat capacities.
 *****************************************************************************/
int
calc_biomass_heat(double         *f_abs_stem,
                  double         *sa_leaf,
                  double         *sa_stem,
                  double         *cp_leaf,
                  double         *cp_stem,
                  double         *stem_biomass,
                  double         *Ra_stem,
                  veg_var_struct *veg_var,
                  veg_lib_struct *veg_lib)
{
    extern option_struct options;
    double Canopy_Upper = veg_lib->Canopy_Upper;
    double NetLAI = veg_var->NetLAI;
    double NetSAI = veg_var->NetSAI;
    double leaf_biomass = 0.0;
    double wood_density = 500.0; // kg/m3
    double rstem = 1000.0; // stem resistance per stem diameter (s/m2)
    double liq_bioms = veg_lib->liq_bioms; // fraction of biomass that is water, for calculating heat capacity
    double carea_stem = 0.0;
    if (options.BIOMASST) {
        (*f_abs_stem) = NetSAI / (NetLAI + NetSAI);
        if (NetLAI > 0.0) {
            (*f_abs_stem) = 0.1 * (*f_abs_stem);
        }
        *sa_leaf = 2.0 * NetLAI;
        *sa_stem = veg_lib->stem_num * (Canopy_Upper * CONST_PI * veg_lib->trunk_dia);
        if (Canopy_Upper <= 1.0 || veg_lib->trunk_dia < 0.05) {
            (*f_abs_stem) = 0.0;
            *sa_stem = 0.0;
            *sa_leaf += NetSAI; 
        }
        else {
            if (NetLAI < MIN_VEG_LAI) {
                *sa_leaf += NetSAI;
            }
        }
        leaf_biomass = (1.0e-3 * 2.0 / veg_lib->slatop) * max(0.01, 0.5 * (*sa_leaf)) / (1.0 - liq_bioms);
        carea_stem = CONST_PI * pow(veg_lib->trunk_dia * 0.5, 2.0);
        *stem_biomass = carea_stem * Canopy_Upper * veg_lib->stem_num *
                       wood_density / (1.0 - liq_bioms);
        // calculate specify heat capacity of vegetation
        // as weighted averaged of dry biomass and water
        // lma_dry has units of kg dry mass/m2 here
        // (Appendix B of Bonan et al., GMD, 2018) 

        *cp_leaf  = leaf_biomass * (1400.0 * (1.0 - liq_bioms) + liq_bioms * CONST_CPFWICE);

        // cp-stem will have units J/k/ground_area
        *cp_stem = *stem_biomass * (1400.0 * (1.0 - liq_bioms) + liq_bioms * CONST_CPFWICE);

        // resistance between internal stem temperature and canopy air 
        *Ra_stem = rstem * veg_lib->trunk_dia;
    }
    else {
        (*f_abs_stem) = 0.0;
        *sa_leaf = NetLAI + NetSAI;
        *sa_stem = 0.0;
        *cp_leaf = 0.0;
        *cp_stem = 0.0;
        *stem_biomass = 0.0;
        *Ra_stem = 0.0;
    }
    return (0);
}