/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine computes the fraction of roots in each soil layer based on the
 * root zone distribution defined in the vegetation parameter file.  Roots are
 * assumed to be linearly distributed within each root zone.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    This routine computes the fraction of roots in each soil layer.
 *****************************************************************************/
void
calc_root_fractions(veg_con_struct  *veg_con,
                    soil_con_struct *soil_con)
{
    extern option_struct options;

    int                  layer;
    size_t               veg, Nveg;
    unsigned short       veg_class;
    double               a;          // Empirical parameter a in eqa(2)
    double               b;          // Empirical parameter b in eqa(2)
    double               d;          // Maximum root depth (m)
    double               Y;
    double               sum_fract;
    double               Zsum;
    double               dum;

    /* Set number of vegetation tiles */
    Nveg = veg_con[0].vegetat_type_num;

    for (veg = 0; veg <= Nveg; veg++) {
        if (veg_con[veg].Cv > 0.0) {
            
            veg_class = veg_con[veg].veg_class;

            if (veg_class != options.GLACIER_ID) {

                // Check that root fractions sum to 1.0
                sum_fract = 0.0;
                Zsum = 0.0;
                a = veg_con[veg].a;
                b = veg_con[veg].b;
                d = veg_con[veg].d;

                int last_node = -1;
                for (layer = 0; layer < options.Nnode; layer++) {
                    Zsum += soil_con->dz_soil[layer];
                    if (Zsum > d) {
                        last_node = layer;
                        break;
                    }
                }
                if (last_node == -1) {
                    last_node = options.Nnode - 1;
                }
                veg_con[veg].Nroot = last_node + 1;
                // Reset for actual root fraction calculation
                Zsum = 0.0;
                // Calculate root fraction for each soil layer
                for (layer = 0; layer < last_node; layer++) {

                    // Current layer bottom depth
                    Zsum += soil_con->dz_soil[layer];

                    double Zused = min(Zsum, d);
                    // Cumulative root fraction to top and bottom of this layer
                    Y = 1.0 - 0.5 * (exp(-a * Zused) + exp(-b * Zused));

                    // Single layer root fraction = difference
                    veg_con[veg].root[layer] = Y - sum_fract;

                    // Accumulate for normalization check
                    sum_fract += veg_con[veg].root[layer];

                }
                Zsum += soil_con->dz_soil[last_node];
                veg_con[veg].root[last_node] = 1.0 - sum_fract;
                // Final check on root fractions. If they don't sum to 1, throw error
                // Otherwise, rescale by sum to eliminate small rounding errors
                dum = 0.;
                for (layer = 0; layer < last_node + 1; layer++) {
                    if (veg_con[veg].root[layer] < 1.e-4) {
                        veg_con[veg].root[layer] = 0.;
                    }
                    dum += veg_con[veg].root[layer];
                }
                if (!assert_close_double(dum, 1, 0, 1e-4)) {
                    log_err("Soil layer root fractions do not sum to 1.0: %f, "
                            "veg class: %d", dum, veg_con[veg].veg_class);
                }
                else {
                    for (layer = 0; layer < last_node + 1; layer++) {
                        veg_con[veg].root[layer] /= dum;
                    }
                }
            }
            else {
                veg_con[veg].Nroot = 0;
            }
        }
    }
}
    

