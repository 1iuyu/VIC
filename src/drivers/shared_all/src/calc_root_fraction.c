/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine computes the fraction of roots in each soil layer based on the
 * root zone distribution defined in the vegetation parameter file.  Roots are
 * assumed to be linearly distributed within each root zone.
 *****************************************************************************/

#include "vic_driver_shared_all.h"

/******************************************************************************
 * @brief    This routine computes the fraction of roots in each soil layer.
 *****************************************************************************/
void
calc_root_fractions(size_t            veg_class,
                    cell_data_struct *cell,
                    veg_var_struct   *veg_var,
                    soil_con_struct  *soil_con,
                    veg_lib_struct   *veg_lib)
{
    /* initialization */
    extern parameters_struct param;
    size_t i;
    int    layer;
    double a;          // Empirical parameter a in eqa(2)
    double b;          // Empirical parameter b in eqa(2)
    double d;          // Maximum root depth (m)
    double Y;
    double sum_fract = 0.0;
    size_t Nsoil = soil_con->Nbedrock - 1;

     /* Set number of vegetation tiles */
    if (cell->IS_VEG) {
        
        a = veg_lib[veg_class].root_a;
        b = veg_lib[veg_class].root_b;
        d = veg_lib[veg_class].root_d;
        int last_node = -1;
        double Zsum = 0.0;

        for (i = 0; i < Nsoil; i++) {
            Zsum += soil_con->dz_soil[i];
            if (Zsum > d) {
                last_node = i;
                break;
            }
        }
        
        if (last_node == -1) {
            last_node = Nsoil - 1;
        }
        veg_var->Nroot = last_node + 1;
        // 计算各层根系分数
        Zsum = 0.0;
        for (layer = 0; layer <= last_node; layer++) {
            Zsum += soil_con->dz_soil[layer];
            double Zused = min(Zsum, d);
            Y = 1.0 - 0.5 * (exp(-a * Zused) + exp(-b * Zused));
            cell->root[layer] = Y - sum_fract;
            sum_fract += cell->root[layer];
        }
        // 检查根系参数是否有效
        if (sum_fract <= 0.0) {
            log_err("Invalid root distribution parameters: "
                    "total root fraction = %.4f (must be > 0), "
                    "veg_class = %zu, a = %.4f, b = %.4f, d = %.4f. ",
                     sum_fract, veg_class, a, b, d);
        }
        if (fabs(sum_fract - 1.0) > param.TOL_A) {
            for (layer = 0; layer <= last_node; layer++) {
                cell->root[layer] /= sum_fract;
            }
        }
        // Final check on root fractions. If they don't sum to 1, throw error
        // Otherwise, rescale by sum to eliminate small rounding errors
        double dum = 0.0;
        for (layer = 0; layer <= last_node; layer++) {
            if (cell->root[layer] < 1.e-4) {
                cell->root[layer] = 0.0;
            }
            dum += cell->root[layer];
        }
        if (!assert_close_double(dum, 1, 0, 1e-4)) {
            log_err("Soil layer root fractions do not sum to 1.0: %f, "
                    "veg class: %zu", dum, veg_class);
        }
        else {
            if (dum != 1.0) {
                for (layer = 0; layer <= last_node; layer++) {
                    cell->root[layer] /= dum;
                }
            }
        }
    }
    else {
        veg_var->Nroot = 0;
    }
}