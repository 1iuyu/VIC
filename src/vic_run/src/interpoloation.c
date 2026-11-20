/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine performs a interpolation
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    This routine performs a linear interpolation
 *****************************************************************************/
double
linear_interp(double x,
              double lx,
              double ux,
              double ly,
              double uy)
{
    return((x - lx) / (ux - lx) * (uy - ly) + ly);
}
/******************************************************************************
 * @brief    This routine performs a piecewise linear interpolation method
 *           for 1-dimensional data.
 *****************************************************************************/
double 
piecewise_linear_interp(size_t  nd,
                        double *xd,
                        double *yd,
                        double  xi)
{
    double yi = 0.0;
    double t;
    size_t k;
    
    // 只有一个数据点的情况
    if (nd == 1) {
        return yd[0];
    }
    
    // 多个数据点的情况
    if (xi < xd[0]) {
        // 左边界外推
        t = (xi - xd[0]) / (xd[1] - xd[0]);
        yi = (1.0 - t) * yd[0] + t * yd[1];
    } else if (xi > xd[nd - 1]) {
        // 右边界外推
        t = (xi - xd[nd - 2]) / (xd[nd - 1] - xd[nd - 2]);
        yi = (1.0 - t) * yd[nd - 2] + t * yd[nd - 1];
    } else {
        // 区间内分段线性插值
        for (k = 1; k < nd; k++) {
            if (xd[k - 1] <= xi && xi <= xd[k]) {
                t = (xi - xd[k - 1]) / (xd[k] - xd[k - 1]);
                yi = (1.0 - t) * yd[k - 1] + t * yd[k];
                break;
            }
        }
    }
    
    return yi;
}

/******************************************************************************
 * @brief    定义sign函数
 *****************************************************************************/
double sign(double a, double b) {
    return (b >= 0) ? fabs(a) : -fabs(a);
}

/******************************************************************************
 * @brief    Get the index in the band-diagonal matrix or full matrix
 *****************************************************************************/
size_t
get_band_idx(size_t row, size_t col)
{
    // 检查是否超出矩阵带宽
    if (row >= col) {
        if (row - col > LAPACK_KL) {
            log_err("ERROR: Matrix index out of lower band! "
                    "row=%zu, col=%zu, KL=%d",
                    row, col, LAPACK_KL);
            return (ERROR);
        }
    }
    else {
        if (col - row > LAPACK_KU) {
            log_err("ERROR: Matrix index out of upper band! "
                    "row=%zu, col=%zu, KU=%d",
                    row, col, LAPACK_KU);
            return (ERROR);
        }
    }

    return (LAPACK_KL + LAPACK_KU + row - col) + col * LAPACK_LDAB;
}

/******************************************************************************
 * @brief    Calculate the gradient: grad = J^T * rVec
 *****************************************************************************/
void
calc_gradient(double *mat_LHS, const double *mat_RHS, double *grad, size_t Nstate)
{
    size_t i, j, jmin, jmax, k;

    for (i = 0; i < Nstate; i++) {
        grad[i] = 0.0;

        jmin = (i > LAPACK_KU) ? i - LAPACK_KU : 0;
        jmax = min(Nstate - 1, i + LAPACK_KL);

        for (j = jmin; j <= jmax; j++) {
            k = get_band_idx(j, i);
            grad[i] += mat_LHS[k] * mat_RHS[j];
        }
    }
}
