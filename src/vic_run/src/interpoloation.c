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
