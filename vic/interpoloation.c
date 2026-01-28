/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine performs a interpolation
 *****************************************************************************/

#include <vic_run.h>

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
 * @brief    定义sign函数
 *****************************************************************************/
double sign(double a, double b) {
    return (b >= 0) ? fabs(a) : -fabs(a);
}
