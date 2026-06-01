/******************************************************************************
 * @section DESCRIPTION
 *
 *  This subroutine calculates the maximum daylight duration based on the 
 *  obliquity of the ecliptic as defined by the IAU 1976 standard.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    This subroutine computes maximum daylight duration.
 *****************************************************************************/
double
calc_max_daylength(double lat)
{
    double coslat;
    double sinlat;
    double decl = 0.0;
    double cosdecl;
    double sindecl;
    double cosegeom;
    double sinegeom;
    double coshss;

    /* calculate cos and sin of latitude */
    coslat = cos(lat * CONST_PI / 180);
    sinlat = sin(lat * CONST_PI / 180);

    /* calculate cos and sin of declination */
    if (lat < 0.0) {
        decl = -CONST_OBLIQUITY;
    }
    cosdecl = cos(decl);
    sindecl = sin(decl);

    /* calculate daylength as a function of lat and decl */
    cosegeom = coslat * cosdecl;
    sinegeom = sinlat * sindecl;
    coshss = -(sinegeom) / cosegeom;
    if (coshss < -1.0) {
        coshss = -1.0; /* 24-hr daylight */
    }
    if (coshss > 1.0) {
        coshss = 1.0; /* 0-hr daylight */
    }
    // max daylength
    return 2.0 * CONST_SECPERRAD * acos(coshss);

}