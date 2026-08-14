/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine computes the cosine of the solar zenith angle, given the
 * current location and date.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    This subroutine computes the cosine of the solar zenith angle.
 *****************************************************************************/
int
compute_coszen(double         lat,
               double         lng,
               double         time_zone_lng,
               unsigned short day_in_year,
               unsigned       second,
               double        *coszen,
               double        *daylen)
{
    extern option_struct options;
    double coslat;
    double sinlat;
    double decl;
    double cosdecl;
    double sindecl;
    double cosegeom;
    double sinegeom;
    double coshss;
    double cosh;
    double hour;
    double local_solar_time = 0.0;

    hour = second / SEC_PER_HOUR;

    /* calculate cos of hour angle */
    if (options.FORCE_TIME == FORCE_LOCAL) {
        local_solar_time = hour;
    }
    else if (options.FORCE_TIME == FORCE_UTC) {
        local_solar_time = hour + (time_zone_lng / 15.0) + 
                           (lng - time_zone_lng) / 15.0;      
    }
    while (local_solar_time < 0.0)  {
        local_solar_time += 24.0;
    }
    while (local_solar_time >= 24.0) {
        local_solar_time -= 24.0;
    }

    /* calculate cos and sin of latitude */
    coslat = cos(lat * CONST_PI / 180);
    sinlat = sin(lat * CONST_PI / 180);

    /* calculate cos and sin of declination */
    decl = CONST_MINDECL * cos(((double) day_in_year + CONST_DAYSOFF) *
                               CONST_RADPERDAY);
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
    (*daylen) = 2.0 * CONST_SECPERRAD * acos(coshss);

    cosh = cos((local_solar_time - 12) * CONST_PI / 12);

    /* calculate cosine of solar zenith angle */
    (*coszen) = cosegeom * cosh + sinegeom;

    return (0);
}
