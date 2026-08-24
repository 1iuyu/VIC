/******************************************************************************
* @section DESCRIPTION
*
* This routine was written to handle the various calls and data
* handling needed to solve the various components of the new VIC
* snow code for both the full_energy and water_balance models.
******************************************************************************/

#include "rout.h"

/******************************************************************************
* @brief        This routine was written to handle the various calls and data
*               handling needed to solve the various components of the new VIC
*               snow code for both the full_energy and water_balance models.
******************************************************************************/
int
Euler_Routing(size_t i,
              double step_dt)
{
    extern parameters_struct param;
    extern rout_struct      *rout;
    extern domain_struct local_domain;

    size_t j, k, m;
    size_t river_steps = rout[i].river_steps;
    size_t sub_steps = rout[i].sub_steps;
    double localDeltaT;

    /********************************
      Hillslope routing
    ********************************/
    rout[i].hillslope.ehout = -CREHT_nosqrt(rout[i].hillslope.hslpsqrt, rout[i].hillslope.nh, 
                            rout[i].drainage_density, rout[i].hillslope.yh);
    if (rout[i].hillslope.ehout < 0. &&
        rout[i].hillslope.wh + (rout[i].runoff + rout[i].hillslope.ehout) * step_dt < TOL_VALUE) {
        rout[i].hillslope.ehout = -(rout[i].runoff + rout[i].hillslope.wh / step_dt);
    }
    rout[i].hillslope.dwh = (rout[i].runoff + rout[i].hillslope.ehout);
    rout[i].hillslope.wh += rout[i].hillslope.dwh * step_dt;

    rout[i].sub_channel.etin = (-rout[i].hillslope.ehout + rout[i].baseflow) * 
                    local_domain.locations[i].area * local_domain.locations[i].frac;
    // Update hydraulic properties of the hillslope
    rout[i].hillslope.yh = rout[i].hillslope.wh;

    /********************************
      Sub network routing
    ********************************/ 
    rout[i].main_channel.erlateral = 0.0;
    double discharge_volume = 0.0;
    for (j = 0; j < DLevelH2R; j++) {
        double erlateral_sum = 0.0;   // 用于累计侧向入流体积
        localDeltaT = step_dt / DLevelH2R / sub_steps;
        for (k = 0; k < sub_steps; k++) {
            if (rout[i].sub_channel.tlen <= rout[i].hillslope.hlen) {
                rout[i].sub_channel.etout = -rout[i].sub_channel.etin;
            }
            else {
                rout[i].sub_channel.vt = CRVRMAN_nosqrt(rout[i].sub_channel.tslpsqrt,
                            rout[i].sub_channel.nt, rout[i].sub_channel.rt);
                rout[i].sub_channel.etout = -rout[i].sub_channel.vt * rout[i].sub_channel.mt;
                if (rout[i].sub_channel.wt + (rout[i].sub_channel.etin + 
                                        rout[i].sub_channel.etout) * localDeltaT < TOL_VALUE) {
                    rout[i].sub_channel.etout = -(rout[i].sub_channel.etin +
                                                rout[i].sub_channel.wt / localDeltaT);
                    if (rout[i].sub_channel.mt > 0.0) {
                        rout[i].sub_channel.vt = -rout[i].sub_channel.etout/rout[i].sub_channel.mt;
                    }
                }
            }
            rout[i].sub_channel.dwt = rout[i].sub_channel.etout + rout[i].sub_channel.etin;
            rout[i].sub_channel.wt += rout[i].sub_channel.dwt * localDeltaT;
            // Update hydraulic properties of the channel
            if (rout[i].sub_channel.tlen > 0.0 && rout[i].sub_channel.wt > 0.0) {
                rout[i].sub_channel.mt = GRMR(rout[i].sub_channel.wt, rout[i].sub_channel.tlen);    // 过水断面面积
                rout[i].sub_channel.yt = GRHT(rout[i].sub_channel.mt, rout[i].sub_channel.twidth);  // 水深
                rout[i].sub_channel.pt = GRPT(rout[i].sub_channel.yt, rout[i].sub_channel.twidth);  // 湿周
                rout[i].sub_channel.rt = GRRR(rout[i].sub_channel.mt, rout[i].sub_channel.pt);      // 水力半径
            }
            else {
                rout[i].sub_channel.mt = 0.0;
                rout[i].sub_channel.yt = 0.0;
                rout[i].sub_channel.pt = 0.0;
                rout[i].sub_channel.rt = 0.0;
            }
            erlateral_sum += (-rout[i].sub_channel.etout) * localDeltaT;
        }
        rout[i].main_channel.erlateral = erlateral_sum / (step_dt / DLevelH2R);   // m³/s

        /********************************
             Main network routing
        ********************************/
        localDeltaT = step_dt / DLevelH2R / river_steps;
        for (m = 0; m < river_steps; m++) {
            // 上游来水
            rout[i].main_channel.erin = rout[i].upstream;
            // 无长度河道：所有水立即出流
            if (rout[i].main_channel.rlen <= 0.0) {
                rout[i].main_channel.vr = 0.0;
                rout[i].main_channel.erout = -(rout[i].main_channel.erin + 
                                rout[i].main_channel.erlateral); // 上游来水 + 侧向入流
            }
            else {
                if (rout[i].acc_area / rout[i].main_channel.rwidth / rout[i].main_channel.rlen > param.MAX_LIMIT) {
                    rout[i].main_channel.erout = -rout[i].main_channel.erin - rout[i].main_channel.erlateral;
                }
                else {
                    rout[i].main_channel.vr = CRVRMAN_nosqrt(rout[i].main_channel.rslpsqrt,
                                rout[i].main_channel.nr, rout[i].main_channel.rr);
                    rout[i].main_channel.erout = -rout[i].main_channel.vr * rout[i].main_channel.mr;
                    if (rout[i].main_channel.erout < 0.0 && rout[i].main_channel.wr +
                        (rout[i].main_channel.erin + rout[i].main_channel.erlateral + 
                            rout[i].main_channel.erout) * localDeltaT < TOL_VALUE) {
                        rout[i].main_channel.erout = -(rout[i].main_channel.erin + rout[i].main_channel.wr +
                            rout[i].main_channel.erlateral / localDeltaT);
                        if (rout[i].main_channel.mr > 0.0) {
                            rout[i].main_channel.vr = -rout[i].main_channel.erout/rout[i].main_channel.mr;
                        }
                    }
                }
            }
            rout[i].main_channel.dwr = rout[i].main_channel.erin + rout[i].main_channel.erout + 
                    rout[i].main_channel.erlateral;
            if (rout[i].main_channel.wr / localDeltaT + rout[i].main_channel.dwr < -TOL_VALUE) {
                log_warn("Negative storage in main channel, setting outflow to zero");
                rout[i].main_channel.erout = 0.0;
                rout[i].main_channel.vr = 0.0;
                rout[i].main_channel.dwr = rout[i].main_channel.erin + rout[i].main_channel.erlateral;
            }
            rout[i].main_channel.wr += rout[i].main_channel.dwr * localDeltaT;

            // Update hydraulic properties of the channel
            if (rout[i].main_channel.rlen > 0.0 && rout[i].main_channel.wr > 0.0) {
                rout[i].main_channel.mr = GRMR(rout[i].main_channel.wr, rout[i].main_channel.rlen);
                rout[i].main_channel.yr = GRHR(rout[i].main_channel.mr, rout[i].main_channel.rwidth, 
                                        rout[i].main_channel.rwidth0, rout[i].main_channel.rdepth);
                rout[i].main_channel.pr = GRPR(rout[i].main_channel.yr, rout[i].main_channel.rwidth,
                                        rout[i].main_channel.rwidth0, rout[i].main_channel.rdepth);
                rout[i].main_channel.rr = GRRR(rout[i].main_channel.mr, rout[i].main_channel.pr);
            }
            else {
                rout[i].main_channel.mr = 0.0;
                rout[i].main_channel.yr = 0.0;
                rout[i].main_channel.pr = 0.0;
                rout[i].main_channel.rr = 0.0;
            }
            //temp_erout += rout.main_channel.erout[igrid];
            discharge_volume += -rout[i].main_channel.erout * localDeltaT;
        }   
    }
    // 更新时段内平均discharge
    rout[i].discharge = discharge_volume / step_dt;

    return (0);
}

/******************************************************************************
 * @brief    Calculating channel velocity according to Manning's equation.
 *****************************************************************************/
double
CRVRMAN_nosqrt(double sqrtslp,
               double nn,
               double rr)
{
    double vr;

    if (rr < 0.0) {
        vr = 0.0;
    }
    else {
        vr = pow((rr*rr), (1.0/3.0)) * sqrtslp / nn;
    }
    return vr;
}

/******************************************************************************
 * @brief   Function for overland from hillslope into the sub-network channels
 *****************************************************************************/
double
CREHT_nosqrt(double sqrthslp,
             double nh,
             double Gxr,
             double yh)
{
    double eht;
    double vh;

    vh = CRVRMAN_nosqrt(sqrthslp, nh, yh);
    eht = Gxr * yh * vh;

    return eht;
}

/******************************************************************************
 * @brief   Function for estimate wetted channel area
 *****************************************************************************/
double
GRMR(double wr,     // storage of water in the channel (mm)
     double rlen)   // channel length
{
    double mr;

    if (rlen > 0.0) {
        mr = wr / rlen;
    }
    else {
        mr = 0.0;
    }

    return mr;
}

/******************************************************************************
 * @brief   Function for estimating water depth assuming rectangular channel
 *****************************************************************************/
double
GRHT(double mt,     // wetted channel area
     double twid)   // channel width
{
    double ht;

    if (mt < TOL_VALUE) {
        ht = 0.0;
    }
    else {
        ht = mt / twid;
    }
    return ht;
}

/******************************************************************************
 * @brief   Function for estimating wetted perimeter assuming rectangular channel
 *****************************************************************************/
double
GRPT(double ht,     // wetted channel depth
     double twid)   // channel width
{
    double pr;      // wetted perimeter

    if (ht < TOL_VALUE) {
        pr = 0.0;
    }
    else {
        pr = 2.0 * ht + twid;
    }
    return pr;
}

/******************************************************************************
 * @brief   Function for estimating wetted perimeter assuming rectangular channel
 *****************************************************************************/
double
GRRR(double mr,     // wetted channel area
     double pr)     // wetted perimeter
{
    double rr;      // hydraulic radius

    if (pr < TOL_VALUE) {
        rr = 0.0;
    }
    else {
        rr = mr / pr;
    }
    return rr;
}

/******************************************************************************
 * @brief   Function for estimating maximum water depth assuming rectangular 
 *          channel and tropezoidal flood plain.
 *****************************************************************************/
double
GRHR(double mr,       // wetted channel area
     double rwidth,   // channel width
     double rwidth0,  // channel width at bankfull
     double rdepth)   // channel length
{
    double hr;      // water depth
    double slope = 0.1;
    double deltamr = 0.0;
    if (mr <= TOL_VALUE) {
        hr = 0.0;
    }
    else {
        if (mr - rwidth * rwidth < TOL_VALUE) {
            hr = mr / rwidth;
        }
        else {
            if (mr > rdepth * rwidth + (rwidth + rwidth0) * slope * ((rwidth0 - rwidth) / 2.0) / 2.0 + TOL_VALUE) {
                deltamr = mr - rdepth * rwidth - (rwidth + rwidth0) * slope * ((rwidth0 - rwidth) / 2.0) / 2.0;
                hr = rdepth + slope * ((rwidth0 - rwidth) / 2.0) + deltamr / (rwidth0);
            }
            else {
               deltamr = mr - rdepth * rwidth;
               hr = rdepth + (-rwidth + sqrt((rwidth * rwidth) + 4.0 * deltamr / slope)) * slope / 2.0;
            }
        }
    }
    return (hr);

}

/******************************************************************************
 * @brief   Function for estimating maximum water depth assuming rectangular 
 *          channel and tropezoidal flood plain.
 *****************************************************************************/
double
GRPR(double hr,       // wetted channel depth
     double rwidth,   // channel width
     double rwidth0,  // channel width at bankfull
     double rdepth)   // channel length
{
    double slope = 0.1;
    double deltahr = 0.0;
    double pr;      // wetted perimeter
    double sinatanSLOPE1defr = 0.0;

    bool first = true;
    if (first) {
        sinatanSLOPE1defr = 1.0 / (sin(atan(slope)));
        first = false;
    }
    if(hr < TOL_VALUE) {
        pr = 0.0;
    }
    else {
        if(hr <= rdepth + TOL_VALUE) { // not flooded
            pr = 2.0 * hr + rwidth;
        }
        else {
            if(hr > rdepth + ((rwidth0-rwidth)/2.0) * slope + TOL_VALUE) { // fully flooded
                deltahr = hr - rdepth - ((rwidth0-rwidth)/2.0)*slope;
                pr =  rwidth + 2.0*(rdepth + ((rwidth0-rwidth)/2.0) * slope * sinatanSLOPE1defr + deltahr);
            }
            else {
                pr = rwidth + 2.0*(rdepth + (hr - rdepth)*sinatanSLOPE1defr);
            }
        }
    }
    return pr;
}
