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
Euler_Routing(size_t igrid,
              double step_dt)
{
    extern parameters_struct param;
    extern rout_struct       rout;
    extern domain_struct local_domain;

    size_t i, j, k;
    size_t river_steps = rout.river_steps[igrid];
    size_t sub_steps = rout.sub_steps[igrid];
    double localDeltaT;

    /********************************
      Hillslope routing
    ********************************/
    rout.hillslope.ehout[igrid] = -CREHT_nosqrt(rout.hillslope.hslpsqrt[igrid], rout.hillslope.nh[igrid], 
                            rout.drainage_density[igrid], rout.hillslope.yh[igrid]);
    if (rout.hillslope.ehout[igrid] < 0. &&
        rout.hillslope.wh[igrid] + (rout.runoff[igrid] + rout.hillslope.ehout[igrid]) * step_dt < TOL_VALUE) {
        rout.hillslope.ehout[igrid] = -(rout.runoff[igrid] + rout.hillslope.wh[igrid] / step_dt);
    }
    rout.hillslope.dwh[igrid] = (rout.runoff[igrid] + rout.hillslope.ehout[igrid]);
    rout.hillslope.wh[igrid] += rout.hillslope.dwh[igrid] * step_dt;

    rout.sub_channel.etin[igrid] = (-rout.hillslope.ehout[igrid] + rout.baseflow[igrid]) * 
                    local_domain.locations[igrid].area * local_domain.locations[igrid].frac;
    // Update hydraulic properties of the hillslope
    rout.hillslope.yh[igrid] = rout.hillslope.wh[igrid];

    /********************************
      Sub network routing
    ********************************/ 
    rout.main_channel.erlateral[igrid] = 0.0;
    double discharge_volume = 0.0;
    for (i = 0; i < DLevelH2R; i++) {
        double erlateral_sum = 0.0;   // 用于累计侧向入流体积
        localDeltaT = step_dt / DLevelH2R / sub_steps;
        for (j = 0; j < sub_steps; j++) {
            if (rout.sub_channel.tlen[igrid] <= rout.hillslope.hlen[igrid]) {
                rout.sub_channel.etout[igrid] = -rout.sub_channel.etin[igrid];
            }
            else {
                rout.sub_channel.vt[igrid] = CRVRMAN_nosqrt(rout.sub_channel.tslpsqrt[igrid],
                            rout.sub_channel.nt[igrid], rout.sub_channel.rt[igrid]);
                rout.sub_channel.etout[igrid] = -rout.sub_channel.vt[igrid] * rout.sub_channel.mt[igrid];
                if (rout.sub_channel.wt[igrid] + (rout.sub_channel.etin[igrid] + 
                                        rout.sub_channel.etout[igrid]) * localDeltaT < TOL_VALUE) {
                    rout.sub_channel.etout[igrid] = -(rout.sub_channel.etin[igrid] +
                                                rout.sub_channel.wt[igrid] / localDeltaT);
                    if (rout.sub_channel.mt[igrid] > 0.0) {
                        rout.sub_channel.vt[igrid] = -rout.sub_channel.etout[igrid]/rout.sub_channel.mt[igrid];
                    }
                }
            }
            rout.sub_channel.dwt[igrid] = rout.sub_channel.etout[igrid] + rout.sub_channel.etin[igrid];
            rout.sub_channel.wt[igrid] += rout.sub_channel.dwt[igrid] * localDeltaT;
            // Update hydraulic properties of the channel
            if (rout.sub_channel.tlen[igrid] > 0.0 && rout.sub_channel.wt[igrid] > 0.0) {
                rout.sub_channel.mt[igrid] = GRMR(rout.sub_channel.wt[igrid], rout.sub_channel.tlen[igrid]);    // 过水断面面积
                rout.sub_channel.yt[igrid] = GRHT(rout.sub_channel.mt[igrid], rout.sub_channel.twidth[igrid]);  // 水深
                rout.sub_channel.pt[igrid] = GRPT(rout.sub_channel.yt[igrid], rout.sub_channel.twidth[igrid]);  // 湿周
                rout.sub_channel.rt[igrid] = GRRR(rout.sub_channel.mt[igrid], rout.sub_channel.pt[igrid]);      // 水力半径
            }
            else {
                rout.sub_channel.mt[igrid] = 0.0;
                rout.sub_channel.yt[igrid] = 0.0;
                rout.sub_channel.pt[igrid] = 0.0;
                rout.sub_channel.rt[igrid] = 0.0;
            }
            erlateral_sum += (-rout.sub_channel.etout[igrid]) * localDeltaT;
        }
        rout.main_channel.erlateral[igrid] = erlateral_sum / (step_dt / DLevelH2R);   // m³/s

        /********************************
             Main network routing
        ********************************/
        localDeltaT = step_dt / DLevelH2R / river_steps;
        for (k = 0; k < river_steps; k++) {
            // 上游来水
            rout.main_channel.erin[igrid] = rout.upstream[igrid];
            // 无长度河道：所有水立即出流
            if (rout.main_channel.rlen[igrid] <= 0.0) {
                rout.main_channel.vr[igrid] = 0.0;
                rout.main_channel.erout[igrid] = -(rout.main_channel.erin[igrid] + 
                                rout.main_channel.erlateral[igrid]); // 上游来水 + 侧向入流
            }
            else {
                if (rout.acc_area[igrid] / rout.main_channel.rwidth[igrid] / rout.main_channel.rlen[igrid] > param.MAX_LIMIT) {
                    rout.main_channel.erout[igrid] = -rout.main_channel.erin[igrid] - rout.main_channel.erlateral[igrid];
                }
                else {
                    rout.main_channel.vr[igrid] = CRVRMAN_nosqrt(rout.main_channel.rslpsqrt[igrid],
                                rout.main_channel.nr[igrid], rout.main_channel.rr[igrid]);
                    rout.main_channel.erout[igrid] = -rout.main_channel.vr[igrid] * rout.main_channel.mr[igrid];
                    if (rout.main_channel.erout[igrid] < 0.0 && rout.main_channel.wr[igrid] +
                        (rout.main_channel.erin[igrid] + rout.main_channel.erlateral[igrid] + 
                            rout.main_channel.erout[igrid]) * localDeltaT < TOL_VALUE) {
                        rout.main_channel.erout[igrid] = -(rout.main_channel.erin[igrid] + rout.main_channel.wr[igrid] +
                            rout.main_channel.erlateral[igrid] / localDeltaT);
                        if (rout.main_channel.mr[igrid] > 0.0) {
                            rout.main_channel.vr[igrid] = -rout.main_channel.erout[igrid]/rout.main_channel.mr[igrid];
                        }
                    }
                }
            }
            rout.main_channel.dwr[igrid] = rout.main_channel.erin[igrid] + rout.main_channel.erout[igrid] + 
                    rout.main_channel.erlateral[igrid];
            if (rout.main_channel.wr[igrid] / localDeltaT + rout.main_channel.dwr[igrid] < -TOL_VALUE) {
                log_warn("Negative storage in main channel, setting outflow to zero");
                rout.main_channel.erout[igrid] = 0.0;
                rout.main_channel.vr[igrid] = 0.0;
                rout.main_channel.dwr[igrid] = rout.main_channel.erin[igrid] + rout.main_channel.erlateral[igrid];
            }
            rout.main_channel.wr[igrid] += rout.main_channel.dwr[igrid] * localDeltaT;

            // Update hydraulic properties of the channel
            if (rout.main_channel.rlen[igrid] > 0.0 && rout.main_channel.wr[igrid] > 0.0) {
                rout.main_channel.mr[igrid] = GRMR(rout.main_channel.wr[igrid], rout.main_channel.rlen[igrid]);
                rout.main_channel.yr[igrid] = GRHR(rout.main_channel.mr[igrid], rout.main_channel.rwidth[igrid], 
                                        rout.main_channel.rwidth0[igrid], rout.main_channel.rdepth[igrid]);
                rout.main_channel.pr[igrid] = GRPR(rout.main_channel.yr[igrid], rout.main_channel.rwidth[igrid],
                                        rout.main_channel.rwidth0[igrid], rout.main_channel.rdepth[igrid]);
                rout.main_channel.rr[igrid] = GRRR(rout.main_channel.mr[igrid], rout.main_channel.pr[igrid]);
            }
            else {
                rout.main_channel.mr[igrid] = 0.0;
                rout.main_channel.yr[igrid] = 0.0;
                rout.main_channel.pr[igrid] = 0.0;
                rout.main_channel.rr[igrid] = 0.0;
            }
            //temp_erout += rout.main_channel.erout[igrid];
            discharge_volume += -rout.main_channel.erout[igrid] * localDeltaT;
        }   
    }
    // 更新时段内平均discharge
    rout.discharge[igrid] = discharge_volume / step_dt;

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
