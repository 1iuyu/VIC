/******************************************************************************
 * @section DESCRIPTION
 *
 * clean up functions for routing extension
 *****************************************************************************/

#include "rout.h"

/******************************************************************************
 * @brief    Finalize RVIC by freeing memory.
 *****************************************************************************/
void
rout_finalize(void)
{
    extern rout_struct rout;

    free(rout.rout_param.downstream);
    free(rout.rout_param.indegree);
    free(rout.rout_param.queue);
    free(rout.rout_param.routing_order);
    free(rout.rout_param.source_torow);
    free(rout.rout_param.source_tocol);
    free(rout.discharge);
    free(rout.baseflow);
    free(rout.runoff);
    free(rout.drainage_density);
    free(rout.total_length);
    free(rout.acc_area);
    free(rout.river_steps);
    free(rout.sub_steps);
    free(rout.upstream);
    /* Free sub_channel memory */
    free(rout.sub_channel.etin);
    free(rout.sub_channel.etout);
    free(rout.sub_channel.vt);
    free(rout.sub_channel.wt);
    free(rout.sub_channel.nt);
    free(rout.sub_channel.mt);
    free(rout.sub_channel.tslpsqrt);
    free(rout.sub_channel.dwt);
    free(rout.sub_channel.tlen);
    free(rout.sub_channel.twidth);
    free(rout.sub_channel.yt);
    free(rout.sub_channel.pt);
    free(rout.sub_channel.rt);
    /* Free hillslope memory */
    free(rout.hillslope.ehout);
    free(rout.hillslope.yh);
    free(rout.hillslope.nh);
    free(rout.hillslope.wh);
    free(rout.hillslope.hslpsqrt);
    free(rout.hillslope.dwh);
    free(rout.hillslope.hlen);
    /* Free main_channel memory */
    free(rout.main_channel.erin);
    free(rout.main_channel.erout);
    free(rout.main_channel.vr);
    free(rout.main_channel.wr);
    free(rout.main_channel.mr);
    free(rout.main_channel.yr);
    free(rout.main_channel.pr);
    free(rout.main_channel.nr);
    free(rout.main_channel.rr);
    free(rout.main_channel.erlateral);
    free(rout.main_channel.rslpsqrt);
    free(rout.main_channel.dwr);
    free(rout.main_channel.rlen);
    free(rout.main_channel.rwidth);
    free(rout.main_channel.rwidth0);
    free(rout.main_channel.rdepth);
}
