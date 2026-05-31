/******************************************************************************
 * @section DESCRIPTION
 *
 * Allocate memory for Routing structures.
 *****************************************************************************/

#include "rout.h"

/******************************************************************************
 * @brief    Allocate memory for Routing structures.
 *****************************************************************************/
void
rout_alloc(void)
{
    extern int mpi_rank;
    if (mpi_rank == VIC_MPI_ROOT) {
        extern domain_struct    global_domain;
        extern rout_struct      rout;

        // Allocate memory in rout param_struct
        rout.rout_param.downstream = malloc(
            global_domain.ncells_active * sizeof(*rout.rout_param.downstream));
        check_alloc_status(rout.rout_param.downstream,
                           "Memory allocation error.");

        rout.rout_param.indegree = malloc(
            global_domain.ncells_active * sizeof(*rout.rout_param.indegree));
        check_alloc_status(rout.rout_param.indegree,
                           "Memory allocation error.");

        rout.rout_param.queue = malloc(
            global_domain.ncells_active * sizeof(*rout.rout_param.queue));
        check_alloc_status(rout.rout_param.queue,
                           "Memory allocation error.");

        rout.rout_param.routing_order = malloc(
            global_domain.ncells_active * sizeof(*rout.rout_param.routing_order));
        check_alloc_status(rout.rout_param.routing_order,
                           "Memory allocation error.");

        rout.rout_param.source_torow = malloc(
            global_domain.ncells_active * sizeof(*rout.rout_param.source_torow));
        check_alloc_status(rout.rout_param.source_torow,
                           "Memory allocation error.");

        rout.rout_param.source_tocol = malloc(
            global_domain.ncells_active * sizeof(*rout.rout_param.source_tocol));
        check_alloc_status(rout.rout_param.source_tocol,
                           "Memory allocation error.");

        rout.discharge =
            malloc(global_domain.ncells_active * sizeof(*rout.discharge));
        check_alloc_status(rout.discharge, "Memory allocation error.");

        rout.runoff =
            malloc(global_domain.ncells_active * sizeof(*rout.runoff));
        check_alloc_status(rout.runoff, "Memory allocation error.");

        rout.baseflow =
            malloc(global_domain.ncells_active * sizeof(*rout.baseflow));
        check_alloc_status(rout.baseflow, "Memory allocation error.");

        rout.drainage_density =
            malloc(global_domain.ncells_active * sizeof(*rout.drainage_density));
        check_alloc_status(rout.drainage_density, "Memory allocation error.");

        rout.total_length =
            malloc(global_domain.ncells_active * sizeof(*rout.total_length));
        check_alloc_status(rout.total_length, "Memory allocation error.");

        rout.acc_area =
            malloc(global_domain.ncells_active * sizeof(*rout.acc_area));
        check_alloc_status(rout.acc_area, "Memory allocation error.");

        rout.river_steps = 
            malloc(global_domain.ncells_active * sizeof(*rout.river_steps));
        check_alloc_status(rout.river_steps, "Memory allocation error.");

        rout.sub_steps = 
            malloc(global_domain.ncells_active * sizeof(*rout.sub_steps));
        check_alloc_status(rout.sub_steps, "Memory allocation error.");

        rout.upstream = 
            malloc(global_domain.ncells_active * sizeof(*rout.upstream));
        check_alloc_status(rout.upstream, "Memory allocation error.");

        // 初始化main_channel
        rout.main_channel.dwr = malloc(
            global_domain.ncells_active * sizeof(*rout.main_channel.dwr));
        check_alloc_status(rout.main_channel.dwr,
                           "Memory allocation error.");
        rout.main_channel.erin = malloc(
            global_domain.ncells_active * sizeof(*rout.main_channel.erin));
        check_alloc_status(rout.main_channel.erin,
                           "Memory allocation error.");
        rout.main_channel.erlateral = malloc(
            global_domain.ncells_active * sizeof(*rout.main_channel.erlateral));
        check_alloc_status(rout.main_channel.erlateral,
                           "Memory allocation error.");
        rout.main_channel.erout = malloc(
            global_domain.ncells_active * sizeof(*rout.main_channel.erout));
        check_alloc_status(rout.main_channel.erout,
                           "Memory allocation error.");
         rout.main_channel.mr = malloc(
            global_domain.ncells_active * sizeof(*rout.main_channel.mr));
        check_alloc_status(rout.main_channel.mr,
                           "Memory allocation error.");
        rout.main_channel.nr = malloc(
            global_domain.ncells_active * sizeof(*rout.main_channel.nr));
        check_alloc_status(rout.main_channel.nr,
                           "Memory allocation error.");
        rout.main_channel.pr = malloc(
            global_domain.ncells_active * sizeof(*rout.main_channel.pr));
        check_alloc_status(rout.main_channel.pr,
                           "Memory allocation error.");
        rout.main_channel.rdepth = malloc(
            global_domain.ncells_active * sizeof(*rout.main_channel.rdepth));
        check_alloc_status(rout.main_channel.rdepth,
                           "Memory allocation error.");
        rout.main_channel.rlen = malloc(
            global_domain.ncells_active * sizeof(*rout.main_channel.rlen));
        check_alloc_status(rout.main_channel.rlen,
                           "Memory allocation error.");
        rout.main_channel.rr = malloc(
            global_domain.ncells_active * sizeof(*rout.main_channel.rr));
        check_alloc_status(rout.main_channel.rr,
                           "Memory allocation error.");
        rout.main_channel.rslpsqrt = malloc(
            global_domain.ncells_active * sizeof(*rout.main_channel.rslpsqrt));
        check_alloc_status(rout.main_channel.rslpsqrt,
                           "Memory allocation error.");
        rout.main_channel.rwidth0 = malloc(
            global_domain.ncells_active * sizeof(*rout.main_channel.rwidth0));
        check_alloc_status(rout.main_channel.rwidth0,
                           "Memory allocation error.");
        rout.main_channel.rwidth = malloc(
            global_domain.ncells_active * sizeof(*rout.main_channel.rwidth));
        check_alloc_status(rout.main_channel.rwidth,
                           "Memory allocation error.");
        rout.main_channel.vr = malloc(
            global_domain.ncells_active * sizeof(*rout.main_channel.vr));
        check_alloc_status(rout.main_channel.vr,
                           "Memory allocation error.");
        rout.main_channel.wr = malloc(
            global_domain.ncells_active * sizeof(*rout.main_channel.wr));
        check_alloc_status(rout.main_channel.wr,
                           "Memory allocation error.");
        rout.main_channel.yr = malloc(
            global_domain.ncells_active * sizeof(*rout.main_channel.yr));
        check_alloc_status(rout.main_channel.yr,
                           "Memory allocation error.");
                           
        /* Memory allocation for sub_channel */
        rout.sub_channel.etin = malloc(
            global_domain.ncells_active * sizeof(*rout.sub_channel.etin));
        check_alloc_status(rout.sub_channel.etin,
                           "Memory allocation error.");
        rout.sub_channel.etout = malloc(
            global_domain.ncells_active * sizeof(*rout.sub_channel.etout));
        check_alloc_status(rout.sub_channel.etout,
                        "Memory allocation error.");

        rout.sub_channel.vt = malloc(
            global_domain.ncells_active * sizeof(*rout.sub_channel.vt));
        check_alloc_status(rout.sub_channel.vt,
                        "Memory allocation error.");

        rout.sub_channel.wt = malloc(
            global_domain.ncells_active * sizeof(*rout.sub_channel.wt));
        check_alloc_status(rout.sub_channel.wt,
                        "Memory allocation error.");

        rout.sub_channel.nt = malloc(
            global_domain.ncells_active * sizeof(*rout.sub_channel.nt));
        check_alloc_status(rout.sub_channel.nt,
                        "Memory allocation error.");

        rout.sub_channel.mt = malloc(
            global_domain.ncells_active * sizeof(*rout.sub_channel.mt));
        check_alloc_status(rout.sub_channel.mt,
                        "Memory allocation error.");

        rout.sub_channel.tslpsqrt = malloc(
            global_domain.ncells_total * sizeof(*rout.sub_channel.tslpsqrt));
        check_alloc_status(rout.sub_channel.tslpsqrt,
                        "Memory allocation error.");

        rout.sub_channel.dwt = malloc(
            global_domain.ncells_active * sizeof(*rout.sub_channel.dwt));
        check_alloc_status(rout.sub_channel.dwt,
                        "Memory allocation error.");

        rout.sub_channel.tlen = malloc(
            global_domain.ncells_active * sizeof(*rout.sub_channel.tlen));
        check_alloc_status(rout.sub_channel.tlen,
                        "Memory allocation error.");

        rout.sub_channel.twidth = malloc(
            global_domain.ncells_active * sizeof(*rout.sub_channel.twidth));
        check_alloc_status(rout.sub_channel.twidth,
                        "Memory allocation error.");

        rout.sub_channel.yt = malloc(
            global_domain.ncells_active * sizeof(*rout.sub_channel.yt));
        check_alloc_status(rout.sub_channel.yt,
                        "Memory allocation error.");

        rout.sub_channel.pt = malloc(
            global_domain.ncells_active * sizeof(*rout.sub_channel.pt));
        check_alloc_status(rout.sub_channel.pt,
                        "Memory allocation error.");

        rout.sub_channel.rt = malloc(
            global_domain.ncells_active * sizeof(*rout.sub_channel.rt));
        check_alloc_status(rout.sub_channel.rt,
                        "Memory allocation error.");

        /* Memory allocation for hillslope */
        rout.hillslope.ehout = malloc(
            global_domain.ncells_active * sizeof(*rout.hillslope.ehout));
        check_alloc_status(rout.hillslope.ehout,
                        "Memory allocation error.");

        rout.hillslope.yh = malloc(
            global_domain.ncells_active * sizeof(*rout.hillslope.yh));
        check_alloc_status(rout.hillslope.yh,
                        "Memory allocation error.");

        rout.hillslope.nh = malloc(
            global_domain.ncells_active * sizeof(*rout.hillslope.nh));
        check_alloc_status(rout.hillslope.nh,
                        "Memory allocation error.");

        rout.hillslope.wh = malloc(
            global_domain.ncells_active * sizeof(*rout.hillslope.wh));
        check_alloc_status(rout.hillslope.wh,
                        "Memory allocation error.");

        rout.hillslope.hslpsqrt = malloc(
            global_domain.ncells_total * sizeof(*rout.hillslope.hslpsqrt));
        check_alloc_status(rout.hillslope.hslpsqrt,
                        "Memory allocation error.");

        rout.hillslope.dwh = malloc(
            global_domain.ncells_active * sizeof(*rout.hillslope.dwh));
        check_alloc_status(rout.hillslope.dwh,
                        "Memory allocation error.");

        rout.hillslope.hlen = malloc(
            global_domain.ncells_active * sizeof(*rout.hillslope.hlen));
        check_alloc_status(rout.hillslope.hlen,
                        "Memory allocation error.");      
    }
}
