/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate the aerodynamic resistances.
 *****************************************************************************/

#include "vic_run.h"

/******************************************************************************
 * @brief    Calculate the aerodynamic resistance for each vegetation layer,
 *           based on Monin-Obukhov (M-O) Similarity Theory (MOST).
 *****************************************************************************/
void
ActiveLayer(cell_data_struct *cell,
            soil_con_struct  *soil_con)
{
    size_t Nsoil = cell->Nsoil;
    int nidx;
    size_t Nthaw = 0;
    size_t Nfrost = 0;
    double *Zsum_soil = soil_con->Zsum_soil;
    double *ice = cell->ice;
    double *liq = cell->liq;
    double *moist = cell->moist;

    bool found_frost = false;

    /*************************************************
     * 从底向上扫描：找冻结层和融化层
     *************************************************/
    for (nidx = Nsoil - 1; nidx >= 0; nidx--) {

        if (ice[nidx] > 0.0) {
            if (!found_frost) {
                Nfrost = nidx;
                found_frost = true;
            }
        }
        else {
            if (found_frost && Nthaw == 0) {
                Nthaw = nidx + 1;
                break;
            }
        }
    }

    /*************************************************
     * 计算冻结深度 FDEPTH
     *************************************************/
    double fdepth = 0.0;

    if (Nfrost > 0 || (Nfrost == 0 && ice[0] > 0.0)) {

        size_t nf = Nfrost;

        double fractn = max(0.0, ice[nf] / moist[nf]);

        if (nf == 0) {
            fdepth = fractn * (Zsum_soil[1] - Zsum_soil[0]) / 2.0;
        }
        else if (nf == Nsoil - 1) {
            fdepth = (Zsum_soil[nf] + Zsum_soil[nf - 1]) / 2.0
                + fractn * (Zsum_soil[nf] - Zsum_soil[nf - 1]) / 2.0;
        }
        else {
            fdepth = (Zsum_soil[nf] + Zsum_soil[nf - 1]) / 2.0
                + fractn * (Zsum_soil[nf + 1] - Zsum_soil[nf - 1]) / 2.0;
        }
    }

    /*************************************************
     * 计算融化深度 TDEPTH
     *************************************************/
    double tdepth = 0.0;

    if (Nthaw > 0) {

        size_t nt = Nthaw;

        double fractn = max(0.0, liq[nt] / moist[nt]);;

        if (nt == 0) {
            tdepth = fractn * (Zsum_soil[1] - Zsum_soil[0]) / 2.0;
        }
        else {

            if (nt == Nfrost) {
                double mid = (Zsum_soil[nt] + Zsum_soil[nt - 1]) / 2.0;
                tdepth = mid + fractn * (fdepth - mid);
            }
            else {
                tdepth = (Zsum_soil[nt] + Zsum_soil[nt - 1]) / 2.0
                    + fractn * (Zsum_soil[nt + 1] - Zsum_soil[nt - 1]) / 2.0;
            }
        }
    }
    // 写回结构体
    cell->fdepth = fdepth;
    cell->tdepth = tdepth;
    cell->Nthaw  = Nthaw;
    cell->Nfrost = Nfrost;

}