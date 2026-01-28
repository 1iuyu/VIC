/******************************************************************************
* @section DESCRIPTION
*
* Calculate soil thermal conduction.
******************************************************************************/

#include <vic_run.h>

/******************************************************************************
* @brief    Soil thermal conductivity calculated using Johansen's method.
*
* @note     Reference: Farouki, O.T., "Thermal Properties of Soils" 1986
*               Chapter 7: Methods for Calculating the Thermal Conductivity
*               of Soils
******************************************************************************/
double
soil_conductivity(double moist,
                  double liq,
                  double Wsat_node,
                  double quartz,
                  double organic)
{
    double Ke;
    double Ksat;
    double Kdry;        /* Dry thermal conductivity of soil (W/mK), including mineral and organic fractions */
    double Kdry_org = 0.05; /* Dry thermal conductivity of organic fraction (W/mK) (Farouki 1981) */
    double Kdry_min;    /* Dry thermal conductivity of mineral fraction (W/mK) */
    double Ks;          /* thermal conductivity of solid (W/mK), including mineral and organic fractions */
    double Ks_org = 0.25; /* thermal conductivity of organic fraction of solid (W/mK) (Farouki 1981) */
    double Ks_min;      /* thermal conductivity of mineral fraction of solid (W/mK) */
    double Sr;          /* fractional degree of saturation */
    double K;

    /* Calculate dry conductivity as weighted average of mineral and organic fractions. */
    double bulk_dens = (1.0 - Wsat_node) * 2700.;
    Kdry_min =
        (0.135 * bulk_dens +
         64.7) / (2700. - 0.947 * bulk_dens);
    Kdry = (1 - organic) * Kdry_min + organic * Kdry_org;

    if (moist > 0.) {

        Sr = moist / Wsat_node;

        // Compute Ks of mineral soil; here "quartz" is the fraction (quartz volume / mineral soil volume)
        if (quartz < 0.2) {
            Ks_min = pow(CONST_KQUARTZ, quartz) * pow(3.0, 1.0 - quartz); // when quartz is less than 0.2
        }
        else {
            Ks_min = pow(CONST_KQUARTZ, quartz) * pow(2.2, 1.0 - quartz); // when quartz is greater than 0.2
        }
        Ks = (1 - organic) * Ks_min + organic * Ks_org;
        Ksat =
            pow(Ks, 1.0 - Wsat_node) * pow(CONST_KICE, Wsat_node - liq) * pow(CONST_KFWICE, liq);

        if (liq + 0.0005 < moist) {
            /** Soil frozen **/
            Ke = Sr;

        }
        else {
            /** Soil unfrozen **/
            if (Sr > 0.1) {
                Ke = 0.7 * log10(Sr) + 1.0;               
            }
            else {
                Ke = 0.0;
            }
        }

        K = (Ksat - Kdry) * Ke + Kdry;
        if (K < Kdry) {
            K = Kdry;
        }
    }
    else {
        K = Kdry;
    }

    return (K);
}

/******************************************************************************
* @brief    This subroutine calculates the soil volumetric heat capacity
            based on the fractional volume of its component parts.
******************************************************************************/
double
volumetric_heat_capacity(double Wsat_node,
                         double liq,
                         double ice,
                         double moist,
                         double organic)
{
    double Cs;

    // Constant values are volumetric heat capacities in J/m^3/K
    Cs =  CONST_CPFWICE * liq; // water
    Cs += CONST_CPSOIL * (1. - Wsat_node) * (1 - organic); // soil
    Cs += CONST_CPDAIR * (Wsat_node - moist); // air
    Cs += CONST_CPICE * ice; // ice
    Cs += CONST_CPORGANIC * (1. - Wsat_node) * organic; // organic

    return (Cs);
}

/******************************************************************************
* @brief    This subroutine sets the thermal node soil parameters to constant
*           values based on those defined for the current grid cells soil type.
*           Thermal node propertiers for the energy balance solution are also
*           set (these constants are used to reduce the solution time required
*           within each iteration).
******************************************************************************/
void
set_node_parameters(soil_con_struct *soil_con)
{
    extern option_struct options;

    char    PAST_BOTTOM;
    size_t  nidx, lidx;
    double  Lsum; /* cumulative depth of moisture layer */
    // layer 
    double *bexp = soil_con->bexp;
    double *bubble = soil_con->bubble;
    double *Dsat = soil_con->Dsat;
    double *expt = soil_con->expt;
    double *Ksat = soil_con->Ksat; 
    double *Wfc = soil_con->Wfc;
    double *Wpwp = soil_con->Wpwp;
    double *Wsat = soil_con->Wsat;
    double *depth = soil_con->depth;
    double *psi_sat = soil_con->psi_sat;
    double *organic = soil_con->organic;
    double *quartz = soil_con->quartz;
    double *resid_moist = soil_con->resid_moist;
    double *soil_dens_min = soil_con->soil_dens_min;
    double *bulk_dens_min = soil_con->bulk_dens_min;
    // node
    double *bexp_node = soil_con->bexp_node;
    double *bubble_node = soil_con->bubble_node;
    double *Dsat_node = soil_con->Dsat_node;
    double *expt_node = soil_con->expt_node;
    double *Ksat_node = soil_con->Ksat_node;
    double *Wfc_node = soil_con->Wfc_node;
    double *Wsat_node = soil_con->Wsat_node;
    double *Wpwp_node = soil_con->Wpwp_node;
    double *Zsum_soil = soil_con->Zsum_soil;
    double *psi_sat_node = soil_con->psi_sat_node;
    double *organic_node = soil_con->organic_node;
    double *quartz_node = soil_con->quartz_node;
    double *resid_moist_node = soil_con->resid_moist_node;
    double *soil_dens_min_node = soil_con->soil_dens_min_node;
    double *bulk_dens_min_node = soil_con->bulk_dens_min_node;


    PAST_BOTTOM = false;
    lidx = 0;
    Lsum = 0.;
    size_t Nnode = options.Nnode;
    size_t Nlayer = options.Nlayer;

    /* set node parameters */
    for (nidx = 0; nidx < Nnode; nidx++) {
        if (Zsum_soil[nidx] == Lsum + depth[lidx] && nidx != 0 && lidx !=
            Nlayer - 1) {
            /* node on layer boundary */
            expt_node[nidx] = (expt[lidx] + expt[lidx + 1]) / 2.;
            bubble_node[nidx] = (bubble[lidx] + bubble[lidx + 1]) / 2.;
            resid_moist_node[nidx] = (resid_moist[lidx] + resid_moist[lidx + 1]) / 2.;
            Ksat_node[nidx] = (Ksat[lidx] + Ksat[lidx + 1]) / 2.;
            Dsat_node[nidx] = (Dsat[lidx] + Dsat[lidx + 1]) / 2.;
            Wfc_node[nidx] = (Wfc[lidx] + Wfc[lidx + 1]) / 2.;
            Wpwp_node[nidx] = (Wpwp[lidx] + Wpwp[lidx + 1]) / 2.;
            Wsat_node[nidx] = (Wsat[lidx] + Wsat[lidx + 1]) / 2.;
            psi_sat_node[nidx] = (psi_sat[lidx] + psi_sat[lidx + 1]) / 2.;
            bexp_node[nidx] = (bexp[lidx] + bexp[lidx + 1]) / 2.;
            organic_node[nidx] = (organic[lidx] + organic[lidx + 1]) / 2.;
            quartz_node[nidx] = (quartz[lidx] + quartz[lidx + 1]) / 2.;
            soil_dens_min_node[nidx] = (soil_dens_min[lidx] + soil_dens_min[lidx + 1]) / 2.;
            bulk_dens_min_node[nidx] = (bulk_dens_min[lidx] + bulk_dens_min[lidx + 1]) / 2.;
        }
        else {
            /* node completely in layer */
            expt_node[nidx] = expt[lidx];
            bubble_node[nidx] = bubble[lidx];
            resid_moist_node[nidx] = resid_moist[lidx];
            Ksat_node[nidx] = Ksat[lidx];
            Dsat_node[nidx] = Dsat[lidx];
            Wfc_node[nidx] = Wfc[lidx];
            Wpwp_node[nidx] = Wpwp[lidx];
            Wsat_node[nidx] = Wsat[lidx];
            psi_sat_node[nidx] = psi_sat[lidx];
            bexp_node[nidx] = bexp[lidx];
            organic_node[nidx] = organic[lidx];
            quartz_node[nidx] = quartz[lidx];
            soil_dens_min_node[nidx] = soil_dens_min[lidx];
            bulk_dens_min_node[nidx] = bulk_dens_min[lidx];
        }
        if (Zsum_soil[nidx] > Lsum + depth[lidx] && !PAST_BOTTOM) {
            Lsum += depth[lidx];
            lidx++;
            if (lidx == Nlayer) {
                PAST_BOTTOM = true;
                lidx = Nlayer - 1;
            }
        }
    }
}

/******************************************************************************
* @brief    This subroutine determines the moisture and ice contents of each
*           soil thermal node based on the current node temperature and layer
*           moisture content.  Thermal conductivity and volumetric heat
*           capacity are then estimated for each node based on the division of
*           moisture contents.
******************************************************************************/
int
distribute_node_moisture_properties(bool    IS_GLAC,
                                    double *moist,
                                    double *ice,
                                    double *liq,
                                    double *soil_T,
                                    double *Wsat_node,
                                    double *expt_node,
                                    double *bubble_node,
                                    char    FS_ACTIVE)
{
    extern option_struct     options;

    size_t    nidx;

    /* node estimates */
    for (nidx = 0; nidx < options.Nnode; nidx++) {
        if (IS_GLAC) {
            ice[nidx] = 1.0;
            moist[nidx] = 1.0;
        }
        else {
            if (soil_T[nidx] < CONST_TKFRZ && (FS_ACTIVE && options.FROZEN_SOIL)) {
                /* compute moisture and ice contents */
                ice[nidx] =
                    moist[nidx] - maximum_unfrozen_water(soil_T[nidx],
                                                         Wsat_node[nidx],
                                                         bubble_node[nidx],
                                                         expt_node[nidx]);
                if (ice[nidx] < 0) {
                    ice[nidx] = 0;
                }
            }
            liq[nidx] = moist[nidx] - ice[nidx];          
        }
    }
    return (0);
}

/******************************************************************************
* @brief    This subroutine reads through the soil thermal nodes and determines
*           the depths of all thawing and freezing fronts that are present.
******************************************************************************/
void
find_0_degree_fronts(energy_bal_struct *energy,
                     double            *Zsum_node,
                     double            *T,
                     int                Nnode)
{
    int    nidx, fidx;
    int    Nthaw; /* number of thawing fronts found */
    int    Nfrost; /* number of frost fronts found */
    double tdepth[MAX_FRONTS]; /* thawing frost depths */
    double fdepth[MAX_FRONTS]; /* freezing front depths */

    /* Initialize parameters */
    Nthaw = Nfrost = 0;
    for (fidx = 0; fidx < MAX_FRONTS; fidx++) {
        fdepth[fidx] = MISSING;
        tdepth[fidx] = MISSING;
    }

    /* find 0 degree fronts */
    for (nidx = Nnode - 2; nidx >= 0; nidx--) {
        if (T[nidx] > CONST_TKFRZ && T[nidx + 1] <= CONST_TKFRZ && Nthaw < MAX_FRONTS) {
            tdepth[Nthaw] =
                linear_interp(0, T[nidx], T[nidx + 1], Zsum_node[nidx],
                              Zsum_node[nidx + 1]);
            Nthaw++;
        }
        else if (T[nidx] < CONST_TKFRZ && T[nidx + 1] >= CONST_TKFRZ && Nfrost < MAX_FRONTS) {
            fdepth[Nfrost] =
                linear_interp(0, T[nidx], T[nidx + 1], Zsum_node[nidx],
                              Zsum_node[nidx + 1]);
            Nfrost++;
        }
    }

    /* store thaw depths */
    for (fidx = 0; fidx < MAX_FRONTS; fidx++) {
        energy->tdepth[fidx] = tdepth[fidx];
    }
    /* store frost depths */
    for (fidx = 0; fidx < MAX_FRONTS; fidx++) {
        energy->fdepth[fidx] = fdepth[fidx];
    }
    energy->Nthaw = Nthaw;
    energy->Nfrost = Nfrost;
}

/******************************************************************************
* @brief    This subroutine computes the maximum amount of unfrozen water that
*           can exist at the current temperature.
******************************************************************************/
double
maximum_unfrozen_water(double T,
                       double max_moist,
                       double bubble,
                       double expt)
{
    double unfrozen;

    if (T < CONST_TKFRZ) {
        unfrozen = max_moist *
                   pow((-CONST_LATICE *
                        (T - CONST_TKFRZ)) / (CONST_TKTRIP) / (CONST_G * bubble / (CM_PER_M)),
                       -(2.0 / (expt - 3.0)));
        if (unfrozen > max_moist) {
            unfrozen = max_moist;
        }
        if (unfrozen < 0.) {
            unfrozen = 0.;
        }
    }
    else {
        unfrozen = max_moist;
    }

    return (unfrozen);
}
