/******************************************************************************
 * @section DESCRIPTION #include <vic_def.h>
 *
 * Header file for vic_run routines
 *****************************************************************************/

#ifndef VIC_RUN_H
#define VIC_RUN_H

#include "vic_def.h"

void ActiveLayer(cell_data_struct *, soil_con_struct *);
void AdvectedEnergy(double, double, double, double, double,
                    double, double, double, double, double,
                    double *, double *, double *);
void AdvectedEnergyGlac(double, double, double, double, double *);
bool assert_close_double(double x, double y, double rtol, double abs_tol);
bool assert_close_float(float x, float y, float rtol, float abs_tol);
void brent_PHS(double, double, double, double, double, double, double, double, 
               double, double *, double, double, double, double, double, double,
               double, double, double, double, double, double, double, double,
               double, double, double, double, double, double, double, double,
               double *, double *, double *, double *, double *, double *,
               cell_data_struct *, soil_con_struct *, veg_var_struct *, veg_lib_struct *);
int calc_stress(double *, double *, double *, double, double, double, double,
                double, double, double, double, cell_data_struct *, 
                soil_con_struct *, veg_var_struct *, veg_lib_struct *);
int calc_energy_bal(size_t, double, double, double, force_data_struct *, 
                    energy_bal_struct *, cell_data_struct *, snow_data_struct *, 
                    soil_con_struct *, veg_var_struct *, veg_lib_struct *);
int calc_energy_bal_glac(size_t, double, double, double,
                         force_data_struct *, energy_bal_struct *, cell_data_struct *,
                         snow_data_struct *, soil_con_struct *);
int calc_water_bal(double, double, energy_bal_struct *, 
                   cell_data_struct *, soil_con_struct *);
int calc_water_bal_glac(size_t, double, double, double, double, double, double,
                        force_data_struct *, energy_bal_struct *, cell_data_struct *, 
                        snow_data_struct *, soil_con_struct *);
void calc_rainonly(double, double, double, double, double, double, double *, double *);
double calc_veg_displacement(double, double, double);
void calc_root_moist_stress(cell_data_struct *, soil_con_struct *, veg_lib_struct *);
void calc_snow_coverage(double, bool, double,
                        snow_data_struct *, soil_con_struct *);
void calc_net_veg(double, double, double, veg_var_struct *);
void calc_soil_infil(double, double, double, double, double, double, double, double *, double *);
void calc_infil_runoff(double, double, double, double, double *, double, double, double *);
void calc_sat_runoff(double, double, double, double, double *);
int calc_surf_humidity(double, double, double, energy_bal_struct *,
                       snow_data_struct *, cell_data_struct *, soil_con_struct *);
int calc_biomass_heat(double *, double *, double *, double *, double *, double *,
                      double *, veg_var_struct *, veg_lib_struct *);
void calc_dynamicVIC(double, double *, size_t, double *, double *, 
                    cell_data_struct *, soil_con_struct *);
void canopy_two_stream(size_t, size_t, double, double *, energy_bal_struct *, 
                       cell_data_struct *, veg_var_struct *, veg_lib_struct *);
double compute_coszen(double, double, double, unsigned short int, unsigned int);
void compute_soil_resis(cell_data_struct *, soil_con_struct *);
void correct_precip(double *, double, double, double, double);
void ci_func_PHS(bool, double, double, double *, double *, double *, double *,
                 double *, double *, double *, double, double, double, double,
                 double, double, double, double, double, double, double,
                 double, double, double, double, double, double, double,
                 double, double, double, double, double, double,
                 cell_data_struct *, soil_con_struct *, veg_var_struct *, veg_lib_struct *);
int distribute_node_moisture_properties(cell_data_struct *, soil_con_struct *);
void distribute_snow_state(snow_data_struct *);
double devries_weight(double, double, double);
double dlplc(double, double, double);
void enthalpy_to_state(cell_data_struct *, energy_bal_struct *);
int FrictionVelocity(double, double, double *, double *, double *, double *,
                      double *, double);
int frozen_soil(size_t, double, double *, double *, double *, soil_con_struct *);
int func_canopy_energy_bal(double, double, double, double,
                           double, double, double, double,
                           energy_bal_struct *, cell_data_struct *, 
                           snow_data_struct *, soil_con_struct *,
                           veg_var_struct *, veg_lib_struct *);
int func_surf_energy_bal(double, double, double, double, double,
                         double, double, energy_bal_struct *, 
                         cell_data_struct *,
                         snow_data_struct *, soil_con_struct *);
double ft(double, double);
double fth(double, double, double, double);
double fth25(double, double);
int GlacierTemperature(double, cell_data_struct *, energy_bal_struct *, 
                        snow_data_struct *, soil_con_struct *);
void GroundAlbedo(double, double, double *, double *, double *, 
                  double *, double *, double *, double *, double *);
void get_qflx(bool, double, double, double, double, double, double, 
              double *, double *, double *, double *, veg_var_struct *);
void getvegwp(double *, double, double *, double *, double, double, 
              double *, double, double, double, cell_data_struct *, 
              soil_con_struct *, veg_var_struct *, veg_lib_struct *);
void hybrid_PHS(double *, double *, double *, double *, double *,
                double, double, double, double, double, double,
                double, double, double, double, double, double,
                double, double, double, double, double, double,
                double, double, double, double, double *, double *,
                cell_data_struct *, soil_con_struct *, 
                veg_var_struct *, veg_lib_struct *);
void initialize_roughness(bool, double, double, double *,
                          cell_data_struct *, veg_var_struct *);
double initialize_MOST(double, double, double, double, double, double *);
double linear_interp(double, double, double, double, double);
double new_snow_density(double);
int PhotoHydroStress(double, double, double, double, double, double, double,
                     double, cell_data_struct *, soil_con_struct *,
                     veg_var_struct *, veg_lib_struct *);
void prepare_full_energy(double, cell_data_struct *, energy_bal_struct *,
                         snow_data_struct *, soil_con_struct *);
double plc(double, double);
int runoff(double, double, cell_data_struct *, soil_con_struct *);
void set_node_parameters(size_t, double *, double *, double *, double *);
void snow_albedo(double, double, double *, double *);
double snow_density(snow_data_struct *, double, double, double, double);
int snow_hydrology(size_t, double, double, double, double, 
                   force_data_struct *, energy_bal_struct *, cell_data_struct *, 
                   snow_data_struct *, soil_con_struct *);
int snow_intercept(double, double, double *, double *, double,
                   snow_data_struct *, veg_var_struct *);
double snow_aging(double, double, double, snow_data_struct *);
void snow_combination(double, cell_data_struct *, snow_data_struct *);
void snow_compaction(double, double, double, energy_bal_struct *, snow_data_struct *);
void snow_division(snow_data_struct *);
double soil_conductivity(double, double, double, double, double, 
                         double, double, double, double, double);
int SoilTemperature(double, double, cell_data_struct *, energy_bal_struct *,
                    snow_data_struct *, soil_con_struct *);
int soil_transp(cell_data_struct *, soil_con_struct *);
int soil_thermal_fluxes(double, double, double, cell_data_struct *, energy_bal_struct *, 
                             snow_data_struct *);
int soil_vapor_flux(double, cell_data_struct *, soil_con_struct *);
void soil_hydraulic_conductivity(cell_data_struct *, soil_con_struct *);
int surface_albedo(double, double, double, energy_bal_struct *, 
                   cell_data_struct *, snow_data_struct *, soil_con_struct *,
                   veg_var_struct *, veg_lib_struct *);
int surface_fluxes(double, double, double, double, force_data_struct *, 
                   energy_bal_struct *, cell_data_struct *, 
                   snow_data_struct *, soil_con_struct *, 
                   veg_var_struct *, veg_lib_struct *);
int surface_fluxes_glac(double, double, double, force_data_struct *, 
                        energy_bal_struct *, global_param_struct *,
                        cell_data_struct *, snow_data_struct *, soil_con_struct *);                        
void svp_flags(double, double, double *, double *, double *, double *, int);
double SoilWaterRetentionCurve(int, size_t, double, double, soil_con_struct *);
double sign(double a, double b);
double StabilityFunc1(double);
double StabilityFunc2(double);
void spacF(double *, double *, double *, double, double, cell_data_struct *, 
           soil_con_struct *, veg_var_struct *, veg_lib_struct *);
double solve_quadratic_min(double, double, double);
void solve_quadratic(double, double, double, double *, double *);
void update_snow(double, double, double, snow_data_struct *);
void update_node(cell_data_struct *, snow_data_struct *, soil_con_struct *);
void update_fluxes(energy_bal_struct *, cell_data_struct *, snow_data_struct *,
                  soil_con_struct *, veg_lib_struct *);
void update_snow_fluxes(double *, double *, double *, double *,
                           double, double, double, double);
int vic_run(force_data_struct *, all_vars_struct *,
            global_param_struct *, soil_con_struct *,
            veg_con_struct *, veg_lib_struct *);
double volumetric_heat_capacity(double, double, double, double, double,
                                double, double, double, double);
double water_curve_deriv(size_t, double, double, double, soil_con_struct *);
double wrap_compute_zwt(double, cell_data_struct *, soil_con_struct *);

#endif
