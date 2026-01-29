/******************************************************************************
 * @section DESCRIPTION
 *
 * Header file for vic_run routines
 *****************************************************************************/

#ifndef VIC_RUN_H
#define VIC_RUN_H

#include <vic_def.h>

void AdvectedEnergy(double, double, double, double, double,
                    double, double, double, double, double,
                    double *, double *, double *);
void AdvectedEnergyGlac(double, double, double, double, double *);
double calc_atmos_energy_bal(double, double, double, double, double, double,
                             double, double, double, double, double *,
                             energy_bal_struct *, cell_data_struct *,
                             veg_var_struct *);
bool assert_close_double(double x, double y, double rtol, double abs_tol);
bool assert_close_float(float x, float y, float rtol, float abs_tol);

int calc_energy_bal(size_t, double, double, double, double, double,
                    double *, double *, force_data_struct *,
                    energy_bal_struct *, cell_data_struct *,
                    snow_data_struct *, soil_con_struct *, 
                    veg_var_struct *, veg_lib_struct *);
int calc_energy_bal_glac(size_t, double, double, double, double, double,
                            double *, double *, force_data_struct *, 
                            energy_bal_struct *, cell_data_struct *,
                            snow_data_struct *, soil_con_struct *);
int calc_water_bal(size_t, double, double, double, double, force_data_struct *,
                   energy_bal_struct *, cell_data_struct *,
                   snow_data_struct *, veg_var_struct *,
                   soil_con_struct *);
int calc_water_bal_glac(size_t, double, double, double, double, double, double,
                        force_data_struct *, cell_data_struct *, 
                        snow_data_struct *, soil_con_struct *);
double calc_longwave(double, double);
double calc_sensible_heat(double, double, double);
double calc_latent_heat(double, double, double, double);
double calc_ground_heat(double, double, double);
void calc_root_fractions(veg_con_struct *, soil_con_struct *);
void calc_rainonly(double, double, double, double, double, double, double *, double *);
int calc_rc(size_t, double, double, double, 
               energy_bal_struct *, veg_var_struct *, veg_lib_struct *);
int CalcPhaseChange(double, energy_bal_struct *, cell_data_struct *,
                    snow_data_struct *, soil_con_struct *);
int calc_soil_thermal_fluxes(int, double *, double *, char *, unsigned int *,
                             double *, double *, double *, double *, double *,
                             double *, double *, double *, double *, double *,
                             double *, int, int, int);
double calc_veg_displacement(double, double);
void calc_snow_coverage(double, unsigned short int, 
                        snow_data_struct *, soil_con_struct *);
void calc_net_veg(double, double, double, veg_var_struct *);
int CalcAerodynamic(size_t, double, double, double *, double *, double *, double *, 
                    double, double, double, double, double);
void calc_soil_infil(double, double, double, double, double, double, double, double *, double *);
void calc_infil_runoff(double, double, double, double, double *, double, double, double *);
void calc_sat_runoff(double, double, double, double, double *);
void calc_soilmoist_Richards(double *, double *, double *, double *,
                             double, double *, cell_data_struct *, soil_con_struct *);
void calc_dynamic_runoff(double, double *, size_t, double *, double *, 
                        cell_data_struct *, soil_con_struct *);
double canopy_evap(double, energy_bal_struct *, cell_data_struct *, 
                   snow_data_struct *, veg_var_struct *);
void canopy_two_stream(size_t, size_t, double, double, double, double, double,
                       double *, double, energy_bal_struct *, veg_lib_struct *);
void colavg(double *, double *, double *, double, double *, int, double,
            double);
double compute_coszen(double, double, double, unsigned short int, unsigned int);
double compute_zwt(double, double, cell_data_struct *, soil_con_struct *);
void correct_precip(double *, double, double, double, double);
int distribute_node_moisture_properties(bool, double *, double *, double *, double *,
                                        double *, double *, double *, char);
void distribute_snow_state(bool, double, snow_data_struct *);
void eddy(int, double, double *, double *, double, int, double, double);
void fdjac3(double *, double *, double *, double *, double *, void (*vecfunc)(
                double *, double *, int, int, ...), int);
int FrictionVelocity(double, double, double, double *, double *,
                      double *, double, double);
void find_0_degree_fronts(energy_bal_struct *, double *, double *, int);
int func_canopy_energy_bal(double, double, double, double, double, double,
                           double, double, double, double, double *, double *,
                           energy_bal_struct *, cell_data_struct *, 
                           snow_data_struct *, soil_con_struct *,
                           veg_var_struct *, veg_lib_struct *);
int func_surf_energy_bal(double, double, double, double, double, double,
                         double, double, double, double, double, 
                         double *, double *, energy_bal_struct *, 
                         cell_data_struct *, snow_data_struct *,
                         soil_con_struct *);
int GlacierTemperature(double, cell_data_struct *, energy_bal_struct *,
                       snow_data_struct *, soil_con_struct *);
void GroundAlbedo(double, double, double *, double *, double *, 
                  double *, double *, double *, double *, double *);
double hiTinhib(double);
double initialize_MOST(double, double, double, double, double, double *);
void latsens(double, double, double, double, double, double, double, double,
             double *, double *, double);
double linear_interp(double, double, double, double, double);
double lkdrag(double, double, double, double, double);
void MassRelease(double *, double *, double *, double *);
double maximum_unfrozen_water(double, double, double, double);
double new_snow_density(double);
int newt_raph(void (*vecfunc)(double *, double *, int, int,
                              ...), double *, int);
int PhaseChangeGlac(double, energy_bal_struct *, cell_data_struct *,
                    snow_data_struct *, soil_con_struct *);
void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
void prepare_full_energy(unsigned short int, double, cell_data_struct *, energy_bal_struct *,
                         snow_data_struct *, soil_con_struct *);
void rhoinit(double *, double);
double rtnewt(double x1, double x2, double xacc, double Ur, double Zr);
int runoff(double, cell_data_struct *, soil_con_struct *);
void set_node_parameters(soil_con_struct *);
void snow_albedo(double, double, double *, double *);
double snow_density(snow_data_struct *, double, double, double, double);
int snow_intercept(double, double, double *, double *, double,
                   snow_data_struct *, veg_var_struct *);
double snow_aging(double, double, double *, double, double, double *);
void snow_combination(double, cell_data_struct *, snow_data_struct *);
void snow_compaction(double, double, double, cell_data_struct *, snow_data_struct *);
void snow_division(snow_data_struct *);
double soil_conductivity(double, double, double, double, double);
int SoilTemperature(double, cell_data_struct *, energy_bal_struct *,
                    snow_data_struct *, soil_con_struct *);
void solve_soilmoist_Richards(double *, double *, double *, double *, size_t, double *,
                              double *, double *, double *, double *, double *, double *);
void solve_temperature(double *, double *, double *, double *, double, size_t, 
                       cell_data_struct *, energy_bal_struct *, snow_data_struct *);
int SuperCooling(size_t, double, double, double, double *, soil_con_struct *);
int surface_fluxes(double, double, double, double, double *, double *,
                   force_data_struct *, energy_bal_struct *, global_param_struct *,
                   cell_data_struct *, snow_data_struct *, soil_con_struct *, 
                   veg_var_struct *, veg_lib_struct *);
int surface_fluxes_glac(double, double, double, double, double *, double *,
                        force_data_struct *, energy_bal_struct *, global_param_struct *,
                        cell_data_struct *, snow_data_struct *, soil_con_struct *);                        
void svp(double, double *, double *);
void svp_slope(double, double *, double *);
void s_humid(double, double, double *, double *);
double sign(double a, double b);
double StabilityFunc1(double);
double StabilityFunc2(double);
void tracer_mixer(double *, int *, double *, int, double, double, double *);
void tridia(int, double *, double *, double *, double *, double *);
void tridiag(double *, double *, double *, double *, unsigned int);
void update_snow(double, double, double, snow_data_struct *);
void update_surface_fluxes(double *, double *, double *, double *,
                           double, double, double, double);
int vic_run(force_data_struct *, all_vars_struct *,
            global_param_struct *, soil_con_struct *,
            veg_con_struct *, veg_lib_struct *);
double volumetric_heat_capacity(double, double, double, double, double);
void wrap_compute_zwt(soil_con_struct *, cell_data_struct *);

#endif
