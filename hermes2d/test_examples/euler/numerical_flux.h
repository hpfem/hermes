#ifndef NUMERICAL_FLUX_H
#define NUMERICAL_FLUX_H
#include "euler_util.h"

class NumericalFlux
{
public:
  NumericalFlux(double kappa);

  /// Calculates all components of the flux.
  /// Stores the result in the array result.
  virtual void numerical_flux(double result[4], double w_L[4], double w_R[4],
          double nx, double ny) = 0;
  
  /// Calculates a specified component of the flux.
  /// Returns the result.
  virtual double numerical_flux_i(int component, double w_L[4], double w_R[4],
          double nx, double ny) = 0;

  virtual void numerical_flux_solid_wall(double result[4], double w_L[4], double nx, double ny) = 0;
  
  virtual double numerical_flux_solid_wall_i(int component, double w_L[4], double nx, double ny) = 0;

  virtual void numerical_flux_inlet(double result[4], double w_L[4], double w_B[4],
          double nx, double ny) = 0;
  
  virtual double numerical_flux_inlet_i(int component, double w_L[4], double w_B[4],
          double nx, double ny) = 0;

  virtual void numerical_flux_outlet(double result[4], double w_L[4], double pressure, double nx, double ny) = 0;
  
  virtual double numerical_flux_outlet_i(int component, double w_L[4], double pressure, double nx, double ny) = 0;

  /// Rotates the state_vector into the local coordinate system.
  void Q(double result[4], double state_vector[4], double nx, double ny);

  /// Rotates the state_vector back from the local coordinate system.
  void Q_inv(double result[4], double state_vector[4], double nx, double ny);

  void f_1(double result[4], double state[4]);

  // Poisson adiabatic ant = c_p/c_v = 1 + R/c_v.
  double kappa;
};

class StegerWarmingNumericalFlux : public NumericalFlux
{
public:
  StegerWarmingNumericalFlux(double kappa);

  virtual void numerical_flux(double result[4], double w_L[4], double w_R[4],
          double nx, double ny);
  
  virtual double numerical_flux_i(int component, double w_L[4], double w_R[4],
          double nx, double ny);

  void P_plus(double* result, double w[4], double param[4],
          double nx, double ny);

  void P_minus(double* result, double w[4], double param[4],
          double nx, double ny);

  // Also calculates the speed of sound.
  void Lambda_plus(double result[4]);

  // Also calculates the speed of sound.
  void Lambda_minus(double result[4]);

  // Calculates all eigenvalues.
  void Lambda(double result[4]);

  void T_1(double result[4][4]);
  void T_2(double result[4][4]);
  void T_3(double result[4][4]);
  void T_4(double result[4][4]);

  void T_inv_1(double result[4][4]);
  void T_inv_2(double result[4][4]);
  void T_inv_3(double result[4][4]);
  void T_inv_4(double result[4][4]);

  virtual void numerical_flux_solid_wall(double result[4], double w_L[4], double nx, double ny);
  
  virtual double numerical_flux_solid_wall_i(int component, double w_L[4], double nx, double ny);

  virtual void numerical_flux_inlet(double result[4], double w_L[4], double w_B[4],
          double nx, double ny);
  
  virtual double numerical_flux_inlet_i(int component, double w_L[4], double w_B[4],
          double nx, double ny);

  virtual void numerical_flux_outlet(double result[4], double w_L[4], double pressure, double nx, double ny);
  
  virtual double numerical_flux_outlet_i(int component, double w_L[4], double pressure, double nx, double ny);

  double* get_q();

protected:
  // States.
  double q[4];
  double q_1[4];
  double q_L[4];
  double q_R[4];
  double q_L_star[4];
  double q_R_star[4];
  double q_B[4]; // Boundary

  // Speeds of sound.
  double a;
  double a_L;
  double a_R;
  double a_L_star;
  double a_R_star;
  double a_B; // Boundary.

  // x-velocity, y-velocity, magnitude.
  double u, v, V;
};

class VijayasundaramNumericalFlux : public StegerWarmingNumericalFlux
{
public:
  VijayasundaramNumericalFlux(double kappa);

  virtual void numerical_flux(double result[4], double w_L[4], double w_R[4],
          double nx, double ny);
};


class OsherSolomonNumericalFlux : public NumericalFlux
{
public:
  OsherSolomonNumericalFlux(double kappa);

  virtual void numerical_flux(double result[4], double w_L[4], double w_R[4],
          double nx, double ny);
  
  virtual double numerical_flux_i(int component, double w_L[4], double w_R[4],
          double nx, double ny);

  virtual void numerical_flux_solid_wall(double result[4], double w_L[4], double nx, double ny);
  
  virtual double numerical_flux_solid_wall_i(int component, double w_L[4], double nx, double ny);

  virtual void numerical_flux_inlet(double result[4], double w_L[4], double w_R[4],
          double nx, double ny);
  
  virtual double numerical_flux_inlet_i(int component, double w_L[4], double w_R[4],
          double nx, double ny);

  virtual void numerical_flux_outlet(double result[4], double w_L[4], double pressure, double nx, double ny);
  
  virtual double numerical_flux_outlet_i(int component, double w_L[4], double pressure, double nx, double ny);

protected:
  void calculate_q_1_a_1_a_3();
  
  void calculate_q_L_star();

  void calculate_q_3();

  void calculate_q_R_star();


  // States.
  double q_L[4];
  double q_R[4];
  double q_L_star[4];
  double q_R_star[4];
  double q_1[4];
  double q_3[4];
  double q_B[4]; // Boundary

  // Speeds of sound.
  double a_L;
  double a_R;
  double a_L_star;
  double a_R_star;
  double a_1;
  double a_3;
  double a_B; // Boundary.

  // Utility quantities.
  double z_L, z_R, s_L, s_R, alpha;

};

#endif