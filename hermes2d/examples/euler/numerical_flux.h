#ifndef NUMERICAL_FLUX_H
#define NUMERICAL_FLUX_H
#include "euler_util.h"

class NumericalFlux
{
public:
  NumericalFlux();

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

protected:
  /// Rotates the state_vector into the local coordinate system.
  void Q(double result[4], double state_vector[4], double nx, double ny);

  /// Rotates the state_vector back from the local coordinate system.
  void Q_inv(double result[4], double state_vector[4], double nx, double ny);
};

class VijayasundaramNumericalFlux : public NumericalFlux
{
public:
  VijayasundaramNumericalFlux();

  virtual void numerical_flux(double result[4], double w_L[4], double w_R[4],
          double nx, double ny);
  
  virtual double numerical_flux_i(int component, double w_L[4], double w_R[4],
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

  void f_1(double result[4], double state[4]);

  // Poisson adiabatic ant = c_p/c_v = 1 + R/c_v.
  double kappa;

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