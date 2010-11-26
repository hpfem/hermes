#ifndef __HERMES2D_EULER_NUMERICAL_FLUX_H
#define __HERMES2D_EULER_NUMERICAL_FLUX_H

#include "hermes2d.h"

class HERMES_API NumericalFlux
{
public:
  NumericalFlux(double kappa) : kappa(kappa) {};

  double f_x(int i, double w0, double w1, double w3, double w4);
  double f_z(int i, double w0, double w1, double w3, double w4);

  double A_x(int i, int j, double w0, double w1, double w3, double w4);
  double A_z(int i, int j, double w0, double w1, double w3, double w4);

  double matrix_R(int i, int j, double w0, double w1, double w3, double w4);
  double matrix_R_inv(int i, int j, double w0, double w1, double w3, double w4);
  double matrix_D_minus(int i, int j, double w0, double w1, double w3, double w4);
  void T_rot(double result[4][4], double beta);
  void dot_vector(double result[4], double A[4][4], double B[4]);
  void riemann_solver(double result[4], double w_l[4], double w_r[4]);
  void riemann_solver_invert(double result[4], double w_l[4], double w_r[4]);
  void numerical_flux(double result[4], double w_l[4], double w_r[4],
          double nx, double ny);
  double numerical_flux_i(int i, double w_l[4], double w_r[4],
          double nx, double ny);

  double kappa;

protected:

  void dot(double result[4][4], double A[4][4], double B[4][4]);
  void A_minus(double result[4][4], double w0, double w1, double w3, double w4);
};

#endif
