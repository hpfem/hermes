#ifndef __HERMES2D_EULER_NUMERICAL_FLUX_H
#define __HERMES2D_EULER_NUMERICAL_FLUX_H

#include "hermes2d.h"

// physical constants
#define R 287.14           // Gas constant [J/(kg*K)]
#define c_v 717.5       // specific heat capacity [J/(kg*K)]
#define kappa (1 + R/c_v)

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

#endif
