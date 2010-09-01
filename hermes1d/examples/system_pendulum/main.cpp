#include "hermes1d.h"

// ********************************************************************
// This example solves the mathematical pendulum equation 
// y'' + k**2 * sin(y) = 0 in an interval (A, B), equipped with the 
// initial conditions y(A) = Init_angle, y'(0) = Init_vel. The 
// system is decomposed into two first order ODE and solved via 
// the Newton's method starting from zero initial condition.
// Note that the method diverges for longer time intervals, 
// depending on the interval length, number of elements, and 
// the initial polynomial degree.
//
// Derivation:
// m*l*u'' = -m*g*sin(u)
// so:
// u'' + k^2 * sin(u) = 0
// with k^2 = g/l
// so we have to solve a system of two nonlinear second-order equations
// v' + k^2 sin u = 0
// u' - v = 0
// in an interval (0, 2*pi) equipped with Dirichlet bdy conditions
// u(0) = 0, v(0) = k
// The approximate (linearized) solution is u(x) = sin(k*x), v(x) = k*cos(k*x)

// General input:
static int N_eq = 2;
int N_elem = 1292;          // number of elements
double A = 0, B = 10;     // domain end points
int P_init = 1;            // initial polynomal degree
double k = 0.5;

// Newton's method
double NEWTON_TOL = 1e-5;
int NEWTON_MAXITER = 150;

// Boundary conditions
double Init_angle = M_PI/2.;      // initial angle
double Init_vel = 0;              // initial velocity

// Weak forms for Jacobi matrix and residual
#include "forms.cpp"

/******************************************************************************/
int main() {
  // Create coarse mesh, set Dirichlet BC, enumerate 
  // basis functions
  Mesh *mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  mesh->set_bc_left_dirichlet(0, Init_angle);
  mesh->set_bc_left_dirichlet(1, Init_vel);
  printf("N_dof = %d\n", mesh->assign_dofs());

  // Register weak forms
  DiscreteProblem *dp = new DiscreteProblem();
  dp->add_matrix_form(0, 0, jacobian_0_0);
  dp->add_matrix_form(0, 1, jacobian_0_1);
  dp->add_matrix_form(1, 0, jacobian_1_0);
  dp->add_matrix_form(1, 1, jacobian_1_1);
  dp->add_vector_form(0, residual_0);
  dp->add_vector_form(1, residual_1);

  // Newton's loop
  newton(dp, mesh, NULL, NEWTON_TOL, NEWTON_MAXITER);

  // Plot the solution
  Linearizer l(mesh);
  l.plot_solution("solution.gp");

  printf("Done.\n");
  return 1;
}
