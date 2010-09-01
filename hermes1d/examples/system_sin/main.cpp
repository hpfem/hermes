#include "hermes1d.h"

// ********************************************************************
// This example solves a system of two linear second-order equations 
// u' + k^2 v = 0
// u - v' = 0
// which is equivalent to u'' + k^2 u = 0
// in an interval (0, 2*pi) equipped with Dirichlet bdy conditions 
// u(0) = 0, v(0) = k
// The exact solution is u(x) = sin(k*x), v(x) = k*cos(k*x)

// General input:
static int N_eq = 2;
int N_elem = 20;          // number of elements
double A = 0, B = 2*M_PI;     // domain end points
int P_init = 2;          // initial polynomal degree
double k = 1.0;          // the constant in the equation

// Newton's method
double NEWTON_TOL = 1e-5;
int NEWTON_MAXITER = 150;

// Boundary conditions
double Val_dir_left_0 = 0;
double Val_dir_left_1 = k;

// Weak forms for Jacobi matrix and residual
#include "forms.cpp"

/******************************************************************************/
int main() {
  // Create coarse mesh, set Dirichlet BC, enumerate 
  // basis functions
  Mesh *mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  mesh->set_bc_left_dirichlet(0, Val_dir_left_0);
  mesh->set_bc_left_dirichlet(1, Val_dir_left_1);
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
