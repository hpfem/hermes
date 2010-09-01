#include "hermes1d.h"

// ********************************************************************
// This example solves a system of two linear second-order equations 
// - u'' + v - f_0 = 0
// - v'' + u - f_1 = 0
// in an interval (A, B) equipped with Dirichlet bdy conditions 
// u(A) = exp(A), u(B) = exp(B), v(A) = exp(-A), v(B) = exp(-B). 
// The exact solution is u(x) = exp(x), v(x) = exp(-x). 

// General input:
static int N_eq = 2;
int N_elem = 2;          // number of elements
double A = 0, B = 1;     // domain end points
int P_init = 2;          // initial polynomal degree

// Newton's method
double NEWTON_TOL = 1e-5;
int NEWTON_MAXITER = 150;

// Boundary conditions
double Val_dir_left_0 = exp(A);
double Val_dir_right_0 = exp(B);
double Val_dir_left_1 = exp(-A);
double Val_dir_right_1 = exp(-B);

// Function f_0(x)
double f_0(double x) {
  return -exp(x) + exp(-x);
}

// Function f_1(x)
double f_1(double x) {
  return -exp(-x) + exp(x);
}

// Weak forms for Jacobi matrix and residual
#include "forms.cpp"

/******************************************************************************/
int main() {
  // Create coarse mesh, set Dirichlet BC, enumerate 
  // basis functions
  Mesh *mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  mesh->set_bc_left_dirichlet(0, Val_dir_left_0);
  mesh->set_bc_right_dirichlet(0, Val_dir_right_0);
  mesh->set_bc_left_dirichlet(1, Val_dir_left_1);
  mesh->set_bc_right_dirichlet(1, Val_dir_right_1);
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

