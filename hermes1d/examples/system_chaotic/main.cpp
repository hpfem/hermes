#include "hermes1d.h"

// ********************************************************************
// This example solves a nonlinear system of four first-order equations
// x1' - DAMPING*(1 - x2^2)*x1   +   x2 = 0
// x2'         -x1   +   x3 = 0
// x3'         -x2   +   x4 = 0
// x4'         -x3          = 0

// in an interval (0, 10) equipped with Dirichlet bdy conditions
// x1(0) = 1, x2(0) = 0, x3(0) = 0, x4(0) = 0

// General input:
static int N_eq = 4;
int N_elem = 500;           // number of elements
double A = 0, B = 10;       // domain end points
int P_init = 2;             // initial polynomal degree

// Damping parameter
int DAMPING_STEPS = 20;  // Number of damping steps. The entire problem
                         // will be run repeatedly, with the DAMPING parameter 
                         // increased from 0 to 1 in DAMPING_STEPS. Every time, 
                         // the last result is used as initial cond. for the 
                         // new computation.   
double DAMPING = 1.0;    // DAMPING is an artificial param. used to 
                         // reduce the strength of the nonlinearity. 
                         // (The nonlinearity is multiplied with it.)

// Newton's method
double NEWTON_TOL = 1e-5;
int NEWTON_MAXITER = 150;

// Boundary conditions
double Val_dir_left_1 = 1;
double Val_dir_left_2 = 0;
double Val_dir_left_3 = 0;
double Val_dir_left_4 = 0;

// Weak forms for Jacobi matrix and residual
#include "forms.cpp"

/******************************************************************************/
int main() {
  // Create coarse mesh, set Dirichlet BC, enumerate 
  // basis functions
  Mesh *mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  mesh->set_bc_left_dirichlet(0, Val_dir_left_1);
  mesh->set_bc_left_dirichlet(1, Val_dir_left_2);
  mesh->set_bc_left_dirichlet(2, Val_dir_left_3);
  mesh->set_bc_left_dirichlet(3, Val_dir_left_4);
  printf("N_dof = %d\n", mesh->assign_dofs());

  // Register weak forms
  DiscreteProblem *dp = new DiscreteProblem();
  dp->add_matrix_form(0, 0, jacobian_1_1);
  dp->add_matrix_form(0, 1, jacobian_1_2);
  dp->add_matrix_form(1, 0, jacobian_2_1);
  dp->add_matrix_form(1, 1, jacobian_2_2);
  dp->add_matrix_form(1, 2, jacobian_2_3);
  dp->add_matrix_form(2, 1, jacobian_3_2);
  dp->add_matrix_form(2, 2, jacobian_3_3);
  dp->add_matrix_form(2, 3, jacobian_3_4);
  dp->add_matrix_form(3, 2, jacobian_4_3);
  dp->add_matrix_form(3, 3, jacobian_4_4);
  dp->add_vector_form(0, residual_1);
  dp->add_vector_form(1, residual_2);
  dp->add_vector_form(2, residual_3);
  dp->add_vector_form(3, residual_4);

  // Newton's loop
  newton(dp, mesh, NULL, NEWTON_TOL, NEWTON_MAXITER);

  // Plot the solution
  Linearizer l(mesh);
  l.plot_solution("solution.gp");

  printf("Done.\n");
  return 1;
}
