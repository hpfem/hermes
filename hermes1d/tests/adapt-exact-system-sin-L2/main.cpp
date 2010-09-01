#include "hermes1d.h"

#include "legendre.h"
#include "lobatto.h"
#include "quad_std.h"

// This test makes sure that an exact function 
// sin(K*x), K*cos(K*x) is approximated adaptively 
// with relative error 1e-1% with less than 70 DOF.
// Adaptivity is done in L2 norm.

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

// General input:
static int N_eq = 2;
int N_elem = 3;                         // Number of elements
double A = 0, B = 2*M_PI;               // Domain end points
int P_init = 1;                         // Initial polynomial degree
double K = 1.0;                         // Equation parameter

// Newton's method
const double NEWTON_TOL_COARSE = 1e-6;        // Coarse mesh
const double NEWTON_TOL_REF = 1e-6;           // Reference mesh
const int NEWTON_MAXITER = 150;

// Adaptivity
const int ADAPT_TYPE = 0;               // 0... hp-adaptivity
                                        // 1... h-adaptivity
                                        // 2... p-adaptivity
const double THRESHOLD = 0.7;           // Refined will be all elements whose error
                                        // is greater than THRESHOLD*max_elem_error
const double TOL_ERR_REL = 1e-1;        // Tolerance for the relative error between 
                                        // the coarse mesh and reference solutions
const int NORM = 0;                     // 1... H1 norm
                                        // 0... L2 norm

// Boundary conditions
double Val_dir_left_0 = 0;
double Val_dir_left_1 = K;

// Exact solution:
// When changing exact solution, do not 
// forget to update interval accordingly
const int EXACT_SOL_PROVIDED = 1;
double exact_sol(double x, double u[MAX_EQN_NUM], double dudx[MAX_EQN_NUM]) {
  u[0] = sin(K*x);
  dudx[0] = K*cos(K*x);
  u[1] = K*cos(K*x);
  dudx[1] = -K*K*sin(K*x);
}

// ********************************************************************

// Jacobi matrix block 0, 0 (equation 0, solution component 0)
// Note: u_prev[c][i] contains the values of solution component c 
// in integration point x[i]. similarly for du_prevdx.
double jacobian_0_0(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += dudx[i]*v[i]*weights[i];
  }
  return val;
};

// Jacobi matrix block 0, 1 (equation 0, solution component 1)
double jacobian_0_1(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += -u[i]*v[i]*weights[i];
  }
  return val;
};

// Jacobi matrix block 1, 0 (equation 1, solution component 0)
double jacobian_1_0(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += K*K*u[i]*v[i]*weights[i];
  }
  return val;
};

// Jacobi matrix block 1, 1 (equation 1, solution component 1)
double jacobian_1_1(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += dudx[i]*v[i]*weights[i];
  }
  return val;
};

// Residual part 0 (equation 0) 
double residual_0(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double *v, double *dvdx, void *user_data)
{
  double val = 0;
  int si = 0;      // solution index (only 0 is relevant for this example)
  for(int i = 0; i<num; i++) {
    val += (du_prevdx[si][0][i] - u_prev[si][1][i])*v[i]*weights[i];
  }
  return val;
};

// Residual part 1 (equation 1)
double residual_1(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double *v, double *dvdx, void *user_data)
{
  double val = 0;
  int si = 0;      // solution index (only 0 is relevant for this example)
  for(int i = 0; i<num; i++) {
    val += (K*K*u_prev[si][0][i] + du_prevdx[si][1][i])*v[i]*weights[i];
  }
  return val;
};

// ********************************************************************
int main() {
  // Create coarse mesh, set Dirichlet BC, enumerate basis functions
  Mesh *mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  mesh->set_bc_left_dirichlet(0, Val_dir_left_0);
  mesh->set_bc_left_dirichlet(1, Val_dir_left_1);
  printf("N_dof = %d\n", mesh->assign_dofs());

  // Create discrete problem on coarse mesh
  DiscreteProblem *dp = new DiscreteProblem();
  dp->add_matrix_form(0, 0, jacobian_0_0);
  dp->add_matrix_form(0, 1, jacobian_0_1);
  dp->add_matrix_form(1, 0, jacobian_1_0);
  dp->add_matrix_form(1, 1, jacobian_1_1);
  dp->add_vector_form(0, residual_0);
  dp->add_vector_form(1, residual_1);

  // scipy umfpack solver
  CommonSolverSciPyUmfpack solver;

  // Initial Newton's loop on coarse mesh
  newton(dp, mesh, &solver, NEWTON_TOL_COARSE, NEWTON_MAXITER);

  // Replicate coarse mesh including dof arrays
  Mesh *mesh_ref = mesh->replicate();

  // Refine entire mesh_ref uniformly in 'h' and 'p'
  int start_elem_id = 0; 
  int num_to_ref = mesh_ref->get_n_active_elem();
  mesh_ref->reference_refinement(start_elem_id, num_to_ref);
  printf("Fine mesh created (%d DOF).\n", mesh_ref->get_n_dof());

  // Convergence graph wrt. the number of degrees of freedom
  GnuplotGraph graph;
  graph.set_log_y();
  graph.set_captions("Convergence History", "Degrees of Freedom", "Error [%]");
  graph.add_row("exact error", "k", "-", "o");
  graph.add_row("error estimate", "k", "--");

  // Main adaptivity loop
  int adapt_iterations = 1;
  while(1) {
    printf("============ Adaptivity step %d ============\n", adapt_iterations); 

    // Newton's loop on fine mesh
    newton(dp, mesh_ref, &solver, NEWTON_TOL_REF, NEWTON_MAXITER);

    // Starting with second adaptivity step, obtain new coarse 
    // mesh solution via Newton's method. Initial condition is 
    // the last coarse mesh solution.
    if (adapt_iterations > 1) {

      // Newton's loop on coarse mesh
      newton(dp, mesh, &solver, NEWTON_TOL_COARSE, NEWTON_MAXITER);
    }

    // In the next step, estimate element errors based on 
    // the difference between the fine mesh and coarse mesh solutions. 
    double err_est_array[MAX_ELEM_NUM]; 
    double err_est_total = calc_error_estimate(NORM, mesh, mesh_ref, 
                           err_est_array);

    // Calculate the norm of the fine mesh solution
    double ref_sol_norm = calc_solution_norm(NORM, mesh_ref);

    // Calculate an estimate of the global relative error
    double err_est_rel = err_est_total/ref_sol_norm;
    printf("Relative error (est) = %g %%\n", 100.*err_est_rel);

    // If exact solution available, also calculate exact error
    double err_exact_rel;  
    if (EXACT_SOL_PROVIDED) {
      // Calculate element errors wrt. exact solution
      double err_exact_total = calc_error_exact(NORM, mesh, exact_sol);
     
      // Calculate the norm of the exact solution
      // (using a fine subdivision and high-order quadrature)
      int subdivision = 500; // heuristic parameter
      int order = 20;        // heuristic parameter
      double exact_sol_norm = calc_solution_norm(NORM, exact_sol, N_eq, A, B,
                                                  subdivision, order);
      // Calculate an estimate of the global relative error
      err_exact_rel = err_exact_total/exact_sol_norm;
      printf("Relative error (exact) = %g %%\n", 100.*err_exact_rel);
      graph.add_values(0, mesh->get_n_dof(), 100 * err_exact_rel);
    }

    // add entry to DOF convergence graph
    graph.add_values(1, mesh->get_n_dof(), 100 * err_est_rel);

    // Decide whether the relative error is sufficiently small
    if(err_est_rel*100 < TOL_ERR_REL) break;

    // debug
    // (adapt_iterations == 2) break;
 
    // Returns updated coarse and fine meshes, with the last 
    // coarse and fine mesh solutions on them, respectively. 
    // The coefficient vectors and numbers of degrees of freedom 
    // on both meshes are also updated. 
    adapt(NORM, ADAPT_TYPE, THRESHOLD, err_est_array,
          mesh, mesh_ref);

    adapt_iterations++;
  }

  // Plot meshes, results, and errors
  adapt_plotting(mesh, mesh_ref, 
                 NORM, EXACT_SOL_PROVIDED, exact_sol);

  // Save convergence graph
  graph.save("conv_dof.gp");

  int success_test = 1; 
  printf("N_dof = %d\n", mesh->get_n_dof());
  if (mesh->get_n_dof() > 70) success_test = 0;

  if (success_test) {
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}



