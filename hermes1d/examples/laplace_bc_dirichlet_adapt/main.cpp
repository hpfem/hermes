#include "hermes1d.h"

// ********************************************************************

// This example uses hp-adaptivity to solve the Poisson equation 
// -u'' - f = 0 in an interval (A, B), equipped with Dirichlet 
// boundary conditions on both end points. Global reference solution 
// is used both to decide what elements will be refined, and how they 
// will be refined. 

// General input:
static int N_eq = 1;
int N_elem = 2;                         // Number of elements
double A = -M_PI, B = M_PI;             // Domain end points
//double A = -1, B = 1;                 // Domain end points
int P_init = 1;                         // Initial polynomal degree

// Stopping criteria for Newton
double NEWTON_TOL_COARSE = 1e-10;       // Coarse mesh
double NEWTON_TOL_REF = 1e-10;          // Reference mesh
int NEWTON_MAXITER = 150;

// Adaptivity
const int ADAPT_TYPE = 0;               // 0... hp-adaptivity
                                        // 1... h-adaptivity
                                        // 2... p-adaptivity
const double THRESHOLD = 0.7;           // Refined will be all elements whose error
                                        // is greater than THRESHOLD*max_elem_error
const double TOL_ERR_REL = 1e-6;        // Tolerance for the relative error between 
                                        // the coarse mesh and reference solutions
const int NORM = 1;                     // To measure errors:
                                        // 1... H1 norm
                                        // 0... L2 norm

// Boundary conditions
double Val_dir_left = 0;                // Dirichlet condition left
double Val_dir_right = 0;               // Dirichlet condition right

// Function f(x)
double alpha = 1.5;
double f(double x) {
  //if (x <= 0) return 0;
  //else return -alpha*(alpha-1.)*pow(x, alpha-2.); 
  return sin(x);
  //return 2;
}

// Exact solution:
// When changing exact solution, do not 
// forget to update interval accordingly
const int EXACT_SOL_PROVIDED = 1;
double exact_sol(double x, double u[MAX_EQN_NUM], double dudx[MAX_EQN_NUM]) {
  //if (x <= 0) u[0] = 0;
  //else u[0] = pow(x, alpha);
  //if (x <= 0) dudx[0] = 0;
  //else dudx[0] = alpha*pow(x, alpha-1.);
  u[0] = sin(x);
  dudx[0] = cos(x);
  //u[0] = 1. - x*x;
  //dudx[0] = -2.*x;
}

// Weak forms for Jacobi matrix and residual
#include "forms.cpp"

/******************************************************************************/
int main() {
  // Create coarse mesh, set Dirichlet BC, enumerate basis functions
  Mesh *mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  mesh->set_bc_left_dirichlet(0, Val_dir_left);
  mesh->set_bc_right_dirichlet(0, Val_dir_right);
  printf("N_dof = %d\n", mesh->assign_dofs());

  // Create discrete problem on coarse mesh
  DiscreteProblem *dp = new DiscreteProblem();
  dp->add_matrix_form(0, 0, jacobian);
  dp->add_vector_form(0, residual);

  // Initial Newton's loop on coarse mesh
  printf("Newton on coarse mesh.\n");
  newton(dp, mesh, NULL, NEWTON_TOL_COARSE, NEWTON_MAXITER);

  // Replicate coarse mesh including solution.
  Mesh *mesh_ref = mesh->replicate();

  // Refine entire mesh_ref uniformly in 'h' and 'p'
  // Solution is transfered to new elements.
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
    printf("Newton on fine mesh.\n");
    newton(dp, mesh_ref, NULL, NEWTON_TOL_REF, NEWTON_MAXITER);

    // Starting with second adaptivity step, obtain new coarse 
    // mesh solution via Newton's method. Initial condition is 
    // the last coarse mesh solution.
    if (adapt_iterations > 1) {

      // Newton's loop on coarse mesh
      printf("Newton on coarse mesh.\n");
      newton(dp, mesh, NULL, NEWTON_TOL_COARSE, NEWTON_MAXITER);
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
      double err_exact_rel = err_exact_total/exact_sol_norm;
      printf("Relative error (exact) = %g %%\n", 100.*err_exact_rel);
      graph.add_values(0, mesh->get_n_dof(), 100 * err_exact_rel);
    }

    // add entry to DOF convergence graph
    graph.add_values(1, mesh->get_n_dof(), 100 * err_est_rel);

    // Decide whether the relative error is sufficiently small
    if(err_est_rel*100 < TOL_ERR_REL) break;

    // debug
    if (adapt_iterations == 8) break;

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

  printf("Done.\n");
  return 1;
}


















