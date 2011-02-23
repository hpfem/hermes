#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "function/function.h"

using namespace RefinementSelectors;

//  This example is the same as 03-newton-heat-ie except that time
//  discretization is done via the SDIRK-22 method.
//
//  The Butcher's table of the SDIRK-22 method:

double GAMMA = 1. - 1./sqrt(2.);
double BUTCHER_A_11 = 1. - 1./sqrt(2.);
double BUTCHER_A_12 = 0.;
double BUTCHER_A_21 = 1. - (1. - 1./sqrt(2.));
double BUTCHER_A_22 = 1. - 1./sqrt(2.);
double BUTCHER_B_1 = 1. - (1. - 1./sqrt(2.));
double BUTCHER_B_2 = 1. - 1./sqrt(2.);
double BUTCHER_C_1 = 1. - 1./sqrt(2.);
double BUTCHER_C_2 = 1.;

//  The method can be found in Butcher's book on page 244. 
//
//  PDE: time-dependent heat transfer equation with nonlinear thermal
//  conductivity, du/dt - div[lambda(u)grad u] = f.
//
//  Domain: square (-10, 10)^2.
//
//  BC: Dirichlet, given by the function dir_lift() below.
//  IC: Same function dir_lift().
//
//  The following parameters can be changed:

const int INIT_GLOB_REF_NUM = 3;                  // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 4;                   // Number of initial refinements towards boundary.
const int P_INIT = 2;                             // Initial polynomial degree.
double time_step = 0.2;                           // Time step.
const double T_FINAL = 5.0;                       // Time interval length.
const double NEWTON_TOL = 1e-5;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

const double ALPHA = 4.0;                         // For the nonlinear thermal ocnductivity.
double current_time = 0.0;

// Model parameters.
#include "model.cpp"

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[]) 
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(BDY_DIRICHLET, INIT_BDY_REF_NUM);

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_DIRICHLET);

  // Enter Dirichlet boundary values.
  BCValues bc_values;
  bc_values.add_function(BDY_DIRICHLET, essential_bc_values);   

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d.", ndof);

  // Previous time level solution (initialized by the initial condition).
  Solution u_prev_time(&mesh, init_cond);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  info("Projecting initial condition to obtain initial vector for the Newton's method.");
  scalar* coeff_vec1 = new scalar[ndof];
  OGProjection::project_global(&space, &u_prev_time, coeff_vec1, matrix_solver);

  // Initialize views.
  ScalarView sview("Solution", new WinGeom(0, 0, 500, 400));
  OrderView oview("Mesh", new WinGeom(520, 0, 450, 400));
  oview.show(&space);

  // Initialize the weak formulation and the FE problem.
  WeakForm wf, wf1, wf2;
  Solution sdirk_stage_sol;
  scalar* coeff_vec2;
  DiscreteProblem *dp, *dp1, *dp2;

  sdirk_stage_sol.set_exact(&mesh, init_cond);

  coeff_vec2 = new scalar[ndof];
  for (int i = 0; i < ndof; i++) coeff_vec2[i] = coeff_vec1[i];

  wf1.add_matrix_form(callback(jac_sdirk));
  wf1.add_vector_form(callback(res_sdirk_stage_1), HERMES_ANY, 
                      Hermes::vector<MeshFunction*>(&u_prev_time));
  wf2.add_matrix_form(callback(jac_sdirk));
  wf2.add_vector_form(callback(res_sdirk_stage_2), HERMES_ANY, 
                      Hermes::vector<MeshFunction*>(&u_prev_time, &sdirk_stage_sol));

  bool is_linear = false;
  dp1 = new DiscreteProblem(&wf1, &space, is_linear);
  dp2 = new DiscreteProblem(&wf2, &space, is_linear);

  // Time stepping loop:
  int ts = 1;
  do {
    // Perform Newton's iteration for sdirk_stage_sol.
    info("SDIRK-22 time step (t = %g, tau = %g)", current_time, time_step);
    info("---- Stage I:");
    bool verbose = true;
    if (!solve_newton(coeff_vec1, dp1, solver, matrix,
		      rhs, NEWTON_TOL, NEWTON_MAX_ITER, verbose))
      error("Newton's iteration did not converge."); 

    // Convert the vector coeff_vec1 into a Solution.
    Solution::vector_to_solution(coeff_vec1, &space, &sdirk_stage_sol);

    // Perform Newton's iteration for the final solution.
    info("---- Stage II:");

    if (!solve_newton(coeff_vec2, dp2, solver, matrix,
	       	      rhs, NEWTON_TOL, NEWTON_MAX_ITER, verbose))
      error("Newton's iteration did not converge."); 

    // Translate Y2 into a Solution.
    Solution::vector_to_solution(coeff_vec2, &space, &u_prev_time);
  
    // Update time.
    current_time += time_step;

    // Show the new time level solution.
    char title[100];
    sprintf(title, "Solution, t = %g", current_time);
    sview.set_title(title);
    sview.show(&u_prev_time);

    ts++;
  } while (current_time < T_FINAL);

  // Clean up.
  delete [] coeff_vec1;
  delete dp1;
  delete dp2;
  delete [] coeff_vec2;
  delete matrix;
  delete solver;

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
