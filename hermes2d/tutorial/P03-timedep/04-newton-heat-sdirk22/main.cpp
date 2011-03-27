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

const std::string BDY_DIRICHLET = "1";

// Weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[]) 
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(BDY_DIRICHLET, INIT_BDY_REF_NUM);

  // Initialize boundary conditions.
  EssentialBCNonConst bc_essential(BDY_DIRICHLET);
  EssentialBCs bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();
  info("ndof = %d.", ndof);

  // Previous time level solution (initialized by the initial condition).
  InitialSolutionHeatTransfer u_prev_time(&mesh);

  // Stage 1 Solution (initialized by the initial condition).
  InitialSolutionHeatTransfer sdirk_stage_sol(&mesh);

  // Initialize the weak formulation
  WeakFormHeatTransferNewtonTimedepSDIRKStage1 wf1(ALPHA, time_step, &u_prev_time, BUTCHER_A_11, GAMMA, BUTCHER_C_1);
  WeakFormHeatTransferNewtonTimedepSDIRKStage2 wf2(ALPHA, time_step, &u_prev_time, &sdirk_stage_sol, BUTCHER_A_11, GAMMA, BUTCHER_B_1, BUTCHER_B_2, BUTCHER_C_1, BUTCHER_C_2);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  info("Projecting initial condition to obtain initial vector for the Newton's method.");
  scalar* coeff_vec1 = new scalar[ndof];
  scalar* coeff_vec2 = new scalar[ndof];
  OGProjection::project_global(&space, &u_prev_time, coeff_vec1, matrix_solver);
  for (int i = 0; i < ndof; i++) 
    coeff_vec2[i] = coeff_vec1[i];

  // Initialize views.
  ScalarView sview("Solution", new WinGeom(0, 0, 500, 400));
  OrderView oview("Mesh", new WinGeom(520, 0, 450, 400));
  oview.show(&space);

  // Initialize the weak formulation and the FE problem.
  bool is_linear = false;
  DiscreteProblem dp1(&wf1, &space, is_linear);
  DiscreteProblem dp2(&wf2, &space, is_linear);

  // Time stepping loop:
  double current_time = 0; int ts = 1;
  do {
    // Set current time for the time-dependent rhs.
    wf1.set_current_time(current_time);
    wf2.set_current_time(current_time);

    // Perform Newton's iteration for sdirk_stage_sol.
    info("SDIRK-22 time step (t = %g, tau = %g)", current_time, time_step);
    info("---- Stage I:");
    bool verbose = true;
    if (!hermes2d.solve_newton(coeff_vec1, &dp1, solver, matrix,
		      rhs, NEWTON_TOL, NEWTON_MAX_ITER, verbose))
      error("Newton's iteration did not converge."); 

    // Convert the vector coeff_vec1 into a Solution.
    Solution::vector_to_solution(coeff_vec1, &space, &sdirk_stage_sol);

    // Perform Newton's iteration for the final solution.
    info("---- Stage II:");

    if (!hermes2d.solve_newton(coeff_vec2, &dp2, solver, matrix,
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
  delete [] coeff_vec2;
  delete matrix;
  delete solver;

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
