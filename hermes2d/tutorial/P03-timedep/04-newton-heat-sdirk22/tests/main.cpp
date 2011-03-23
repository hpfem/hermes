#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "function/function.h"

using namespace RefinementSelectors;

// This test makes sure that the example "19-rk-comparison" works correctly.
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

// The following parameters can be changed:
const int TIME_INTEGRATION = 2;                   // 1 = implicit Euler, 2 = SDIRK-2
const int INIT_GLOB_REF_NUM = 3;                  // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 0;                   // Number of initial refinements towards boundary.
const int P_INIT = 2;                             // Initial polynomial degree.
double time_step = 0.2;                                 // Time step.
const double T_FINAL = 1.0;                       // Time interval length.
const double NEWTON_TOL = 1e-5;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

const double ALPHA = 4.0;                         // For the nonlinear thermal ocnductivity.

const std::string BDY_DIRICHLET = "1";

// Weak forms.
#include "../forms.cpp"

// Initial condition.
#include "initial_condition.cpp"

int main(int argc, char* argv[]) 
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../square.mesh", &mesh);

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
    bool verbose = true;
    if (!hermes2d.solve_newton(coeff_vec1, &dp1, solver, matrix,
			rhs, NEWTON_TOL, NEWTON_MAX_ITER, verbose))
        error("Newton's iteration did not converge."); 

      // Convert the vector coeff_vec1 into a Solution.
      Solution::vector_to_solution(coeff_vec1, &space, &sdirk_stage_sol);

    if (!hermes2d.solve_newton(coeff_vec2, &dp2, solver, matrix,
			rhs, NEWTON_TOL, NEWTON_MAX_ITER, verbose))
        error("Newton's iteration did not converge."); 

      // Translate Y2 into a Solution.
      Solution::vector_to_solution(coeff_vec2, &space, &u_prev_time);
  
    // Update time.
    current_time += time_step;

    ts++;
  } while (current_time < T_FINAL);

  info("Coordinate (-8.0, -8.0) value = %lf", u_prev_time.get_pt_value(-8.0, -8.0));
  info("Coordinate (-5.0, -5.0) value = %lf", u_prev_time.get_pt_value(-5.0, -5.0));
  info("Coordinate (-3.0, -3.0) value = %lf", u_prev_time.get_pt_value(-3.0, -3.0));
  info("Coordinate ( 0.0,  0.0) value = %lf", u_prev_time.get_pt_value(0.0,  0.0));
  info("Coordinate ( 3.0,  3.0) value = %lf", u_prev_time.get_pt_value(3.0,  3.0));
  info("Coordinate ( 5.0,  5.0) value = %lf", u_prev_time.get_pt_value(5.0,  5.0));
  info("Coordinate ( 8.0,  8.0) value = %lf", u_prev_time.get_pt_value(8.0,  8.0));

  double coor_x_y[7] = {-8.0, -5.0, -3.0, 0.0, 3.0, 5.0, 8.0};
  bool success = true;

  if (fabs(u_prev_time.get_pt_value(coor_x_y[0], coor_x_y[0]) - 1.009815) > 1E-6) success = false;
  if (fabs(u_prev_time.get_pt_value(coor_x_y[1], coor_x_y[1]) - 1.267173) > 1E-6) success = false;
  if (fabs(u_prev_time.get_pt_value(coor_x_y[2], coor_x_y[2]) - 1.680294) > 1E-6) success = false;
  if (fabs(u_prev_time.get_pt_value(coor_x_y[3], coor_x_y[3]) - 2.367316) > 1E-6) success = false;
  if (fabs(u_prev_time.get_pt_value(coor_x_y[4], coor_x_y[4]) - 2.749174) > 1E-6) success = false;
  if (fabs(u_prev_time.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.966687) > 1E-6) success = false;
  if (fabs(u_prev_time.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.442497) > 1E-6) success = false;

  // Clean up.
  delete [] coeff_vec1;
    delete [] coeff_vec2;
  delete matrix;
  delete solver;

  if (success) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}
