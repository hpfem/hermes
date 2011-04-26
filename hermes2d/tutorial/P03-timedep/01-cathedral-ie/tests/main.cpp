#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

// This test makes sure that example 09-timedep-basic works correctly.

const int P_INIT = 2;                             // Polynomial degree of all mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 3;                   // Number of initial uniform mesh refinements towards the boundary.
const double time_step = 3e+2;                    // Time step in seconds.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
const double TEMP_INIT = 10;         // Temperature of the ground (also initial temperature).
const double ALPHA = 10;             // Heat flux coefficient for Newton's boundary condition.
const double LAMBDA = 1e2;           // Thermal conductivity of the material.
const double HEATCAP = 1e2;          // Heat capacity.
const double RHO = 3000;             // Material density.
const double T_FINAL = 5*time_step;

// Weak forms.
#include "../definitions.cpp"

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../cathedral.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary("Boundary air", INIT_REF_NUM_BDY);
  mesh.refine_towards_boundary("Boundary ground", INIT_REF_NUM_BDY);

  // Previous time level solution (initialized by the external temperature).
  Solution tsln(&mesh, TEMP_INIT);

  // Initialize the weak formulation.
  double current_time = 0;
  CustomWeakFormHeatRK1 wf("Boundary air", ALPHA, LAMBDA, HEATCAP, RHO, time_step,
                           &current_time, TEMP_INIT, T_FINAL, &tsln);

  // Initialize boundary conditions.
  DefaultEssentialBCConst bc_essential("Boundary ground", TEMP_INIT);
  EssentialBCs bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();
  info("ndof = %d", ndof);

  // Initialize the FE problem.
  DiscreteProblem dp(&wf, &space);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
  solver->set_factorization_scheme(HERMES_REUSE_FACTORIZATION_COMPLETELY);

  // Initial coefficient vector for the Newton's method.
  scalar* coeff_vec = new scalar[ndof];
  memset(coeff_vec, 0, ndof*sizeof(scalar));

  // Initialize views.
  //ScalarView Tview("Temperature", new WinGeom(0, 0, 450, 600));
  //Tview.set_min_max_range(0,20);
  //Tview.fix_scale_width(30);

  // Time stepping:
  int ts = 1;
  bool jacobian_changed = true;
  do
  {
    info("---- Time step %d, time %3.5f s", ts, current_time);

    // Perform Newton's iteration.
    if (!hermes2d.solve_newton(coeff_vec, &dp, solver, matrix, rhs,
         jacobian_changed)) error("Newton's iteration failed.");
    jacobian_changed = false;

    // Translate the resulting coefficient vector into the Solution sln.
    Solution::vector_to_solution(coeff_vec, &space, &tsln);

    // Visualize the solution.
    //char title[100];
    //sprintf(title, "Time %3.2f s", current_time);
    //Tview.set_title(title);
    //Tview.show(&tsln);

    // Increase current time and time step counter.
    current_time += time_step;
    ts++;
  }
  while (current_time < T_FINAL);

  info("Coordinate (-2.0, 2.0) value = %lf", tsln.get_pt_value(-2.0, 2.0));
  info("Coordinate (-1.0, 2.0) value = %lf", tsln.get_pt_value(-1.0, 2.0));
  info("Coordinate ( 0.0, 2.0) value = %lf", tsln.get_pt_value(0.0, 2.0));
  info("Coordinate ( 1.0, 2.0) value = %lf", tsln.get_pt_value(1.0, 2.0));
  info("Coordinate ( 2.0, 2.0) value = %lf", tsln.get_pt_value(2.0, 2.0));

  bool success = true;

  if (fabs(tsln.get_pt_value(-2.0, 2.0) - 10.001043) > 1E-6) success = false;
  if (fabs(tsln.get_pt_value(-1.0, 2.0) - 10.000171) > 1E-6) success = false;
  if (fabs(tsln.get_pt_value(0.0, 2.0) - 10.000085) > 1E-6) success = false;
  if (fabs(tsln.get_pt_value(1.0, 2.0) - 10.000171) > 1E-6) success = false;
  if (fabs(tsln.get_pt_value(2.0, 2.0) - 10.001043) > 1E-6) success = false;

  if (success) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}

