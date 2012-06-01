#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

const bool STOKES = false;                        // For application of Stokes flow (creeping flow).

// If this is defined, the pressure is approximated using
// discontinuous L2 elements (making the velocity discreetely
// divergence-free, more accurate than using a continuous
// pressure approximation). Otherwise the standard continuous
// elements are used. The results are striking - check the
// tutorial for comparisons.
#define PRESSURE_IN_L2

const int P_INIT_VEL = 2;                         // Initial polynomial degree for velocity components.

// Initial polynomial degree for pressure.
// Note: P_INIT_VEL should always be greater than
// P_INIT_PRESSURE because of the inf-sup condition.
const int P_INIT_PRESSURE = 1;

const double RE = 200.0;                          // Reynolds number.
const double VEL_INLET = 1.0;                     // Inlet velocity (reached after STARTUP_TIME).

// During this time, inlet velocity increases gradually
// from 0 to VEL_INLET, then it stays constant.
const double STARTUP_TIME = 1.0;

const double TAU = 0.1;                           // Time step.
const double T_FINAL = 0.21;                      // Time interval length.
const double NEWTON_TOL = 1e-3;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 10;                   // Maximum allowed number of Newton iterations.

// Domain height (necessary to define the parabolic
// velocity profile at inlet).
const double H = 5;

// Possibilities: Hermes::SOLVER_AMESOS, Hermes::SOLVER_AZTECOO, Hermes::SOLVER_MUMPS,
// Hermes::SOLVER_PETSC, Hermes::SOLVER_SUPERLU, Hermes::SOLVER_UMFPACK.
Hermes::MatrixSolverType matrix_solver_type = Hermes::SOLVER_UMFPACK;

// Boundary markers.
const std::string BDY_BOTTOM = "1";
const std::string BDY_RIGHT = "2";
const std::string BDY_TOP = "3";
const std::string BDY_LEFT = "4";
const std::string BDY_OBSTACLE = "5";

// Current time (used in weak forms).
double current_time = 0;

// Weak forms.
#include "../definitions.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("../domain.mesh", &mesh);

  // Initial mesh refinements.
  //mesh.refine_all_elements();
  mesh.refine_towards_boundary(BDY_OBSTACLE, 4, false);
  mesh.refine_towards_boundary(BDY_TOP, 4, true);     // '4' is the number of levels,
  mesh.refine_towards_boundary(BDY_BOTTOM, 4, true);  // 'true' stands for anisotropic refinements.

  // Initialize boundary conditions.
  EssentialBCNonConst bc_left_vel_x(BDY_LEFT, VEL_INLET, H, STARTUP_TIME);
  Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_other_vel_x(Hermes::vector<std::string>(BDY_BOTTOM, BDY_TOP, BDY_OBSTACLE), 0.0);
  Hermes::Hermes2D::EssentialBCs<double> bcs_vel_x(Hermes::vector<EssentialBoundaryCondition<double> *>(&bc_left_vel_x, &bc_other_vel_x));
  Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_vel_y(Hermes::vector<std::string>(BDY_LEFT, BDY_BOTTOM, BDY_TOP, BDY_OBSTACLE), 0.0);
  Hermes::Hermes2D::EssentialBCs<double> bcs_vel_y(&bc_vel_y);
  Hermes::Hermes2D::EssentialBCs<double> bcs_pressure;

  // Spaces for velocity components and pressure.
  H1Space<double> xvel_space(&mesh, &bcs_vel_x, P_INIT_VEL);
  H1Space<double> yvel_space(&mesh, &bcs_vel_y, P_INIT_VEL);
#ifdef PRESSURE_IN_L2
  L2Space<double> p_space(&mesh, P_INIT_PRESSURE);
#else
  H1Space<double> p_space(&mesh, &bcs_pressure, P_INIT_PRESSURE);
#endif

  // Calculate and report the number of degrees of freedom.
  int ndof = Space<double>::get_num_dofs(Hermes::vector<const Space<double> *>(&xvel_space, &yvel_space, &p_space));
  info("ndof = %d.", ndof);

  // Define projection norms.
  ProjNormType vel_proj_norm = HERMES_H1_NORM;
#ifdef PRESSURE_IN_L2
  ProjNormType p_proj_norm = HERMES_L2_NORM;
#else
  ProjNormType p_proj_norm = HERMES_H1_NORM;
#endif

  // Solutions for the Newton's iteration and time stepping.
  info("Setting initial conditions.");
  ZeroSolution<double> xvel_prev_time(&mesh), yvel_prev_time(&mesh), p_prev_time(&mesh);

  // Initialize weak formulation.
  WeakForm<double>* wf = new WeakFormNSNewton(STOKES, RE, TAU, &xvel_prev_time, &yvel_prev_time);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(wf, Hermes::vector<const Space<double> *>(&xvel_space, &yvel_space, &p_space));

  // Initialize the Newton solver.
  Hermes::Hermes2D::NewtonSolver<double> newton(&dp, matrix_solver_type);

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  double* coeff_vec = new double[Space<double>::get_num_dofs(Hermes::vector<const Space<double> *>(&xvel_space, &yvel_space, &p_space))];
  info("Projecting initial condition to obtain initial vector for the Newton's method.");
  OGProjection<double>::project_global(Hermes::vector<const Space<double> *>(&xvel_space, &yvel_space, &p_space),
    Hermes::vector<MeshFunction<double> *>(&xvel_prev_time, &yvel_prev_time, &p_prev_time),
    coeff_vec, matrix_solver_type,
    Hermes::vector<ProjNormType>(vel_proj_norm, vel_proj_norm, p_proj_norm));

  // Time-stepping loop:
  int num_time_steps = T_FINAL / TAU;
  for (int ts = 1; ts <= num_time_steps; ts++)
  {
    current_time += TAU;
    info("---- Time step %d, time = %g:", ts, current_time);

    // Update time-dependent essential BCs.
    if (current_time <= STARTUP_TIME)
    {
      info("Updating time-dependent essential BC.");
      Space<double>::update_essential_bc_values(Hermes::vector<Space<double> *>(&xvel_space, &yvel_space, &p_space), current_time);
    }

    // Perform Newton's iteration.
    bool verbose = true;
    // Perform Newton's iteration and translate the resulting coefficient vector into previous time level solutions.
    newton.set_verbose_output(verbose);
    try{
      newton.solve(coeff_vec, NEWTON_TOL, NEWTON_MAX_ITER);
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.printMsg();
      error("Newton's iteration failed.");
    }
    Hermes::vector<Solution<double> *> tmp(&xvel_prev_time, &yvel_prev_time, &p_prev_time);
    Hermes::Hermes2D::Solution<double>::vector_to_solutions(newton.get_sln_vector(), Hermes::vector<const Space<double> *>(&xvel_space, &yvel_space, &p_space), tmp);
  }

  delete [] coeff_vec;

  info("Coordinate (   0, 2.5) xvel value = %lf", xvel_prev_time.get_pt_value(0.0, 2.5));
  info("Coordinate (   5, 2.5) xvel value = %lf", xvel_prev_time.get_pt_value(5.0, 2.5));
  info("Coordinate ( 7.5, 2.5) xvel value = %lf", xvel_prev_time.get_pt_value(7.5, 2.5));
  info("Coordinate (  10, 2.5) xvel value = %lf", xvel_prev_time.get_pt_value(10.0, 2.5));
  info("Coordinate (12.5, 2.5) xvel value = %lf", xvel_prev_time.get_pt_value(12.5, 2.5));
  info("Coordinate (  15, 2.5) xvel value = %lf", xvel_prev_time.get_pt_value(15.0, 2.5));

  info("Coordinate (   0, 2.5) yvel value = %lf", yvel_prev_time.get_pt_value(0.0, 2.5));
  info("Coordinate (   5, 2.5) yvel value = %lf", yvel_prev_time.get_pt_value(5.0, 2.5));
  info("Coordinate ( 7.5, 2.5) yvel value = %lf", yvel_prev_time.get_pt_value(7.5, 2.5));
  info("Coordinate (  10, 2.5) yvel value = %lf", yvel_prev_time.get_pt_value(10.0, 2.5));
  info("Coordinate (12.5, 2.5) yvel value = %lf", yvel_prev_time.get_pt_value(12.5, 2.5));
  info("Coordinate (  15, 2.5) yvel value = %lf", yvel_prev_time.get_pt_value(15.0, 2.5));

  int success = 1;
  double eps = 1e-5;
  if (fabs(xvel_prev_time.get_pt_value(0.0, 2.5) - 0.200000) > eps) {
    printf("Coordinate (   0, 2.5) xvel value is %g\n", xvel_prev_time.get_pt_value(0.0, 2.5));
    success = 0;
  }
  if (fabs(xvel_prev_time.get_pt_value(5, 2.5) - 0.134291) > eps) {
    printf("Coordinate (   5, 2.5) xvel value is %g\n", xvel_prev_time.get_pt_value(5, 2.5));
    success = 0;
  }
  if (fabs(xvel_prev_time.get_pt_value(7.5, 2.5) - 0.135088) > eps) {
    printf("Coordinate ( 7.5, 2.5) xvel value is %g\n", xvel_prev_time.get_pt_value(7.5, 2.5));
    success = 0;
  }
  if (fabs(xvel_prev_time.get_pt_value(10, 2.5) - 0.134944) > eps) {
    printf("Coordinate (  10, 2.5) xvel value is %g\n", xvel_prev_time.get_pt_value(10, 2.5));
    success = 0;
  }
  if (fabs(xvel_prev_time.get_pt_value(12.5, 2.5) - 0.134888) > eps) {
    printf("Coordinate (12.5, 2.5) xvel value is %g\n", xvel_prev_time.get_pt_value(12.5, 2.5));
    success = 0;
  }
  if (fabs(xvel_prev_time.get_pt_value(15, 2.5) - 0.134864) > eps) {
    printf("Coordinate (  15, 2.5) xvel value is %g\n", xvel_prev_time.get_pt_value(15, 2.5));
    success = 0;
  }

  if (fabs(yvel_prev_time.get_pt_value(0.0, 2.5) - 0.000000) > eps) {
    printf("Coordinate (   0, 2.5) yvel value is %g\n", yvel_prev_time.get_pt_value(0.0, 2.5));
    success = 0;
  }
  if (fabs(yvel_prev_time.get_pt_value(5, 2.5) - 0.000493) > eps) {
    printf("Coordinate (   5, 2.5) yvel value is %g\n", yvel_prev_time.get_pt_value(5, 2.5));
    success = 0;
  }
  if (fabs(yvel_prev_time.get_pt_value(7.5, 2.5) - 0.000070) > eps) {
    printf("Coordinate ( 7.5, 2.5) yvel value is %g\n", yvel_prev_time.get_pt_value(7.5, 2.5));
    success = 0;
  }
  if (fabs(yvel_prev_time.get_pt_value(10, 2.5) - 0.000008) > eps) {
    printf("Coordinate (  10, 2.5) yvel value is %g\n", yvel_prev_time.get_pt_value(10, 2.5));
    success = 0;
  }
  if (fabs(yvel_prev_time.get_pt_value(12.5, 2.5) + 0.000003) > eps) {
    printf("Coordinate (12.5, 2.5) yvel value is %g\n", yvel_prev_time.get_pt_value(12.5, 2.5));
    success = 0;
  }
  if (fabs(yvel_prev_time.get_pt_value(15, 2.5) + 0.000006) > eps) {
    printf("Coordinate (  15, 2.5) yvel value is %g\n", yvel_prev_time.get_pt_value(15, 2.5));
    success = 0;
  }

  if (success == 1) {
    printf("Success!\n");
    return TEST_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return TEST_FAILURE;
  }
}
