#include "hermes2d.h"

#define USE_PARALUTION

using namespace Hermes;
using namespace Hermes::Hermes2D;

// The time-dependent laminar incompressible Navier-Stokes equations are
// discretized in time via the implicit Euler method. If NEWTON == true,
// the Newton's method is used to solve the nonlinear problem at each time
// step. If NEWTON == false, the convective term is only linearized using the
// velocities from the previous time step. Obviously the latter approach is wrong,
// but people do this frequently because it is faster and simpler to implement.
// Therefore we include this case for comparison purposes. We also show how
// to use discontinuous ($L^2$) elements for pressure and thus make the
// velocity discreetely divergence free. Comparison to approximating the
// pressure with the standard (continuous) Taylor-Hood elements is enabled.
// The Reynolds number Re = 200 which is embarrassingly low. You
// can increase it but then you will need to make the mesh finer, and the
// computation will take more time.
//
// PDE: incompressible Navier-Stokes equations in the form
// \partial v / \partial t - \Delta v / Re + (v \cdot \nabla) v + \nabla p = 0,
// div v = 0.
//
// BC: u_1 is a time-dependent constant and u_2 = 0 on Gamma_4 (inlet),
//     u_1 = u_2 = 0 on Gamma_1 (bottom), Gamma_3 (top) and Gamma_5 (obstacle),
//     "do nothing" on Gamma_2 (outlet).
//
// Geometry: Rectangular channel containing an off-axis circular obstacle. The
//           radius and position of the circle, as well as other geometry
//           parameters can be changed in the mesh file "domain.mesh".
//
// The following parameters can be changed:

// For application of Stokes flow (creeping flow).
// If this is defined, the pressure is approximated using
// discontinuous L2 elements (making the velocity discreetely
// divergence-free, more accurate than using a continuous
// pressure approximation). Otherwise the standard continuous
// elements are used. The results are striking - check the
// tutorial for comparisons.
const bool STOKES = false;

const bool HERMES_VISUALIZATION = true;

#define PRESSURE_IN_L2

// Initial polynomial degree for velocity components.
const int P_INIT_VEL = 2;

// Initial polynomial degree for pressure.
// Note: P_INIT_VEL should always be greater than
// P_INIT_PRESSURE because of the inf-sup condition.
const int P_INIT_PRESSURE = 1;

// Reynolds number.
const double RE = 200.0;

// Inlet velocity (reached after STARTUP_TIME).
const double VEL_INLET = 1.0;

// Current time (used in weak forms).
double current_time = 0;

// During this time, inlet velocity increases gradually
// from 0 to VEL_INLET, then it stays constant.
const double STARTUP_TIME = 1.0;

const double TAU = 0.1;                           // Time step.
const double T_FINAL = 100.0;                      // Time interval length.
const double NEWTON_TOL = 1e-3;                   // Stopping criterion for the Newton's method.
const double H = 5;                               // Domain height (necessary to define the parabolic velocity profile at inlet).

// Boundary markers.
const std::string BDY_BOTTOM = "1";
const std::string BDY_RIGHT = "2";
const std::string BDY_TOP = "3";
const std::string BDY_LEFT = "4";
const std::string BDY_OBSTACLE = "5";

// Weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  // This either enables PARALUTION or leaves the default UMFPACK.
#ifdef USE_PARALUTION
  HermesCommonApi.set_integral_param_value(matrixSolverType, SOLVER_PARALUTION_ITERATIVE);
#endif

#include "setup.cpp"

  // Initialize the Newton solver.
  Hermes::Hermes2D::NewtonSolver<double> newton(&wf, spaces);
  newton.output_matrix();
  newton.set_matrix_export_format(Hermes::Algebra::MatrixExportFormat::EXPORT_FORMAT_MATRIX_MARKET);
  newton.set_rhs_export_format(Hermes::Algebra::MatrixExportFormat::EXPORT_FORMAT_MATRIX_MARKET);
  newton.output_rhs();
  // Verbose.
  newton.get_linear_matrix_solver()->set_verbose_output(true);

#ifdef USE_PARALUTION
  // Relative tolerance of PARALUTION.
  newton.get_linear_matrix_solver()->as_IterSolver()->set_tolerance(1e-8, Solvers::LoopSolverToleranceType::RelativeTolerance);
  // Use GMRES.
  newton.get_linear_matrix_solver()->as_IterSolver()->set_solver_type(Solvers::IterSolverType::GMRES);
  // Use Saddle-point preconditioner.
  newton.get_linear_matrix_solver()->as_IterSolver()->set_precond(new Preconditioners::ParalutionPrecond<double>(Preconditioners::PreconditionerType::SaddlePoint));
#endif
 
  // Newton method setup:
  // - max allowed iterations
  newton.set_max_allowed_iterations(10);
  // - no damping
  newton.set_manual_damping_coeff(true, 1.0);
  // - nonlinear tolerance (absolute)
  newton.set_tolerance(1e-3, Hermes::Solvers::ResidualNormAbsolute);

  // Time-stepping loop:
  for (int ts = 1; ts <= T_FINAL / TAU; ts++)
  {
    current_time += TAU;
    Hermes::Mixins::Loggable::Static::info("Time step %i, time %f.", ts, current_time);

    // Update time-dependent essential BCs.
    newton.set_time(current_time);

    // Solve Newton.
    newton.solve(coeff_vec);

    // Get the solutions.
    Hermes::Hermes2D::Solution<double>::vector_to_solutions(newton.get_sln_vector(), spaces, sln_prev_time);

    // Visualization.
    vview.set_title("Velocity, time %g", current_time);
    vview.show(xvel_prev_time, yvel_prev_time);
    pview.set_title("Pressure, time %g", current_time);
    pview.get_linearizer()->set_criterion(Views::LinearizerCriterionFixed(1));
    pview.show(p_prev_time);
  }

  return 0;
}
