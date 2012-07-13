#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"

#include "hermes2d.h"

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

// During this time, inlet velocity increases gradually
// from 0 to VEL_INLET, then it stays constant.
const double STARTUP_TIME = 1.0;

const double TAU = 0.1;                           // Time step.
const double T_FINAL = 30000.0;                   // Time interval length.
const double NEWTON_TOL = 1e-3;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 10;                   // Maximum allowed number of Newton iterations.
const double H = 5;                               // Domain height (necessary to define the parabolic
// velocity profile at inlet).

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
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  Hermes::Hermes2D::Hermes2DApi.setParamValue(Hermes::Hermes2D::numThreads, 8);
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", &mesh);

  // Initial mesh refinements.
  mesh.refine_towards_boundary(BDY_OBSTACLE, 1, false);
  mesh.refine_towards_boundary(BDY_TOP, 1, true);     // '4' is the number of levels,
  mesh.refine_towards_boundary(BDY_BOTTOM, 1, true);  // 'true' stands for anisotropic refinements.

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

  // Define projection norms.
  ProjNormType vel_proj_norm = HERMES_H1_NORM;
#ifdef PRESSURE_IN_L2
  ProjNormType p_proj_norm = HERMES_L2_NORM;
#else
  ProjNormType p_proj_norm = HERMES_H1_NORM;
#endif

  // Solutions for the Newton's iteration and time stepping.
  ZeroSolution<double> xvel_prev_time(&mesh), yvel_prev_time(&mesh), p_prev_time(&mesh);

  // Initialize weak formulation.
  WeakForm<double>* wf;
  wf = new WeakFormNSNewton(STOKES, RE, TAU, &xvel_prev_time, &yvel_prev_time);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(wf, Hermes::vector<const Space<double> *>(&xvel_space, &yvel_space, &p_space));

  // Initialize the Newton solver.
  Hermes::Hermes2D::NewtonSolver<double> newton(&dp);

  // Initialize views.
  Views::VectorView vview("velocity[m/s]", new Views::WinGeom(0, 0, 750, 240));
  Views::ScalarView pview("pressure[Pa]", new Views::WinGeom(0, 290, 750, 240));
  vview.set_min_max_range(0, 1.6);
  vview.fix_scale_width(80);
  //pview.set_min_max_range(-0.9, 1.0);
  pview.fix_scale_width(80);
if(HERMES_VISUALIZATION)
  pview.show_mesh(true);

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  double* coeff_vec = new double[Space<double>::get_num_dofs(Hermes::vector<const Space<double> *>(&xvel_space, &yvel_space, &p_space))];
  OGProjection<double> ogProjection;

  ogProjection.project_global(Hermes::vector<const Space<double> *>(&xvel_space, &yvel_space, &p_space),
    Hermes::vector<MeshFunction<double> *>(&xvel_prev_time, &yvel_prev_time, &p_prev_time),
    coeff_vec, Hermes::vector<ProjNormType>(vel_proj_norm, vel_proj_norm, p_proj_norm));

  // Time-stepping loop:
  char title[100];
  int num_time_steps = T_FINAL / TAU;
  for (int ts = 1; ts <= num_time_steps; ts++)
  {
    std::cout << ts << std::endl;
    current_time += TAU;

    // Update time-dependent essential BCs.
    if(current_time <= STARTUP_TIME)
    {
      Space<double>::update_essential_bc_values(Hermes::vector<Space<double> *>(&xvel_space, &yvel_space, &p_space), current_time);
    }

    // Perform Newton's iteration and translate the resulting coefficient vector into previous time level solutions.
    try{
      newton.solve(coeff_vec, NEWTON_TOL, NEWTON_MAX_ITER);
    }
    catch(Hermes::Exceptions::Exception& e)
    {
      e.printMsg();
    }
    Hermes::vector<Solution<double> *> tmp(&xvel_prev_time, &yvel_prev_time, &p_prev_time);
    Hermes::Hermes2D::Solution<double>::vector_to_solutions(newton.get_sln_vector(), Hermes::vector<const Space<double> *>(&xvel_space, &yvel_space, &p_space), tmp);

    // Show the solution at the end of time step.
    if(HERMES_VISUALIZATION)
    {
      sprintf(title, "Velocity, time %g", current_time);
      vview.set_title(title);
      vview.show(&xvel_prev_time, &yvel_prev_time);
      sprintf(title, "Pressure, time %g", current_time);
      pview.set_title(title);
      pview.show(&p_prev_time);
    }
  }

  delete [] coeff_vec;

  // Wait for all views to be closed.
  if(HERMES_VISUALIZATION)
    Views::View::wait();
  return 0;
}