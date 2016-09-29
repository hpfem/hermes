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

// Current time (used in weak forms).
double current_time = 0;

// During this time, inlet velocity increases gradually
// from 0 to VEL_INLET, then it stays constant.
const double STARTUP_TIME = 1.0;

// Time step.
const double TAU = 0.1;
// Time interval length.
const double T_FINAL = 1000.0;
// Stopping criterion for the Newton's method.
const double NEWTON_TOL = 1e-3;
// Domain height (necessary to define the parabolic velocity profile at inlet).
const double H = 5;

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
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", mesh);

  // Initial mesh refinements.
  mesh->refine_towards_boundary(BDY_OBSTACLE, 2, false);
  mesh->refine_towards_boundary(BDY_TOP, 2, true);     // '4' is the number of levels,
  mesh->refine_towards_boundary(BDY_BOTTOM, 2, true);  // 'true' stands for anisotropic refinements.
  mesh->refine_all_elements();
  mesh->refine_all_elements();

  // Initialize boundary conditions.
  EssentialBCNonConst bc_left_vel_x(BDY_LEFT, VEL_INLET, H, STARTUP_TIME);
  Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_other_vel_x({ BDY_BOTTOM, BDY_TOP, BDY_OBSTACLE }, 0.0);
  Hermes::Hermes2D::EssentialBCs<double> bcs_vel_x({ &bc_left_vel_x, &bc_other_vel_x });
  Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_vel_y({ BDY_LEFT, BDY_BOTTOM, BDY_TOP, BDY_OBSTACLE }, 0.0);
  Hermes::Hermes2D::EssentialBCs<double> bcs_vel_y(&bc_vel_y);
  Hermes::Hermes2D::EssentialBCs<double> bcs_pressure;

  // Spaces for velocity components and pressure.
  SpaceSharedPtr<double> xvel_space(new H1Space<double>(mesh, &bcs_vel_x, P_INIT_VEL));
  SpaceSharedPtr<double> yvel_space(new H1Space<double>(mesh, &bcs_vel_y, P_INIT_VEL));
#ifdef PRESSURE_IN_L2
  SpaceSharedPtr<double> p_space(new L2Space<double>(mesh, P_INIT_PRESSURE));
#else
  SpaceSharedPtr<double> p_space(new H1Space<double>(mesh, &bcs_pressure, P_INIT_PRESSURE));
#endif
  std::vector<SpaceSharedPtr<double> > spaces({ xvel_space, yvel_space, p_space });

  // Calculate and report the number of degrees of freedom.
  int ndof = Space<double>::get_num_dofs(spaces);

  // Define projection norms.
  NormType vel_proj_norm = HERMES_H1_NORM;
#ifdef PRESSURE_IN_L2
  NormType p_proj_norm = HERMES_L2_NORM;
#else
  NormType p_proj_norm = HERMES_H1_NORM;
#endif
  std::vector<NormType> proj_norms({ vel_proj_norm, vel_proj_norm, p_proj_norm });

  // Solutions for the Newton's iteration and time stepping.
  MeshFunctionSharedPtr<double> xvel_prev_time(new ConstantSolution<double>(mesh, 0.0));
  MeshFunctionSharedPtr<double> yvel_prev_time(new ConstantSolution<double>(mesh, 0.0));
  MeshFunctionSharedPtr<double> p_prev_time(new ConstantSolution<double>(mesh, 0.0));
  std::vector<MeshFunctionSharedPtr<double> > sln_prev_time = { xvel_prev_time, yvel_prev_time, p_prev_time };
  double* coeff_vec = new double[Space<double>::get_num_dofs(spaces)];

  // Project the initial condition on the FE space to obtain initial coefficient vector for the Newton's method.
  OGProjection<double>::project_global(spaces, sln_prev_time, coeff_vec, proj_norms);

  // Initialize weak formulation.
  WeakFormSharedPtr<double> wf(new WeakFormNSNewton(STOKES, RE, TAU, xvel_prev_time, yvel_prev_time));
  UExtFunctionSharedPtr<double> fn_0(new CustomUExtFunction(0));
  UExtFunctionSharedPtr<double> fn_1(new CustomUExtFunction(1));
  wf->set_ext({ xvel_prev_time, yvel_prev_time });
  wf->set_u_ext_fn({ fn_0, fn_1 });

  // Initialize views.
  Views::VectorView vview("velocity[m/s]", new Views::WinGeom(0, 0, 750, 240));
  Views::ScalarView pview("pressure[Pa]", new Views::WinGeom(0, 290, 750, 240));
  vview.set_min_max_range(0, 1.6);
  vview.fix_scale_width(80);
  //pview.set_min_max_range(-0.9, 1.0);
  pview.fix_scale_width(80);
  pview.show_mesh(true);

  // Initialize the Newton solver.
  Hermes::Hermes2D::NewtonSolver<double> newton(wf, spaces);
  newton.get_linear_matrix_solver()->set_verbose_output(true);
  //newton.output_matrix();
  newton.set_matrix_export_format(Hermes::Algebra::MatrixExportFormat::EXPORT_FORMAT_MATRIX_MARKET);
  //newton.output_rhs();
  newton.set_rhs_export_format(Hermes::Algebra::MatrixExportFormat::EXPORT_FORMAT_MATRIX_MARKET);

  // Newton method setup:
  // - max allowed iterations
  newton.set_max_allowed_iterations(10);
  // - no damping
  newton.set_manual_damping_coeff(true, 1.0);
  // - nonlinear tolerance (absolute)
  newton.set_tolerance(1e-3, Hermes::Solvers::ResidualNormAbsolute);

  // Time-stepping loop:
  for (int time_step = 1; time_step <= T_FINAL / TAU; time_step++, current_time += TAU)
  {
    // Update time-dependent essential BCs.
    newton.set_time(current_time);

    // Solve Newton.
    newton.solve(coeff_vec);

    // Get the solutions.
    Hermes::Hermes2D::Solution<double>::vector_to_solutions(newton.get_sln_vector(), spaces, sln_prev_time);

    // Visualization.
    vview.set_title("Velocity, time %g", current_time);
    vview.show(xvel_prev_time, yvel_prev_time);
    if(!(time_step % 100))
      vview.save_numbered_screenshot("Velocity", (time_step / 100), true);

    pview.set_title("Pressure, time %g", current_time);
    pview.show(p_prev_time);
    if (!(time_step % 100))
      pview.save_numbered_screenshot("Pressure", (time_step / 100), true);
  }

  return 0;
}