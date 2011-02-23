#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

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

#define PRESSURE_IN_L2                            // If this is defined, the pressure is approximated using
                                                  // discontinuous L2 elements (making the velocity discreetely
                                                  // divergence-free, more accurate than using a continuous
                                                  // pressure approximation). Otherwise the standard continuous
                                                  // elements are used. The results are striking - check the
                                                  // tutorial for comparisons.
const bool NEWTON = true;                         // If NEWTON == true then the Newton's iteration is performed.
                                                  // in every time step. Otherwise the convective term is linearized
                                                  // using the velocities from the previous time step.
const int P_INIT_VEL = 2;                         // Initial polynomial degree for velocity components.
const int P_INIT_PRESSURE = 1;                    // Initial polynomial degree for pressure.
                                                  // Note: P_INIT_VEL should always be greater than
                                                  // P_INIT_PRESSURE because of the inf-sup condition.
const double RE = 200.0;                          // Reynolds number.
const double VEL_INLET = 1.0;                     // Inlet velocity (reached after STARTUP_TIME).
const double STARTUP_TIME = 1.0;                  // During this time, inlet velocity increases gradually
                                                  // from 0 to VEL_INLET, then it stays constant.
const double time_step = 0.1;                     // Time step.
const double T_FINAL = 30000.0;                   // Time interval length.
const double NEWTON_TOL = 1e-3;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 10;                   // Maximum allowed number of Newton iterations.
const double H = 5;                               // Domain height (necessary to define the parabolic
                                                  // velocity profile at inlet).
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Boundary markers.
const int BDY_BOTTOM = 1;
const int BDY_RIGHT = 2;
const int BDY_TOP = 3;
const int BDY_LEFT = 4;
const int BDY_OBSTACLE = 5;

// Current time (used in weak forms).
double current_time = 0;

// Essential (Dirichlet) boundary condition values for x-velocity.
scalar essential_bc_values_xvel(double x, double y, double time) {
  // Time-dependent parabolic profile at inlet.
  double val_y = VEL_INLET * y*(H-y) / (H/2.)/(H/2.); // Peak value VEL_INLET at y = H/2.
  if (time <= STARTUP_TIME) return val_y * time/STARTUP_TIME;
  else return val_y;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Initial mesh refinements.
  mesh.refine_all_elements();
  mesh.refine_towards_boundary(BDY_OBSTACLE, 4, false);
  mesh.refine_towards_boundary(BDY_TOP, 4, true);     // '4' is the number of levels,
  mesh.refine_towards_boundary(BDY_BOTTOM, 4, true);  // 'true' stands for anisotropic refinements.

  // Boundary condition types for x-velocity and y-velocity.
  BCTypes xvel_bc_type, yvel_bc_type;
  xvel_bc_type.add_bc_none(BDY_RIGHT);
  xvel_bc_type.add_bc_dirichlet(Hermes::vector<int>(BDY_BOTTOM, BDY_TOP, BDY_LEFT, BDY_OBSTACLE));
  yvel_bc_type.add_bc_none(BDY_RIGHT);
  yvel_bc_type.add_bc_dirichlet(Hermes::vector<int>(BDY_BOTTOM, BDY_TOP, BDY_LEFT, BDY_OBSTACLE));

  // Enter Dirichlet boundary values.
  BCValues bc_values_x(&current_time);
  bc_values_x.add_timedep_function(BDY_LEFT, essential_bc_values_xvel);
  bc_values_x.add_zero(Hermes::vector<int>(BDY_BOTTOM, BDY_TOP, BDY_OBSTACLE));
  
  BCValues bc_values_y;
  bc_values_y.add_zero(Hermes::vector<int>(BDY_BOTTOM, BDY_TOP, BDY_LEFT, BDY_OBSTACLE));

  // Spaces for velocity components and pressure.
  H1Space xvel_space(&mesh, &xvel_bc_type, &bc_values_x, P_INIT_VEL);
  H1Space yvel_space(&mesh, &yvel_bc_type, &bc_values_y, P_INIT_VEL);
#ifdef PRESSURE_IN_L2
  L2Space p_space(&mesh, P_INIT_PRESSURE);
#else
  H1Space p_space(&mesh, NULL, NULL, P_INIT_PRESSURE);
#endif

  // Calculate and report the number of degrees of freedom.
  int ndof = Space::get_num_dofs(Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space));
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
  Solution xvel_prev_time, yvel_prev_time, p_prev_time; 
  xvel_prev_time.set_zero(&mesh);
  yvel_prev_time.set_zero(&mesh);
  p_prev_time.set_zero(&mesh);

  // Initialize weak formulation.
  WeakForm wf(3);
  if (NEWTON) {
    wf.add_matrix_form(0, 0, callback(bilinear_form_sym_0_0_1_1), HERMES_SYM);
    wf.add_matrix_form(0, 0, callback(newton_bilinear_form_unsym_0_0), HERMES_NONSYM, HERMES_ANY);
    wf.add_matrix_form(0, 1, callback(newton_bilinear_form_unsym_0_1), HERMES_NONSYM, HERMES_ANY);
    wf.add_matrix_form(0, 2, callback(bilinear_form_unsym_0_2), HERMES_ANTISYM);
    wf.add_matrix_form(1, 0, callback(newton_bilinear_form_unsym_1_0), HERMES_NONSYM, HERMES_ANY);
    wf.add_matrix_form(1, 1, callback(bilinear_form_sym_0_0_1_1), HERMES_SYM);
    wf.add_matrix_form(1, 1, callback(newton_bilinear_form_unsym_1_1), HERMES_NONSYM, HERMES_ANY);
    wf.add_matrix_form(1, 2, callback(bilinear_form_unsym_1_2), HERMES_ANTISYM);
    wf.add_vector_form(0, callback(newton_F_0), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&xvel_prev_time, &yvel_prev_time));
    wf.add_vector_form(1, callback(newton_F_1), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&xvel_prev_time, &yvel_prev_time));
    wf.add_vector_form(2, callback(newton_F_2), HERMES_ANY);
  }
  else {
    wf.add_matrix_form(0, 0, callback(bilinear_form_sym_0_0_1_1), HERMES_SYM);
    wf.add_matrix_form(0, 0, callback(simple_bilinear_form_unsym_0_0_1_1), 
                  HERMES_NONSYM, HERMES_ANY, Hermes::vector<MeshFunction*>(&xvel_prev_time, &yvel_prev_time));
    wf.add_matrix_form(1, 1, callback(bilinear_form_sym_0_0_1_1), HERMES_SYM);
    wf.add_matrix_form(1, 1, callback(simple_bilinear_form_unsym_0_0_1_1), 
                  HERMES_NONSYM, HERMES_ANY, Hermes::vector<MeshFunction*>(&xvel_prev_time, &yvel_prev_time));
    wf.add_matrix_form(0, 2, callback(bilinear_form_unsym_0_2), HERMES_ANTISYM);
    wf.add_matrix_form(1, 2, callback(bilinear_form_unsym_1_2), HERMES_ANTISYM);
    wf.add_vector_form(0, callback(simple_linear_form), HERMES_ANY, &xvel_prev_time);
    wf.add_vector_form(1, callback(simple_linear_form), HERMES_ANY, &yvel_prev_time);
  }

  // Initialize the FE problem.
  bool is_linear;
  if (NEWTON) is_linear = false;
  else is_linear = true;
  DiscreteProblem dp(&wf, Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space), is_linear);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Initialize views.
  VectorView vview("velocity [m/s]", new WinGeom(0, 0, 750, 240));
  ScalarView pview("pressure [Pa]", new WinGeom(0, 290, 750, 240));
  vview.set_min_max_range(0, 1.6);
  vview.fix_scale_width(80);
  //pview.set_min_max_range(-0.9, 1.0);
  pview.fix_scale_width(80);
  pview.show_mesh(true);

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  scalar* coeff_vec = new scalar[Space::get_num_dofs(Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space))];
  if (NEWTON) {
    info("Projecting initial condition to obtain initial vector for the Newton's method.");
    OGProjection::project_global(Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space), 
                   Hermes::vector<MeshFunction *>(&xvel_prev_time, &yvel_prev_time, &p_prev_time), 
                   coeff_vec, matrix_solver, 
                   Hermes::vector<ProjNormType>(vel_proj_norm, vel_proj_norm, p_proj_norm));
  }

  // Time-stepping loop:
  char title[100];
  int num_time_steps = T_FINAL / time_step;
  for (int ts = 1; ts <= num_time_steps; ts++)
  {
    current_time += time_step;
    info("---- Time step %d, time = %g:", ts, current_time);

    // Update time-dependent essential BC are used.
    if (current_time <= STARTUP_TIME) {
      info("Updating time-dependent essential BC.");
      update_essential_bc_values(Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space));
    }

    if (NEWTON) 
    {
      // Perform Newton's iteration.
      info("Solving nonlinear problem:");
      bool verbose = true;
      if (!solve_newton(coeff_vec, &dp, solver, matrix, rhs, 
          NEWTON_TOL, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

      // Update previous time level solutions.
      Solution::vector_to_solutions(coeff_vec, Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space), 
                                    Hermes::vector<Solution *>(&xvel_prev_time, &yvel_prev_time, &p_prev_time));
    }
    else {
      // Linear solve.
      info("Assembling and solving linear problem.");
      dp.assemble(matrix, rhs, false);
      if(solver->solve()) 
        Solution::vector_to_solutions(solver->get_solution(), 
                  Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space), 
                  Hermes::vector<Solution *>(&xvel_prev_time, &yvel_prev_time, &p_prev_time));
      else 
        error ("Matrix solver failed.\n");
    }

    // Show the solution at the end of time step.
    sprintf(title, "Velocity, time %g", current_time);
    vview.set_title(title);
    vview.show(&xvel_prev_time, &yvel_prev_time, HERMES_EPS_LOW);
    sprintf(title, "Pressure, time %g", current_time);
    pview.set_title(title);
    pview.show(&p_prev_time);
 }

  delete [] coeff_vec;
  delete matrix;
  delete rhs;
  delete solver;
  
  // Wait for all views to be closed.
  View::wait();
  return 0;
}
