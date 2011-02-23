#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This test makes sure that example 21-newton-timedep-ns works correctly.

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
const double TAU = 0.1;                           // Time step.
const double T_FINAL = 0.21;                      // Time interval length.
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
double TIME = 0;

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
  BCValues bc_values_x(&TIME);
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

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  scalar* coeff_vec = new scalar[Space::get_num_dofs(Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space))];
  if (NEWTON) {
    info("Projecting initial condition to obtain initial vector for the Newton's method.");
    OGProjection::project_global(Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space), 
                   Hermes::vector<MeshFunction *>(&xvel_prev_time, &yvel_prev_time, &p_prev_time), 
                   coeff_vec, 
                   matrix_solver, 
                   Hermes::vector<ProjNormType>(vel_proj_norm, vel_proj_norm, p_proj_norm));
  }

  // Time-stepping loop:
  char title[100];
  int num_time_steps = T_FINAL / TAU;
  for (int ts = 1; ts <= num_time_steps; ts++)
  {
    TIME += TAU;
    info("---- Time step %d, time = %g:", ts, TIME);

    // Update time-dependent essential BC are used.
    if (TIME <= STARTUP_TIME) {
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
 }

  delete [] coeff_vec;
  delete matrix;
  delete rhs;
  delete solver;

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
  if (fabs(xvel_prev_time.get_pt_value(5, 2.5) - 0.130866) > eps) {
    printf("Coordinate (   5, 2.5) xvel value is %g\n", xvel_prev_time.get_pt_value(5, 2.5));
    success = 0;
  }
  if (fabs(xvel_prev_time.get_pt_value(7.5, 2.5) - 0.134637) > eps) {
    printf("Coordinate ( 7.5, 2.5) xvel value is %g\n", xvel_prev_time.get_pt_value(7.5, 2.5));
    success = 0;
  }
  if (fabs(xvel_prev_time.get_pt_value(10, 2.5) - 0.134801) > eps) {
    printf("Coordinate (  10, 2.5) xvel value is %g\n", xvel_prev_time.get_pt_value(10, 2.5));
    success = 0;
  }
  if (fabs(xvel_prev_time.get_pt_value(12.5, 2.5) - 0.134826) > eps) {
    printf("Coordinate (12.5, 2.5) xvel value is %g\n", xvel_prev_time.get_pt_value(12.5, 2.5));
    success = 0;
  }
  if (fabs(xvel_prev_time.get_pt_value(15, 2.5) - 0.134840) > eps) {
    printf("Coordinate (  15, 2.5) xvel value is %g\n", xvel_prev_time.get_pt_value(15, 2.5));
    success = 0;
  }

  if (fabs(yvel_prev_time.get_pt_value(0.0, 2.5) - 0.000000) > eps) {
    printf("Coordinate (   0, 2.5) yvel value is %g\n", yvel_prev_time.get_pt_value(0.0, 2.5));
    success = 0;
  }
  if (fabs(yvel_prev_time.get_pt_value(5, 2.5) - 0.000584) > eps) {
    printf("Coordinate (   5, 2.5) yvel value is %g\n", yvel_prev_time.get_pt_value(5, 2.5));
    success = 0;
  }
  if (fabs(yvel_prev_time.get_pt_value(7.5, 2.5) - 0.000101) > eps) {
    printf("Coordinate ( 7.5, 2.5) yvel value is %g\n", yvel_prev_time.get_pt_value(7.5, 2.5));
    success = 0;
  }
  if (fabs(yvel_prev_time.get_pt_value(10, 2.5) - 0.000029) > eps) {
    printf("Coordinate (  10, 2.5) yvel value is %g\n", yvel_prev_time.get_pt_value(10, 2.5));
    success = 0;
  }
  if (fabs(yvel_prev_time.get_pt_value(12.5, 2.5) - 0.000013) > eps) {
    printf("Coordinate (12.5, 2.5) yvel value is %g\n", yvel_prev_time.get_pt_value(12.5, 2.5));
    success = 0;
  }
  if (fabs(yvel_prev_time.get_pt_value(15, 2.5) - 0.000009) > eps) {
    printf("Coordinate (  15, 2.5) yvel value is %g\n", yvel_prev_time.get_pt_value(15, 2.5));
    success = 0;
  }

  if (success == 1) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}
