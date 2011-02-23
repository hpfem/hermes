#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace Teuchos;
using namespace RefinementSelectors;

// This test makes sure that example 42-trilinos-timedep works correctly.

const int INIT_REF_NUM = 4;       // Number of initial uniform mesh refinements.
const int P_INIT = 1;             // Initial polynomial degree of all mesh elements.
const double ALPHA = 10.0;        // Coefficient for the Nwwton boundary condition.
const double LAMBDA = 1e5;
const double HEATCAP = 1e6;
const double RHO = 3000.0;
const double TEMP_EXT = 20.0;
const double TEMP_INIT = 10.0;

const double TAU = 50.0;          // Time step.        

const bool JFNK = true;
const bool PRECOND = true;

// Boundary markers.
const int BDY_BOTTOM = 1, BDY_RIGHT = 2, BDY_TOP = 3, BDY_LEFT = 4;

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Perform initial mesh refinemets.
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_BOTTOM);
  bc_types.add_bc_newton(Hermes::vector<int>(BDY_RIGHT, BDY_TOP, BDY_LEFT));

  // Enter Dirichlet boundary values.
  BCValues bc_values;
  bc_values.add_const(BDY_BOTTOM, TEMP_INIT);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof: %d", ndof);

  // Define constant initial condition. 
  Solution t_prev_time;
  t_prev_time.set_const(&mesh, TEMP_INIT);

  // Initialize the weak formulation.
  WeakForm wf(1, JFNK ? true : false);
  wf.add_matrix_form(callback(jacobian));
  wf.add_matrix_form_surf(callback(jacobian_surf));
  wf.add_vector_form(callback(residual), HERMES_ANY, &t_prev_time);
  wf.add_vector_form_surf(callback(residual_surf));

  // Initialize the finite element problem.
  DiscreteProblem dp(&wf, &space);

  // Project the function "titer" on the FE space 
  // in order to obtain initial vector for NOX. 
  info("Projecting initial solution on the FE mesh.");
  scalar* coeff_vec = new scalar[ndof];
  OGProjection::project_global(&space, &t_prev_time, coeff_vec);

  // Initialize NOX solver.
  NoxSolver solver(&dp);

  // Select preconditioner.
  RCP<Precond> pc = rcp(new MlPrecond("sa"));
  if (PRECOND)
  {
    if (JFNK) solver.set_precond(pc);
    else solver.set_precond("ML");
  }

  // Time stepping loop:
  double total_time = 0.0;
  for (int ts = 1; total_time <= 2000.0; ts++)
  {
    info("---- Time step %d, t = %g s", ts, total_time += TAU);

    info("Assembling by DiscreteProblem, solving by NOX.");
    solver.set_init_sln(coeff_vec);
    if (solver.solve())
      Solution::vector_to_solution(solver.get_solution(), &space, &t_prev_time);
    else
      error("NOX failed.");

    info("Number of nonlin iterations: %d (norm of residual: %g)", 
      solver.get_num_iters(), solver.get_residual());
    info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
      solver.get_num_lin_iters(), solver.get_achieved_tol());
  }

  info("Coordinate ( 0.6,  0.6) t_prev_time value = %lf", t_prev_time.get_pt_value( 0.6,  0.6));
  info("Coordinate ( 0.4,  0.6) t_prev_time value = %lf", t_prev_time.get_pt_value( 0.4,  0.6));
  info("Coordinate ( 0.4,  0.4) t_prev_time value = %lf", t_prev_time.get_pt_value( 0.4,  0.4));
  info("Coordinate ( 0.6,  0.0) t_prev_time value = %lf", t_prev_time.get_pt_value( 0.6,  0.0));
  info("Coordinate ( 0.5,  0.5) t_prev_time value = %lf", t_prev_time.get_pt_value( 0.5,  0.5));

#define ERROR_SUCCESS                                0
#define ERROR_FAILURE                               -1
  int success = 1;
  double eps = 1e-5;
  if (fabs(t_prev_time.get_pt_value(0.6, 0.6) - 13.826133) > eps) {
    printf("Coordinate (0.6, 0.6) t_prev_time value is %g\n", t_prev_time.get_pt_value(0.6, 0.6));
    success = 0;
  }
  if (fabs(t_prev_time.get_pt_value( 0.4, 0.6) - 13.826133) > eps) {
    printf("Coordinate ( 0.4, 0.6) t_prev_time value is %g\n", t_prev_time.get_pt_value( 0.4, 0.6));
    success = 0;
  }
  if (fabs(t_prev_time.get_pt_value( 0.4,  0.4) - 12.599711) > eps) {
    printf("Coordinate ( 0.4,  0.4) t_prev_time value is %g\n", t_prev_time.get_pt_value( 0.4,  0.4));
    success = 0;
  }
  if (fabs(t_prev_time.get_pt_value(0.6,  0.0) - 10.000000) > eps) {
    printf("Coordinate (0.6,  0.0) t_prev_time value is %g\n", t_prev_time.get_pt_value(0.6,  0.0));
    success = 0;
  }
  if (fabs(t_prev_time.get_pt_value( 0.5,  0.5) - 12.916592) > eps) {
    printf("Coordinate ( 0.5,  0.5) t_prev_time value is %g\n", t_prev_time.get_pt_value( 0.5,  0.5));
    success = 0;
  }

  if (success == 1) {
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}
