#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
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

const int BDY_BOTTOM = 1;

// Boundary condition types.
BCType bc_types(int marker)
{
  if (marker == BDY_BOTTOM) return BC_ESSENTIAL;
  else return BC_NATURAL;
}

// Essential (Dirichlet) boundary conditions values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return TEMP_INIT;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Perform initial mesh refinemets.
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);
  int ndof = get_num_dofs(&space);
  info("ndof: %d", ndof);

  // Define constant initial condition. 
  Solution t_prev_time;
  t_prev_time.set_const(&mesh, 20.0);

  // Initialize the weak formulation.
  WeakForm wf(1, JFNK ? true : false);
  wf.add_matrix_form(callback(jacobian));
  wf.add_matrix_form_surf(callback(jacobian_surf));
  wf.add_vector_form(callback(residual), H2D_ANY, &t_prev_time);
  wf.add_vector_form_surf(callback(residual_surf));

  // Initialize the finite element problem.
  FeProblem fep(&wf, &space);

  // Project the function "titer" on the FE space 
  // in order to obtain initial vector for NOX. 
  info("Projecting initial solution on the FE mesh.");
  Vector* coeff_vec = new AVector(ndof);
  project_global(&space, H2D_H1_NORM, &t_prev_time, &t_prev_time, coeff_vec);

  // Measure the projection time.
  double proj_time = cpu_time.tick().last();

  // Initialize NOX solver.
  NoxSolver solver(&fep);

  // Select preconditioner.
  RCP<Precond> pc = rcp(new MlPrecond("sa"));
  if (PRECOND)
  {
    if (JFNK) solver.set_precond(pc);
    else solver.set_precond("ML");
  }

  // Time stepping loop:
  double total_time = 0.0;
  cpu_time.tick_reset();
  for (int ts = 1; total_time <= 2000.0; ts++)
  {
    info("---- Time step %d, t = %g s", ts, total_time += TAU);

    info("Assembling by FeProblem, solving by NOX.");
    solver.set_init_sln(coeff_vec->get_c_array());
    bool solved = solver.solve();
    if (solved)
    {
      double *s = solver.get_solution_vector();
      AVector *tmp_vector = new AVector(ndof);
      tmp_vector->set_c_array(s, ndof);
      t_prev_time.set_fe_solution(&space, tmp_vector);
      delete tmp_vector;
    }
    else
      error("NOX failed.");

    info("Number of nonlin iterations: %d (norm of residual: %g)", 
      solver.get_num_iters(), solver.get_residual());
    info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
      solver.get_num_lin_iters(), solver.get_achieved_tol());
  }

  info("Total running time: %g", cpu_time.accumulated());

  info("Coordinate ( 0.6,  0.6) t_prev_time value = %lf", t_prev_time.get_pt_value( 0.6,  0.6));
  info("Coordinate ( 0.4,  0.6) t_prev_time value = %lf", t_prev_time.get_pt_value( 0.4,  0.6));
  info("Coordinate ( 0.4,  0.4) t_prev_time value = %lf", t_prev_time.get_pt_value( 0.4,  0.4));
  info("Coordinate ( 0.6,  0.0) t_prev_time value = %lf", t_prev_time.get_pt_value( 0.6,  0.0));
  info("Coordinate ( 0.5,  0.5) t_prev_time value = %lf", t_prev_time.get_pt_value( 0.5,  0.5));

#define ERROR_SUCCESS                                0
#define ERROR_FAILURE                               -1
  int success = 1;
  double eps = 1e-5;
  if (fabs(t_prev_time.get_pt_value(0.6, 0.6) - 15.241357) > eps) {
    printf("Coordinate (0.6, 0.6) t_prev_time value is %g\n", t_prev_time.get_pt_value(0.6, 0.6));
    success = 0;
  }
  if (fabs(t_prev_time.get_pt_value( 0.4, 0.6) - 15.241357) > eps) {
    printf("Coordinate ( 0.4, 0.6) t_prev_time value is %g\n", t_prev_time.get_pt_value( 0.4, 0.6));
    success = 0;
  }
  if (fabs(t_prev_time.get_pt_value( 0.4,  0.4) - 13.825044) > eps) {
    printf("Coordinate ( 0.4,  0.4) t_prev_time value is %g\n", t_prev_time.get_pt_value( 0.4,  0.4));
    success = 0;
  }
  if (fabs(t_prev_time.get_pt_value(0.6,  0.0) - 10.000000) > eps) {
    printf("Coordinate (0.6,  0.0) t_prev_time value is %g\n", t_prev_time.get_pt_value(0.6,  0.0));
    success = 0;
  }
  if (fabs(t_prev_time.get_pt_value( 0.5,  0.5) - 14.345194) > eps) {
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
