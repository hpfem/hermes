#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using Teuchos::RCP;
using Teuchos::rcp;

// This test makes sure that example 44-trilinos-coupled works correctly.

const int INIT_REF_NUM = 2;                 // Number of initial uniform mesh refinements.
const int P_INIT = 2;                       // Initial polynomial degree of all mesh elements.
const bool JFNK = true;                     // true = jacobian-free method,
                                            // false = Newton
const int PRECOND = 2;                      // Preconditioning by jacobian (1) (less GMRES iterations, more time to create precond)
                                            // or by approximation of jacobian (2) (less time for precond creation, more GMRES iters).
                                            // in case of jfnk,
                                            // default Ifpack proconditioner in case of Newton.
const double TAU = 0.05;                    // Time step.
const double T_FINAL = 2*TAU + 1e-4;        // Time interval length.
const bool TRILINOS_OUTPUT = true;          // Display more details about nonlinear and linear solvers.

// Problem parameters.
const double Le    = 1.0;
const double alpha = 0.8;
const double beta  = 10.0;
const double kappa = 0.1;
const double x1    = 9.0;

// Initial conditions.
scalar temp_ic(double x, double y, scalar& dx, scalar& dy)
  { return (x <= x1) ? 1.0 : exp(x1 - x); }

scalar conc_ic(double x, double y, scalar& dx, scalar& dy)
  { return (x <= x1) ? 0.0 : 1.0 - exp(Le*(x1 - x)); }

// Weak forms. 
# include "../definitions.cpp"

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../domain.mesh", &mesh);

  // Perform initial mesh refinemets.
  for (int i=0; i < INIT_REF_NUM; i++)  mesh.refine_all_elements();

  // Initialize boundary conditions.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet("Left");
  bc_types.add_bc_neumann("Neumann");
  bc_types.add_bc_newton("Cooled");

  // Enter Dirichlet boundary values.
  BCValues bc_values_t;
  bc_values_t.add_const("Left", 1.0);

  BCValues bc_values_c;
  bc_values_c.add_zero("Left");

  // Create H1 spaces with default shapesets.
  H1Space* t_space = new H1Space(&mesh, &bc_types, &bc_values_t, P_INIT);
  H1Space* c_space = new H1Space(&mesh, &bc_types, &bc_values_c, P_INIT);
  int ndof = Space::get_num_dofs(Hermes::vector<Space *>(t_space, c_space));
  info("ndof = %d.", ndof);

  // Define initial conditions.
  Solution t_prev_time_1, c_prev_time_1, t_prev_time_2, 
           c_prev_time_2, t_iter, c_iter, t_prev_newton, c_prev_newton;
  t_prev_time_1.set_exact(&mesh, temp_ic);  
  c_prev_time_1.set_exact(&mesh, conc_ic);
  t_prev_time_2.set_exact(&mesh, temp_ic);  
  c_prev_time_2.set_exact(&mesh, conc_ic);
  t_iter.set_exact(&mesh, temp_ic);   
  c_iter.set_exact(&mesh, conc_ic);

  // Filters for the reaction rate omega and its derivatives.
  DXDYFilter omega(omega_fn, Hermes::vector<MeshFunction*>(&t_prev_time_1, &c_prev_time_1));
  DXDYFilter omega_dt(omega_dt_fn, Hermes::vector<MeshFunction*>(&t_prev_time_1, &c_prev_time_1));
  DXDYFilter omega_dc(omega_dc_fn, Hermes::vector<MeshFunction*>(&t_prev_time_1, &c_prev_time_1));

  // Initialize weak formulation.
  WeakForm wf(2, JFNK ? true : false);
  if (!JFNK || (JFNK && PRECOND == 1))
  {
    wf.add_matrix_form(0, 0, callback(newton_bilinear_form_0_0), HERMES_NONSYM, HERMES_ANY, &omega_dt);
    wf.add_matrix_form_surf(0, 0, callback(newton_bilinear_form_0_0_surf), 3);
    wf.add_matrix_form(1, 1, callback(newton_bilinear_form_1_1), HERMES_NONSYM, HERMES_ANY, &omega_dc);
    wf.add_matrix_form(0, 1, callback(newton_bilinear_form_0_1), HERMES_NONSYM, HERMES_ANY, &omega_dc);
    wf.add_matrix_form(1, 0, callback(newton_bilinear_form_1_0), HERMES_NONSYM, HERMES_ANY, &omega_dt);
  }
  else if (PRECOND == 2)
  {
    wf.add_matrix_form(0, 0, callback(precond_0_0));
    wf.add_matrix_form(1, 1, callback(precond_1_1));
  }
  wf.add_vector_form(0, callback(newton_linear_form_0), HERMES_ANY, 
                     Hermes::vector<MeshFunction*>(&t_prev_time_1, &t_prev_time_2, &omega));
  wf.add_vector_form_surf(0, callback(newton_linear_form_0_surf), 3);
  wf.add_vector_form(1, callback(newton_linear_form_1), HERMES_ANY, 
                     Hermes::vector<MeshFunction*>(&c_prev_time_1, &c_prev_time_2, &omega));

  // Project the functions "t_iter" and "c_iter" on the FE space 
  // in order to obtain initial vector for NOX. 
  info("Projecting initial solutions on the FE meshes.");
  scalar* coeff_vec = new scalar[ndof];
  OGProjection::project_global(Hermes::vector<Space *>(t_space, c_space), Hermes::vector<MeshFunction*>(&t_prev_time_1, &c_prev_time_1),
    coeff_vec);

  // Measure the projection time.
  double proj_time = cpu_time.tick().last();

  // Initialize finite element problem.
  DiscreteProblem dp(&wf, Hermes::vector<Space*>(t_space, c_space));

  // Initialize NOX solver and preconditioner.
  NoxSolver solver(&dp);
  RCP<Precond> pc = rcp(new MlPrecond("sa"));
  if (PRECOND)
  {
    if (JFNK) solver.set_precond(pc);
    else solver.set_precond("Ifpack");
  }
  if (TRILINOS_OUTPUT)
    solver.set_output_flags(NOX::Utils::Error | NOX::Utils::OuterIteration |
                            NOX::Utils::OuterIterationStatusTest |
                            NOX::Utils::LinearSolverDetails);

  // Time stepping loop:
  double total_time = 0.0;
  cpu_time.tick_reset();
  for (int ts = 1; total_time <= T_FINAL; ts++)
  {
    info("---- Time step %d, t = %g s", ts, total_time + TAU);

    cpu_time.tick(HERMES_SKIP);
    solver.set_init_sln(coeff_vec);
    if (solver.solve())
    {
      Solution::vector_to_solutions(solver.get_solution(), Hermes::vector<Space *>(t_space, c_space), Hermes::vector<Solution *>(&t_prev_newton, &c_prev_newton));

      cpu_time.tick();
      info("Number of nonlin iterations: %d (norm of residual: %g)",
          solver.get_num_iters(), solver.get_residual());
      info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)",
          solver.get_num_lin_iters(), solver.get_achieved_tol());

      // Time measurement.
      cpu_time.tick(HERMES_SKIP);

      // Update global time.
      total_time += TAU;

      // Saving solutions for the next time step.
      t_prev_time_2.copy(&t_prev_time_1);
      c_prev_time_2.copy(&c_prev_time_1);
      t_prev_time_1 = t_prev_newton;
      c_prev_time_1 = c_prev_newton;
    }
    else
      error("NOX failed.");

    info("Total running time for time level %d: %g s.", ts, cpu_time.tick().last());
  }

  info("T Coordinate (  0,   8) value = %lf", t_prev_time_1.get_pt_value(0.0, 8.0));
  info("T Coordinate (  8,   8) value = %lf", t_prev_time_1.get_pt_value(8.0, 8.0));
  info("T Coordinate ( 15,   8) value = %lf", t_prev_time_1.get_pt_value(15.0, 8.0));
  info("T Coordinate ( 24,   8) value = %lf", t_prev_time_1.get_pt_value(24.0, 8.0));
  info("T Coordinate ( 30,   8) value = %lf", t_prev_time_1.get_pt_value(30.0, 8.0));
  info("T Coordinate ( 40,   8) value = %lf", t_prev_time_1.get_pt_value(40.0, 8.0));
  info("T Coordinate ( 50,   8) value = %lf", t_prev_time_1.get_pt_value(50.0, 8.0));
  info("T Coordinate ( 60,   8) value = %lf", t_prev_time_1.get_pt_value(60.0, 8.0));

  info("C Coordinate (  0,   8) value = %lf", c_prev_time_1.get_pt_value(0.0, 8.0));
  info("C Coordinate (  8,   8) value = %lf", c_prev_time_1.get_pt_value(8.0, 8.0));
  info("C Coordinate ( 15,   8) value = %lf", c_prev_time_1.get_pt_value(15.0, 8.0));
  info("C Coordinate ( 24,   8) value = %lf", c_prev_time_1.get_pt_value(24.0, 8.0));
  info("C Coordinate ( 30,   8) value = %lf", c_prev_time_1.get_pt_value(30.0, 8.0));
  info("C Coordinate ( 40,   8) value = %lf", c_prev_time_1.get_pt_value(40.0, 8.0));
  info("C Coordinate ( 50,   8) value = %lf", c_prev_time_1.get_pt_value(50.0, 8.0));
  info("C Coordinate ( 60,   8) value = %lf", c_prev_time_1.get_pt_value(60.0, 8.0));

  double coor_x[8] = {0.0, 8.0, 15.0, 24.0, 30.0, 40.0, 50.0, 60.0};
  double coor_y = 8.0;
  double t_value[8] = {1.000000, 1.000078, 0.002819, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  double c_value[8] = {0.000000, -0.000078, 0.997181, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000};

  bool success = true;
  for (int i = 0; i < 8; i++)
  {
    if ((abs(t_value[i] - t_prev_time_1.get_pt_value(coor_x[i], coor_y)) > 1E-4) || 
        (abs(c_value[i] - c_prev_time_1.get_pt_value(coor_x[i], coor_y)) > 1E-4))
      success = false;
  }

  if (success) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}

