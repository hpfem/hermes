#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "runge_kutta.h"

// This test makes sure that example 09-timedep works correctly.


const int P_INIT = 2;                             // Polynomial degree of all mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 3;                   // Number of initial uniform mesh refinements towards the boundary.
const double time_step = 3e+2;                    // Time step in seconds.
const double NEWTON_TOL = 1e-5;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Time integration. Choose one of the following methods, or define your own Butcher's table:
// Explicit_RK_1, Implicit_RK_1, Explicit_RK_2, Implicit_Crank_Nicolson_2_2, Implicit_SIRK_2_2, 
// Implicit_ESIRK_2_2, Implicit_SDIRK_2_2, Implicit_Lobatto_IIIA_2_2, Implicit_Lobatto_IIIB_2_2, 
// Implicit_Lobatto_IIIC_2_2, Explicit_RK_3, Explicit_RK_4, Implicit_Lobatto_IIIA_3_4, 
// Implicit_Lobatto_IIIB_3_4, Implicit_Lobatto_IIIC_3_4, Implicit_Radau_IIA_3_5, Implicit_SDIRK_4_5. 

ButcherTableType butcher_table_type = Implicit_SDIRK_2_2;

// Boundary markers.
const std::string BDY_GROUND = "Boundary ground";
const std::string BDY_AIR = "Boundary air";

// Problem parameters.
const double TEMP_INIT = 10;       // Temperature of the ground (also initial temperature).
const double ALPHA = 10;           // Heat flux coefficient for Newton's boundary condition.
const double LAMBDA = 1e5;         // Thermal conductivity of the material.
const double HEATCAP = 1e6;        // Heat capacity.
const double RHO = 3000;           // Material density.
const double T_FINAL = 86400;      // Length of time interval (24 hours) in seconds.

// Time-dependent exterior temperature.
template<typename Real>
Real temp_ext(Real t) {
  return TEMP_INIT + 10. * sin(2*M_PI*t/T_FINAL);
}

// Heat sources (can be a general function of 'x' and 'y').
template<typename Real>
Real heat_src(Real x, Real y) { return 0.0;}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table_type);

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("cathedral.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(BDY_AIR, INIT_REF_NUM_BDY);
  mesh.refine_towards_boundary(BDY_GROUND, INIT_REF_NUM_BDY);

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(Hermes::vector<std::string>(BDY_GROUND));
  bc_types.add_bc_newton(BDY_AIR);

  // Enter Dirichlet boundary values.
  BCValues bc_values;
  bc_values.add_const(BDY_GROUND, TEMP_INIT);

  // Initialize an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d.", ndof);
 
  // Previous time level solution (initialized by the external temperature).
  Solution u_prev_time(&mesh, TEMP_INIT);

  // Initialize weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(stac_jacobian));
  wf.add_vector_form(callback(stac_residual));
  wf.add_matrix_form_surf(callback(bilinear_form_surf), BDY_AIR);
  wf.add_vector_form_surf(callback(linear_form_surf), BDY_AIR);

  // Project the initial condition on the FE space to obtain initial solution coefficient vector.
  info("Projecting initial condition to translate initial condition into a vector.");
  scalar* coeff_vec = new scalar[ndof];
  OGProjection::project_global(&space, &u_prev_time, coeff_vec, matrix_solver);

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem dp(&wf, &space, is_linear);

  // Time stepping loop:
  double current_time = 0.0; int ts = 1;
  do 
  {
    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    info("Runge-Kutta time step (t = %g, tau = %g, stages: %d).", 
         current_time, time_step, bt.get_size());
    bool verbose = true;
    if (!rk_time_step(current_time, time_step, &bt, coeff_vec, &dp, matrix_solver,
		      verbose, NEWTON_TOL, NEWTON_MAX_ITER)) {
      error("Runge-Kutta time step failed, try to decrease time step size.");
    }

    // Convert coeff_vec into a new time level solution.
    Solution::vector_to_solution(coeff_vec, &space, &u_prev_time);

    // Update time.
    current_time += time_step;

    // Increase counter of time steps.
    ts++;
  } 
  while (current_time < time_step*5);

  double sum = 0;
  for (int i = 0; i < ndof; i++) 
    sum += coeff_vec[i];

  printf("sum = %g\n", sum);

  int success = 1;
  if (fabs(sum - 3617.55) > 1e-1) success = 0;

  // Cleanup.
  delete [] coeff_vec;

  if (success == 1) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}

