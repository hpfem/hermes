#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "runge_kutta.h"

using namespace RefinementSelectors;

// This test makes sure that example "butcher" works correctly.

const int INIT_GLOB_REF_NUM = 3;                   // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 0;                    // Number of initial refinements towards boundary.
const int P_INIT = 2;                              // Initial polynomial degree.
const double time_step = 0.001;                    // Time step.
const double T_FINAL = 2*time_step;                // Time interval length.
const double NEWTON_TOL = 1e-5;                    // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                   // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;   // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                   // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Time integration. Choose one of the following methods, or define your own Butcher's table:
// Explicit_RK_1, Implicit_RK_1, Explicit_RK_2, Implicit_Crank_Nicolson_2_2, Implicit_SDIRK_2_2, 
// Implicit_Lobatto_IIIA_2_2, Implicit_Lobatto_IIIB_2_2, Implicit_Lobatto_IIIC_2_2, Explicit_RK_3, Explicit_RK_4,
// Implicit_Lobatto_IIIA_3_4, Implicit_Lobatto_IIIB_3_4, Implicit_Lobatto_IIIC_3_4, Implicit_Radau_IIA_3_5. 

ButcherTableType butcher_table_type = Implicit_SDIRK_2_2;

// Thermal conductivity (temperature-dependent).
// Note: for any u, this function has to be positive.
template<typename Real>
Real lam(Real u) 
{ 
  return 1 + pow(u, 4);
}

// Derivative of the thermal conductivity with respect to 'u'.
template<typename Real>
Real dlam_du(Real u) 
{ 
  return 4*pow(u, 3);
}

// This function is used to define Dirichlet boundary conditions.
double dir_lift(double x, double y, double& dx, double& dy) {
  dx = (y+10)/10.;
  dy = (x+10)/10.;
  return (x+10)*(y+10)/100.;
}

// Initial condition. It will be projected on the FE mesh 
// to obtain initial coefficient vector for the Newton's method.
scalar init_cond(double x, double y, double& dx, double& dy)
{ return dir_lift(x, y, dx, dy);}

// Boundary markers.
const int BDY_DIRICHLET = 1;

// Essential (Dirichlet) boundary condition markers.
scalar essential_bc_values(double x, double y)
{
  double dx, dy;
  return dir_lift(x, y, dx, dy);
}

// Heat sources (can be a general function of 'x' and 'y').
template<typename Real>
Real heat_src(Real x, Real y) { return 1.0;}

// Weak forms.
#include "forms.cpp"

// Main function.
int main(int argc, char* argv[])
{
  // Check number of command-line parameters.
  if (argc < 2) 
    error("Not enough parameters: Provide a Butcher's table type.");

  int b_type = atoi(argv[1]);
  info ("%d", b_type);

  switch (b_type)
  {
    case 1: butcher_table_type = Explicit_RK_1; break;
    case 2: butcher_table_type = Implicit_RK_1; break;

    case 3: butcher_table_type = Explicit_RK_2; break;
    case 4: butcher_table_type = Implicit_Crank_Nicolson_2_2; break;
    case 5: butcher_table_type = Implicit_SDIRK_2_2; break;
    case 6: butcher_table_type = Implicit_Lobatto_IIIA_2_2; break;
    case 7: butcher_table_type = Implicit_Lobatto_IIIB_2_2; break;
    case 8: butcher_table_type = Implicit_Lobatto_IIIC_2_2; break;

    case 9: butcher_table_type = Explicit_RK_3; break;

    case 10: butcher_table_type = Explicit_RK_4; break;
    case 11: butcher_table_type = Implicit_Lobatto_IIIA_3_4; break;
    case 12: butcher_table_type = Implicit_Lobatto_IIIB_3_4; break;
    case 13: butcher_table_type = Implicit_Lobatto_IIIC_3_4; break;

    case 14: butcher_table_type = Implicit_Radau_IIA_3_5; break;

    default: error("Admissible command-line options are from 1 to 14.");
  }

  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table_type);

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(BDY_DIRICHLET, INIT_BDY_REF_NUM);

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_DIRICHLET);

  // Enter Dirichlet boundary values.
  BCValues bc_values;
  bc_values.add_function(BDY_DIRICHLET, essential_bc_values);   

  // Create an H1 space with default shapeset.
  H1Space* space = new H1Space(&mesh, &bc_types, &bc_values, P_INIT);

  int ndof = Space::get_num_dofs(space);
  info("ndof = %d.", ndof);

  // Previous time level solution (initialized by the initial condition).
  Solution* u_prev_time = new Solution(&mesh, init_cond);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(stac_jacobian));
  wf.add_vector_form(callback(stac_residual));

  // Project the initial condition on the FE space to obtain initial solution coefficient vector.
  info("Projecting initial condition to translate initial condition into a vector.");
  scalar* coeff_vec = new scalar[ndof];
  OGProjection::project_global(space, u_prev_time, coeff_vec, matrix_solver);

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem dp(&wf, space, is_linear);

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
    Solution::vector_to_solution(coeff_vec, space, u_prev_time);

    // Update time.
    current_time += time_step;

    // Increase counter of time steps.
    ts++;
  } 
  while (current_time < T_FINAL);

  info("Coordinate (-8.0, -8.0) value = %lf", u_prev_time->get_pt_value(-8.0, -8.0));
  info("Coordinate (-5.0, -5.0) value = %lf", u_prev_time->get_pt_value(-5.0, -5.0));
  info("Coordinate (-3.0, -3.0) value = %lf", u_prev_time->get_pt_value(-3.0, -3.0));
  info("Coordinate ( 0.0,  0.0) value = %lf", u_prev_time->get_pt_value(0.0,  0.0));
  info("Coordinate ( 3.0,  3.0) value = %lf", u_prev_time->get_pt_value(3.0,  3.0));
  info("Coordinate ( 5.0,  5.0) value = %lf", u_prev_time->get_pt_value(5.0,  5.0));
  info("Coordinate ( 8.0,  8.0) value = %lf", u_prev_time->get_pt_value(8.0,  8.0));

  double coor_x_y[7] = {-8.0, -5.0, -3.0, 0.0, 3.0, 5.0, 8.0};
  bool success = true;

  switch (b_type)
  {
    case 1: if (fabs(u_prev_time->get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042560) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251886) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492025) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002150) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693363) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255564) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.266989) > 1E-6) success = false;
            break;

    case 2: if (fabs(u_prev_time->get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251887) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002152) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693351) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255921) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.264214) > 1E-6) success = false;
            break;

    case 3: if (fabs(u_prev_time->get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251887) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693354) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255822) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265515) > 1E-6) success = false;
            break;

    case 4: if (fabs(u_prev_time->get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251886) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693356) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255784) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265491) > 1E-6) success = false;
            break;

    case 5: if (fabs(u_prev_time->get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251886) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693355) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255791) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265487) > 1E-6) success = false;
            break;

    case 6: if (fabs(u_prev_time->get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251886) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693356) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255784) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265491) > 1E-6) success = false;
            break;

    case 7: if (fabs(u_prev_time->get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251886) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693356) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255784) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265491) > 1E-6) success = false;
            break;

    case 8: if (fabs(u_prev_time->get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251887) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693354) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255819) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265435) > 1E-6) success = false;
            break;

    case 9: if (fabs(u_prev_time->get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251886) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693355) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255795) > 1E-6) success = false;
            if (fabs(u_prev_time->get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265458) > 1E-6) success = false;
            break;

    case 10: if (fabs(u_prev_time->get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251886) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693355) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255797) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265487) > 1E-6) success = false;
             break;

    case 11: if (fabs(u_prev_time->get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251886) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693355) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255798) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265482) > 1E-6) success = false;
             break;

    case 12: if (fabs(u_prev_time->get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251886) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693355) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255798) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265482) > 1E-6) success = false;
             break;

    case 13: if (fabs(u_prev_time->get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251886) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693355) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255798) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265482) > 1E-6) success = false;
             break;

    case 14: if (fabs(u_prev_time->get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251886) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693355) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255798) > 1E-6) success = false;
             if (fabs(u_prev_time->get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265482) > 1E-6) success = false;
             break;
    default: error("Admissible command-line options are from 1 to 14.");
  }

  // Cleanup.
  delete [] coeff_vec;
  delete space;
  delete u_prev_time;

  if (success) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}
