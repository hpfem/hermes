#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "runge_kutta.h"

// This test makes sure that example 09-timedep-rk works correctly.


const int P_INIT = 2;                             // Polynomial degree of all mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 3;                   // Number of initial uniform mesh refinements towards the boundary.
const double time_step = 3e+2;                    // Time step in seconds.
const double NEWTON_TOL = 1e-5;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Choose one of the following time-integration methods, or define your own Butcher's table. The last number 
// in the name of each method is its order. The one before last, if present, is the number of stages.
// Explicit methods:
//   Explicit_RK_1, Explicit_RK_2, Explicit_RK_3, Explicit_RK_4.
// Implicit methods: 
//   Implicit_RK_1, Implicit_Crank_Nicolson_2_2, Implicit_SIRK_2_2, Implicit_ESIRK_2_2, Implicit_SDIRK_2_2, 
//   Implicit_Lobatto_IIIA_2_2, Implicit_Lobatto_IIIB_2_2, Implicit_Lobatto_IIIC_2_2, Implicit_Lobatto_IIIA_3_4, 
//   Implicit_Lobatto_IIIB_3_4, Implicit_Lobatto_IIIC_3_4, Implicit_Radau_IIA_3_5, Implicit_SDIRK_4_5.
// Embedded explicit methods:
//   Explicit_HEUN_EULER_2_12_embedded, Explicit_BOGACKI_SHAMPINE_4_23_embedded, Explicit_FEHLBERG_6_45_embedded,
//   Explicit_CASH_KARP_6_45_embedded, Explicit_DORMAND_PRINCE_7_45_embedded.
// Embedded implicit methods:
//   Implicit_SDIRK_CASH_3_23_embedded, Implicit_ESDIRK_TRBDF2_3_23_embedded, Implicit_ESDIRK_TRX2_3_23_embedded, 
//   Implicit_SDIRK_CASH_5_24_embedded, Implicit_SDIRK_CASH_5_34_embedded, Implicit_DIRK_7_45_embedded. 
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
const double T_FINAL = time_step*5;// Length of time interval (24 hours) in seconds.

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
  if (bt.is_explicit()) info("Using a %d-stage explicit R-K method.", bt.get_size());
  if (bt.is_diagonally_implicit()) info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
  if (bt.is_fully_implicit()) info("Using a %d-stage fully implicit R-K method.", bt.get_size());

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
 
  // Previous and next time level solutions.
  Solution* sln_time_prev = new Solution(&mesh, TEMP_INIT);
  Solution* sln_time_new = new Solution(&mesh);

  // Initialize weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(stac_jacobian_vol));
  wf.add_vector_form(callback(stac_residual_vol), HERMES_ANY, sln_time_prev);
  wf.add_matrix_form_surf(callback(stac_jacobian_surf), BDY_AIR, sln_time_prev);
  wf.add_vector_form_surf(callback(stac_residual_surf), BDY_AIR, sln_time_prev);

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem dp(&wf, &space, is_linear);

  // Time stepping loop:
  double current_time = time_step; int ts = 1;
  do 
  {
    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    info("Runge-Kutta time step (t = %g s, tau = %g s, stages: %d).", 
         current_time, time_step, bt.get_size());
    bool verbose = true;
    bool is_linear = true;
    if (!rk_time_step(current_time, time_step, &bt, sln_time_prev, sln_time_new, &dp, matrix_solver,
		      verbose, is_linear)) {
      error("Runge-Kutta time step failed, try to decrease time step size.");
    }

    // Convert coeff_vec into a new time level solution.
    //Solution::vector_to_solution(coeff_vec, &space, &u_prev_time);

    // Increase current time and time step counter.
    current_time += time_step;
    ts++;
  } 
  while (current_time < T_FINAL);

  info("Coordinate (-2.0, 2.0) value = %lf", sln_time_new->get_pt_value(-2.0, 2.0));
  info("Coordinate (-1.0, 2.0) value = %lf", sln_time_new->get_pt_value(-1.0, 2.0));
  info("Coordinate ( 0.0, 2.0) value = %lf", sln_time_new->get_pt_value(0.0, 2.0));
  info("Coordinate ( 1.0, 2.0) value = %lf", sln_time_new->get_pt_value(1.0, 2.0));
  info("Coordinate ( 2.0, 2.0) value = %lf", sln_time_new->get_pt_value(2.0, 2.0));

  bool success = true;

  if (fabs(sln_time_new->get_pt_value(-2.0, 2.0) - 10.000033) > 1E-6) success = false;
  if (fabs(sln_time_new->get_pt_value(-1.0, 2.0) - 9.999999) > 1E-6) success = false;
  if (fabs(sln_time_new->get_pt_value(0.0, 2.0) - 10.000005) > 1E-6) success = false;
  if (fabs(sln_time_new->get_pt_value(1.0, 2.0) - 9.999999) > 1E-6) success = false;
  if (fabs(sln_time_new->get_pt_value(2.0, 2.0) - 10.000033) > 1E-6) success = false;


  if (success) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}

