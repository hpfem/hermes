#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "function/function.h"

using namespace RefinementSelectors;

//  This example is derived from example "butcher" that handles 
//  time integration with general Butcher's tables and thus has some 
//  overhead. This example only implements the implicit Euler method 
//  and an SDIRK-2 method in a straightforward fashion, and thus is 
//  faster than example "butcher" for the same methods. The sole 
//  purpose of this example is a performance comparison with example 
//  "butcher". The Butcher's table of the SDIRK-2 method is:

double GAMMA = 1. - 1./sqrt(2.);
double BUTCHER_A_11 = 1. - 1./sqrt(2.);
double BUTCHER_A_12 = 0.;
double BUTCHER_A_21 = 1. - (1. - 1./sqrt(2.));
double BUTCHER_A_22 = 1. - 1./sqrt(2.);
double BUTCHER_B_1 = 1. - (1. - 1./sqrt(2.));
double BUTCHER_B_2 = 1. - 1./sqrt(2.);
double BUTCHER_C_1 = 1. - 1./sqrt(2.);
double BUTCHER_C_2 = 1.;

//  The method can be found in Butcher's book on page 244. 
//
//  PDE: time-dependent heat transfer equation with nonlinear thermal
//  conductivity, du/dt - div[lambda(u)grad u] = f.
//
//  Domain: square (-10, 10)^2.
//
//  BC: Dirichlet, given by the function dir_lift() below.
//  IC: Same function dir_lift().
//
//  The following parameters can be changed:

const int TIME_INTEGRATION = 2;                   // 1 = implicit Euler, 2 = SDIRK-2
const int INIT_GLOB_REF_NUM = 3;                  // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 4;                   // Number of initial refinements towards boundary.
const int P_INIT = 2;                             // Initial polynomial degree.
double TAU = 0.2;                                 // Time step.
const double T_FINAL = 5.0;                       // Time interval length.
const double NEWTON_TOL = 1e-5;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

const double ALPHA = 4.0;                         // For the nonlinear thermal ocnductivity.
double TIME = 0.0;

// Thermal conductivity (temperature-dependent).
// Note: for any u, this function has to be positive.
template<typename Real>
Real lam(Real u) { return 1 + pow(u, ALPHA);}

// Derivative of the thermal conductivity with respect to 'u'.
template<typename Real>
Real dlam_du(Real u) { return ALPHA*pow(u, ALPHA-1);}

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

// Heat sources (forcing term in accordance with exact solution).
template<typename Real>
Real heat_src(Real x, Real y, double t) { return 1.0;}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[]) 
{
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
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d.", ndof);

  // Previous time level solution (initialized by the initial condition).
  Solution u_prev_time(&mesh, init_cond);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  info("Projecting initial condition to obtain initial vector for the Newton's method.");
  scalar* coeff_vec1 = new scalar[ndof];
  OGProjection::project_global(&space, &u_prev_time, coeff_vec1, matrix_solver);

  // Initialize views.
  ScalarView sview("Solution", new WinGeom(0, 0, 500, 400));
  OrderView oview("Mesh", new WinGeom(520, 0, 450, 400));
  oview.show(&space);

  // Initialize the weak formulation and the FE problem.
  WeakForm wf, wf1, wf2;
  Solution sdirk_stage_sol;
  scalar* coeff_vec2;
  DiscreteProblem *dp, *dp1, *dp2;
  if (TIME_INTEGRATION == 1) {
    wf.add_matrix_form(callback(jac_implicit_euler));
    wf.add_vector_form(callback(res_implicit_euler), HERMES_ANY, &u_prev_time);
 
    bool is_linear = false; 
    dp = new DiscreteProblem(&wf, &space, is_linear); 
  }
  else {
    sdirk_stage_sol.set_exact(&mesh, init_cond);

    coeff_vec2 = new scalar[ndof];
    for (int i = 0; i < ndof; i++) coeff_vec2[i] = coeff_vec1[i];

    wf1.add_matrix_form(callback(jac_sdirk));
    wf1.add_vector_form(callback(res_sdirk_stage_1), HERMES_ANY, 
                        Hermes::vector<MeshFunction*>(&u_prev_time));
    wf2.add_matrix_form(callback(jac_sdirk));
    wf2.add_vector_form(callback(res_sdirk_stage_2), HERMES_ANY, 
                        Hermes::vector<MeshFunction*>(&u_prev_time, &sdirk_stage_sol));

    bool is_linear = false;
    dp1 = new DiscreteProblem(&wf1, &space, is_linear);
    dp2 = new DiscreteProblem(&wf2, &space, is_linear);
  }

  // Time stepping loop:
  int ts = 1;
  do {
    if (TIME_INTEGRATION == 1) {
      // Perform Newton's iteration.
      info("Implicit Euler time step (t = %g, tau = %g):", TIME, TAU);
      bool verbose = true;
      if (!solve_newton(coeff_vec1, dp, solver, matrix, rhs, 
          NEWTON_TOL, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");  

      // Update previous time level solution.
      Solution::vector_to_solution(coeff_vec1, &space, &u_prev_time);
    } 
    else {
      // Perform Newton's iteration for sdirk_stage_sol.
      info("SDIRK-2 time step, stage I (t = %g, tau = %g):", TIME, TAU);
      bool verbose = true;
      if (!solve_newton(coeff_vec1, dp1, solver, matrix,
			rhs, NEWTON_TOL, NEWTON_MAX_ITER, verbose))
        error("Newton's iteration did not converge."); 

      // Convert the vector coeff_vec1 into a Solution.
      Solution::vector_to_solution(coeff_vec1, &space, &sdirk_stage_sol);

      // Perform Newton's iteration for the final solution.
      info("SDIRK-2 time step, stage II (t = %g, tau = %g):", TIME, TAU);

      if (!solve_newton(coeff_vec2, dp2, solver, matrix,
			rhs, NEWTON_TOL, NEWTON_MAX_ITER, verbose))
        error("Newton's iteration did not converge."); 

      // Translate Y2 into a Solution.
      Solution::vector_to_solution(coeff_vec2, &space, &u_prev_time);
    }
  
    // Update time.
    TIME = TIME + TAU;

    // Show the new time level solution.
    char title[100];
    sprintf(title, "Solution, t = %g", TIME);
    sview.set_title(title);
    sview.show(&u_prev_time);

    ts++;
  } while (TIME < T_FINAL);

  // Clean up.
  delete [] coeff_vec1;
  if (TIME_INTEGRATION == 1) {
    delete dp;
  }
  else {
    delete dp1;
    delete dp2;
    delete [] coeff_vec2;
  }
  delete matrix;
  delete solver;

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
