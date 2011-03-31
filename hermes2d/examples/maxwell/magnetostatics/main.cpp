#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "function/function.h"

using namespace RefinementSelectors;

//  This example shows how to set an arbitrary initial guess for the
//  Newton's method, and nonzero Dirichlet boundary conditions.
//
//  PDE: stationary heat transfer equation with nonlinear thermal
//  conductivity, - div[lambda(u)grad u] = heat_src.
//
//  Domain: unit square (-10,10)^2.
//
//  BC: Dirichlet, see function dir_lift() below.
//
//  The following parameters can be changed:

const int P_INIT = 2;                             // Initial polynomial degree.
const double NEWTON_TOL = 1e-6;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.
const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
double MU_VACUUM = 4 * M_PI * 1e-7;
double INIT_COND = 0.0;                           // Initial condition for the magnetic potential.
double CURRENT_DENSITY = 1e6;                     // Volume source term.

// Material and boundary markers.
const std::string MAT_AIR = "0";
const std::string MAT_IRON = "1";
const std::string MAT_COPPER = "3";
const std::string BDY_DIRICHLET = "1";

// Weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Define nonlinear thermal conductivity lambda(u) via a cubic spline.
  Hermes::vector<double> lambda_pts(0.0,    0.5,    0.7,    0.8,    0.9,    1.0,    1.1,    1.2,    1.3,   1.4,   1.6,   1.7,   1.8,   1.9,   3.0);
  Hermes::vector<double> lambda_val(1500.0, 1480.0, 1460.0, 1450.0, 1440.0, 1400.0, 1300.0, 1150.0, 950.0, 750.0, 250.0, 180.0, 175.0, 150.0, 20.0);
  for (unsigned int i=0; i < lambda_val.size(); i++) lambda_val[i] *= MU_VACUUM; 

  // Create the cubic spline (and plot it for visual control). 
  double second_der_left = 0.0;
  double second_der_right = 0.0;
  bool first_der_left = false;
  bool first_der_right = false;
  bool extrapolate_der_left = false;
  bool extrapolate_der_right = false;
  CubicSpline mu_iron(lambda_pts, lambda_val, 0.0, 0.0, first_der_left, first_der_right,
                      extrapolate_der_left, extrapolate_der_right);
  bool success = mu_iron.calculate_coeffs(); 
  if (!success) error("There was a problem constructing a cubic spline.");
  info("Saving cubic spline into a Pylab file spline.dat.");
  double interval_extension = 0.0; // The interval of definition of the spline will be 
                                   // extended by "interval_extension" on both sides.
  mu_iron.plot("spline.dat", interval_extension);

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("actuator.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  MeshView mv("Mesh", new WinGeom(0, 0, 400, 400));
  mv.show(&mesh);

  // Initialize boundary conditions.
  DefaultEssentialBCConst bc_essential(BDY_DIRICHLET, 0.0);
  EssentialBCs bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  info("ndof: %d", Space::get_num_dofs(&space));

  // Initialize the weak formulation
  CustomWeakFormMagnetostatics wf(&mu_iron, MAT_AIR, MU_VACUUM, MAT_COPPER, 
                                  MU_VACUUM, CURRENT_DENSITY);

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem dp(&wf, &space, is_linear);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Initialize the solution.
  Solution sln(&mesh, INIT_COND);

  // Project the initial condition on the FE space to obtain initial 
  // coefficient vector for the Newton's method.
  info("Projecting to obtain initial vector for the Newton's method.");
  scalar* coeff_vec = new scalar[Space::get_num_dofs(&space)] ;
  OGProjection::project_global(&space, &sln, coeff_vec, matrix_solver); 

  // Perform Newton's iteration.
  bool verbose = true;
  if (!hermes2d.solve_newton(coeff_vec, &dp, solver, matrix, rhs, 
      NEWTON_TOL, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

  // Translate the resulting coefficient vector into the Solution sln.
  Solution::vector_to_solution(coeff_vec, &space, &sln);

  // Cleanup.
  delete []coeff_vec;
  delete matrix;
  delete rhs;
  delete solver;

  // Visualise the solution and mesh.
  ScalarView s_view("Solution", new WinGeom(0, 0, 440, 350));
  s_view.show_mesh(false);
  s_view.show(&sln);
  OrderView o_view("Mesh", new WinGeom(450, 0, 400, 350));
  o_view.show(&space);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

