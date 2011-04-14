#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "function/function.h"

using namespace RefinementSelectors;

//  This example shows how to set an arbitrary initial guess for the
//  Newton's method, and nonzero Dirichlet boundary conditions.
//
//  PDE: magnetostatics with nonlinear magnetic permeability
//  curl[1/mu curl u] = current_density.
//
//  Domain: unit square (-10,10)^2.
//
//  BC: Dirichlet, see function dir_lift() below.
//
//  The following parameters can be changed:

const int P_INIT = 3;                             // Initial polynomial degree.
const double NEWTON_TOL = 1e-10;                  // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 1000;                 // Maximum allowed number of Newton iterations.
const double NEWTON_DAMPING = 1.0;                // Number between 0 and 1 to damp Newton's iterations.
const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
double MU_VACUUM = 4. * M_PI * 1e-7;
double INIT_COND = 0.0;                           // Initial condition for the magnetic potential.
double CURRENT_DENSITY = 1e9;                     // Volume source term.

// Material and boundary markers.
const std::string MAT_AIR = "2";
const std::string MAT_IRON_1 = "0";
const std::string MAT_IRON_2 = "3";
const std::string MAT_COPPER = "1";
const std::string BDY_DIRICHLET = "1";

// Weak forms.
#include "../definitions.cpp"

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Define nonlinear magnetic permeability via a cubic spline.
  Hermes::vector<double> mu_inv_pts(0.0,      0.5,      0.9,      1.0,      1.1,      1.2,      1.3,   
                                    1.4,      1.6,      1.7,      1.8,      1.9,      3.0,      5.0,     10.0);
  Hermes::vector<double> mu_inv_val(1/1500.0, 1/1480.0, 1/1440.0, 1/1400.0, 1/1300.0, 1/1150.0, 1/950.0,  
                                    1/750.0,  1/250.0,  1/180.0,  1/175.0,  1/150.0,  1/20.0,   1/10.0,  1/5.0);

  /* This is for debugging (iron is assumed linear with mu_r = 300.0
  Hermes::vector<double> mu_inv_pts(0.0,      10.0);
  Hermes::vector<double> mu_inv_val(1/300.0,   1/300.0);
  */

  // Create the cubic spline (and plot it for visual control). 
  double second_der_left = 0.0;
  double second_der_right = 0.0;
  bool first_der_left = false;
  bool first_der_right = false;
  bool extrapolate_der_left = false;
  bool extrapolate_der_right = false;
  CubicSpline mu_inv_iron(mu_inv_pts, mu_inv_val, 0.0, 0.0, first_der_left, first_der_right,
                          extrapolate_der_left, extrapolate_der_right);
  bool success0 = mu_inv_iron.calculate_coeffs(); 
  if (!success0) error("There was a problem constructing a cubic spline.");
  info("Saving cubic spline into a Pylab file spline.dat.");
  double interval_extension = 1.0; // The interval of definition of the spline will be 
                                   // extended by "interval_extension" on both sides.
  mu_inv_iron.plot("spline.dat", interval_extension, true);
  mu_inv_iron.plot("spline_der.dat", interval_extension, false);

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../actuator.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  //MeshView mv("Mesh", new WinGeom(0, 0, 400, 400));
  //mv.show(&mesh);

  // Initialize boundary conditions.
  DefaultEssentialBCConst bc_essential(BDY_DIRICHLET, 0.0);
  EssentialBCs bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  info("ndof: %d", Space::get_num_dofs(&space));

  // Initialize the weak formulation
  int order_inc = 3; // This determines the increase of integration order 
                     // for the axisymmetric term containing 1/r. Default is 3.
  CustomWeakFormMagnetostatics wf(MAT_IRON_1, MAT_IRON_2, &mu_inv_iron, MAT_AIR, 
                                  MAT_COPPER, MU_VACUUM, CURRENT_DENSITY, order_inc);

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
  bool residual_as_function = false;
  if (!hermes2d.solve_newton(coeff_vec, &dp, solver, matrix, rhs, 
			     NEWTON_TOL, NEWTON_MAX_ITER, verbose, residual_as_function, 
                             NEWTON_DAMPING)) error("Newton's iteration failed.");

  // Translate the resulting coefficient vector into the Solution sln.
  Solution::vector_to_solution(coeff_vec, &space, &sln);

  // Cleanup.
  delete [] coeff_vec;
  delete matrix;
  delete rhs;
  delete solver;

  // Visualise the solution and mesh.
  //ScalarView s_view1("Solution (vector potencial)", new WinGeom(0, 0, 350, 450));
  //s_view1.show_mesh(false);
  //s_view1.show(&sln);

  //ScalarView s_view2("Gradient (flux density)", new WinGeom(360, 0, 350, 450));
  //MagFilter grad(Hermes::vector<MeshFunction *>(&sln, &sln), Hermes::vector<int>(H2D_FN_DX, H2D_FN_DY));
  //s_view2.show_mesh(false);
  //s_view2.show(&grad);

  //OrderView o_view("Mesh", new WinGeom(720, 0, 350, 450));
  //o_view.show(&space);

  double coord_x[4] = {0.03, 0.06, 0.06, 0.03};
  double coord_y[4] = {0, 0, 0.05, 0.05};

  info("Coordinate (0.03,    0) value = %lf", sln.get_pt_value(coord_x[0], coord_y[0]));
  info("Coordinate (0.06,    0) value = %lf", sln.get_pt_value(coord_x[1], coord_y[1]));
  info("Coordinate (0.06, 0.05) value = %lf", sln.get_pt_value(coord_x[2], coord_y[2]));
  info("Coordinate (0.03, 0.05) value = %lf", sln.get_pt_value(coord_x[3], coord_y[3]));

  double t_value[4] = {0.031507, 0.012295, 0.026515, 0.153505};
  bool success = true;

  for (int i = 0; i < 4; i++)
  {
    if (fabs(t_value[i] - sln.get_pt_value(coord_x[i], coord_y[i])) > 1E-6) success = false;
  }

  if (success) {  // should pass with NEWTON_MAX_ITER = 9 and fail with NEWTON_MAX_ITER = 8
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }




  // Wait for all views to be closed.
  View::wait();
  return 0;
}

