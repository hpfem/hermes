#define HERMES_REPORT_ALL
#include "definitions.h"

// This example shows how to use Newton (Robin) boundary conditions.
// These conditions are used, for example, in heat transfer problems 
// when the heat transfer coefficient is known but the heat flux 
// itself is unknown.
//
// PDE: Poisson equation -div(LAMBDA grad u) - VOLUME_HEAT_SRC = 0.
//
// Boundary conditions: 
//   Boundary part "Outer": Newton condition LAMBDA * du/dn = ALPHA * (T_EXTERIOR - u).
//   Boundary parts "Bottom", "Inner" and "Left": 
//                          Dirichlet u(x, y) = A*x + B*y + C where the 
//                          constants A, B, C are arbitrary constants (defined below).
//
// The following parameters can be changed:

const bool HERMES_VISUALIZATION = true;           // Set to "false" to suppress Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = true;              // Set to "true" to enable VTK output.
const int P_INIT = 5;                             // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
const double LAMBDA_AL = 236.0;            // Thermal cond. of Al for temperatures around 20 deg Celsius.
const double LAMBDA_CU = 386.0;            // Thermal cond. of Cu for temperatures around 20 deg Celsius.
const double VOLUME_HEAT_SRC = 0.0;        // Volume heat sources generated by electric current.        
const double ALPHA = 5.0;                  // Heat transfer coefficient.
const double T_EXTERIOR = 50.0;            // Exterior temperature.
const double BDY_A_PARAM = 0.0;
const double BDY_B_PARAM = 0.0;
const double BDY_C_PARAM = 20.0;

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D<double> hermes2d;

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize the weak formulation.
  CustomWeakFormPoissonNewton wf("Aluminum", new HermesFunction<double>(LAMBDA_AL), 
                                 "Copper", new HermesFunction<double>(LAMBDA_CU), 
                                 new HermesFunction<double>(-VOLUME_HEAT_SRC),
                                  "Outer", ALPHA, T_EXTERIOR);
  
  // Initialize boundary conditions.
  CustomDirichletCondition bc_essential(Hermes::vector<std::string>("Bottom", "Inner", "Left"),
                                        BDY_A_PARAM, BDY_B_PARAM, BDY_C_PARAM);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();
  info("ndof = %d", ndof);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, &space);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix<double>* matrix = create_matrix<double>(matrix_solver);
  Vector<double>* rhs = create_vector<double>(matrix_solver);
  Solver<double>* solver = create_linear_solver<double>(matrix_solver, matrix, rhs);

  // Initial coefficient vector for the Newton's method.  
  double* coeff_vec = new double[ndof];
  memset(coeff_vec, 0, ndof*sizeof(double));

  // Perform Newton's iteration.
  if (!hermes2d.solve_newton(coeff_vec, &dp, solver, matrix, rhs)) error("Newton's iteration failed.");

  // Translate the resulting coefficient vector into the Solution sln.
  Solution<double> sln;
  Solution<double>::vector_to_solution(coeff_vec, &space, &sln);

  // VTK output.
  if (VTK_VISUALIZATION) {
    // Output solution in VTK format.
    Views::Linearizer<double> lin;
    bool mode_3D = true;
    lin.save_solution_vtk(&sln, "sln.vtk", "Temperature", mode_3D);
    info("Solution in VTK format saved to file %s.", "sln.vtk");

    // Output mesh and element orders in VTK format.
    Views::Orderizer ord;
    ord.save_orders_vtk(&space, "ord.vtk");
    info("Element orders in VTK format saved to file %s.", "ord.vtk");
  }

  // Visualize the solution.
  if (HERMES_VISUALIZATION) {
    Views::ScalarView<double> view("Solution", new Views::WinGeom(0, 0, 440, 350));
    view.show(&sln, Views::HERMES_EPS_VERYHIGH);
    Views::View::wait();
  }

  // Clean up.
  delete solver;
  delete matrix;
  delete rhs;
  delete [] coeff_vec;

  return 0;
}
