#define HERMES_REPORT_ALL
#include "../definitions.h"

// This test makes sure that example 03-poisson works correctly.
// CAUTION: This test will fail when any changes to the shapeset
// are made, but it is easy to fix (see below).

const int P_INIT = 5;                             // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.
Hermes::MatrixSolverType matrix_solver = Hermes::SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
const double LAMBDA_AL = 236.0;            // Thermal cond. of Al for temperatures around 20 deg Celsius.
const double LAMBDA_CU = 386.0;            // Thermal cond. of Cu for temperatures around 20 deg Celsius.
const double VOLUME_HEAT_SRC = 5e2;        // Volume heat sources generated (for example) by electric current.        
const double FIXED_BDY_TEMP = 20.0;        // Fixed temperature on the boundary.

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes::Hermes2D::Global<double> hermes2d;

  // Load the mesh.
  Hermes::Hermes2D::Mesh mesh;
  Hermes::Hermes2D::H2DReader mloader;
  mloader.load("../domain.mesh", &mesh);

  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) 
    mesh.refine_all_elements();

  // Initialize the weak formulation.
  CustomWeakFormPoisson wf("Aluminum", new Hermes::Hermes2D::HermesFunction<double>(LAMBDA_AL), "Copper", 
    new Hermes::Hermes2D::HermesFunction<double>(LAMBDA_CU), new Hermes::Hermes2D::HermesFunction<double>(-VOLUME_HEAT_SRC));

  // Initialize essential boundary conditions.
  Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential(Hermes::vector<std::string>("Bottom", "Inner", "Outer", "Left"), 
    FIXED_BDY_TEMP);
  Hermes::Hermes2D::EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  Hermes::Hermes2D::H1Space<double> space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();
  info("ndof = %d", ndof);

  // Initialize the FE problem.
  Hermes::Hermes2D::DiscreteProblem<double> dp(&wf, &space);

  // Set up the solver, matrix, and rhs according to the solver selection.
  Hermes::Algebra::SparseMatrix<double>* matrix = create_matrix<double>(matrix_solver);
  Hermes::Algebra::Vector<double>* rhs = create_vector<double>(matrix_solver);
  Hermes::Solvers::Solver<double>* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Initial coefficient vector for the Newton's method.  
  double* coeff_vec = new double[ndof];
  memset(coeff_vec, 0, ndof*sizeof(double));

  // Perform Newton's iteration.
  if (!hermes2d.solve_newton(coeff_vec, &dp, solver, matrix, rhs, true, 1E-8, 100, true)) 
    error("Newton's iteration failed.");

  // Translate the resulting coefficient vector into a Solution.
  Hermes::Hermes2D::Solution<double> sln;
  Hermes::Hermes2D::Solution<double>::vector_to_solution(coeff_vec, &space, &sln);

  // Actual test. The values of 'sum' depend on the
  // current shapeset. If you change the shapeset,
  // you need to correct these numbers.
  double sum = 0;
  for (int i=0; i < ndof; i++) sum += coeff_vec[i];
  printf("coefficient sum = %g\n", sum);

  // Clean up.
  delete [] coeff_vec;
  delete solver;
  delete matrix;
  delete rhs;

  bool success = true;
  if (std::abs(sum + 0.357318) > 1e-4) success = false;

  if (success == true) 
  {
    printf("Success!\n");
    return TEST_SUCCESS;
  }
  else 
  {
    printf("Failure!\n");
    return TEST_FAILURE;
  }
}

