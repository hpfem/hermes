#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

// This test makes sure that example linear-advection-dg works correctly.

const int P_INIT = 0;                             // Polynomial degree of mesh elements.
const int INIT_REF = 1;                           // Number of initial uniform mesh refinements.
const char* iterative_method = "cg";              // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "jacobi";            // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers). 
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  // preconditioner from IFPACK (see hermes_common/solver/aztecoo.h).
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Boundary markers.
const std::string BDY_BOTTOM_LEFT = "1";

// Weak forms.
#include "../definitions.cpp"

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../square.mesh", &mesh);
  //mloader.load("square-tri.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF; i++) 
    mesh.refine_all_elements();
  
  // Create an L2 space with default shapeset.
  L2Space<double> space(&mesh, P_INIT);
  int ndof = Space<double>::get_num_dofs(&space);
  info("ndof = %d", ndof);

  // Initialize the weak formulation.
  // Initialize the weak formulation.
  CustomWeakForm wf(BDY_BOTTOM_LEFT);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, &space);
  
  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix<double>* matrix = create_matrix<double>(matrix_solver);
  Vector<double>* rhs = create_vector<double>(matrix_solver);
  Solver<double>* solver = create_linear_solver<double>(matrix_solver, matrix, rhs);
  
  // Initialize the preconditioner in the case of SOLVER_AZTECOO.
  if (matrix_solver == SOLVER_AZTECOO) 
  {
    dynamic_cast<AztecOOSolver<double>*>(solver)->set_solver(iterative_method);
    dynamic_cast<AztecOOSolver<double>*>(solver)->set_precond(preconditioner);
    // Using default iteration parameters (see solver/aztecoo.h).
  }
    
  // Initialize the solution.
  Solution<double> sln;
  
  // Assemble the stiffness matrix and right-hand side vector.
  info("Assembling the stiffness matrix and right-hand side vector.");
  dp.assemble(matrix, rhs);
  
  // Solve the linear system and if successful, obtain the solution.
  info("Solving the matrix problem.");
  if(solver->solve())
    Solution<double>::vector_to_solution(solver->get_solution(), &space, &sln);
  else
    error ("Matrix solver failed.\n");
  
  // Time measurement.
  cpu_time.tick();
  
  // Clean up.
  delete solver;
  delete matrix;
  delete rhs;
  
  info("ndof = %d", ndof);
  info("Coordinate ( 0.1, 0.1) value = %lf", sln.get_pt_value(0.1, 0.1));
  info("Coordinate ( 0.3, 0.3) value = %lf", sln.get_pt_value(0.3, 0.3));
  info("Coordinate ( 0.5, 0.5) value = %lf", sln.get_pt_value(0.5, 0.5));
  info("Coordinate ( 0.7, 0.7) value = %lf", sln.get_pt_value(0.7, 0.7));

  double coor_x_y[4] = {0.1, 0.3, 0.5, 0.7};
  double value[4] = {0.999885, 0.844340, 0.000000, 0.000000};
  
  bool success = true;
  for (int i = 0; i < 4; i++)
    if (abs(value[i] - sln.get_pt_value(coor_x_y[i], coor_x_y[i])) > 1E-6) 
      success = false;
  
  if (success) {
    printf("Success!\n");
    return TEST_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return TEST_FAILURE;
  }
}

