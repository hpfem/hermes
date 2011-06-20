#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

//  This example solves a linear advection equation using Discontinuous Galerkin (DG) method.
//  Its purpose is to show how to evaluate surface matrix forms that take basis functions defined
//  on different elements.
//
//  PDE: \nabla \cdot (\Beta u) = 0, where \Beta = (-x_2, x_1) / |x| represents a circular counterclockwise flow field.
//
//  Domain: Square (0, 1)x(0, 1).
//
//  BC: Dirichlet,  u = g where \Beta(x) \cdot n(x) < 0; g = 1 on [0,0.5] x {0}, g = 0 anywhere else.
//
//  The following parameters can be changed:

const bool HERMES_VISUALIZATION = true;           // Set to "false" to suppress Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = true;              // Set to "true" to enable VTK output.
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
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);
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
  
  // VTK output.
  if (VTK_VISUALIZATION) 
  {
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
  if (HERMES_VISUALIZATION) 
  {
    Views::ScalarView<double> view("Solution", new Views::WinGeom(0, 0, 440, 350));
    view.show(&sln, Views::HERMES_EPS_VERYHIGH);
    Views::View::wait();
  }
  
  return 0;
}

