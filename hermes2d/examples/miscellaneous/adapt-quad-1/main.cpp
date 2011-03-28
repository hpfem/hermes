#define HERMES_REPORT_ALL
#include "hermes2d.h"

// This example shows how to turn on adaptive quadrature in the 
// evaluation of weak forms. It is derived from example P01-linear/03-poisson.
//
// PDE: Poisson equation -Laplace u = CONST_F with homogeneous (zero)
//      Dirichlet boundary conditions.
//

const bool ADAPTIVE_QUADRATURE = true;            // Evaluate weak forms using adaptive quadrature.
const bool HERMES_VISUALIZATION = true;           // Set to "false" to suppress Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = true;              // Set to "true" to enable VTK output.
const int P_INIT = 2;                             // Uniform polynomial degree of mesh elements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Boundary markers.
const std::string BDY_BOTTOM = "1", BDY_OUTER = "2", BDY_LEFT = "3", BDY_INNER = "4";

// Problem parameters.
const double CONST_F = 2.0;  

// Weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinements (optional).
  //mesh.refine_all_elements();

  // Initialize boundary conditions
  DefaultEssentialBCConst bc_essential(Hermes::vector<std::string>(BDY_BOTTOM, BDY_OUTER, BDY_LEFT, BDY_INNER), 0.0);
  EssentialBCs bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d", ndof);

  // Initialize the weak formulation. Not providing the order determination form 
  // (or callback) turns on adaptive numerical quadrature. The quadrature begins 
  // with using a first-order rule in the entire element. Then the element is split 
  // uniformly in space and the quadrature order is increased by "adapt_order_increase".
  // Then the form is calculated again by employing the new quadrature in subelements. 
  // This provides a more accurate result. If relative error is less than 
  // "adapt_rel_error_tol", the computation stops, otherwise the same procedure is 
  // applied recursively to all four subelements. 
  int adapt_order_increase = 1;
  double adapt_rel_error_tol = 1e1;
  WeakFormPoisson wf(CONST_F, ADAPTIVE_QUADRATURE, adapt_order_increase, adapt_rel_error_tol);
  
  if (ADAPTIVE_QUADRATURE)
    info("Adaptive quadrature ON.");    
  else
    info("Adaptive quadrature OFF.");    

  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, &space, is_linear);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Initialize the solution.
  Solution sln;

  // Assemble the stiffness matrix and right-hand side vector.
  info("Assembling the stiffness matrix and right-hand side vector.");
  dp.assemble(matrix, rhs);

  // Solve the linear system and if successful, obtain the solution.
  info("Solving the matrix problem.");
  if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), &space, &sln);
  else error ("Matrix solver failed.\n");

  // VTK output.
  if (VTK_VISUALIZATION) {
    // Output solution in VTK format.
    Linearizer lin;
    bool mode_3D = true;
    lin.save_solution_vtk(&sln, "sln.vtk", "Temperature", mode_3D);
    info("Solution in VTK format saved to file %s.", "sln.vtk");

    // Output mesh and element orders in VTK format.
    Orderizer ord;
    ord.save_orders_vtk(&space, "ord.vtk");
    info("Element orders in VTK format saved to file %s.", "ord.vtk");
  }

  // Visualize the solution.
  if (HERMES_VISUALIZATION) {
    ScalarView view("Solution", new WinGeom(0, 0, 440, 350));
    view.show(&sln);
    View::wait();
  }

  // Clean up.
  delete solver;
  delete matrix;
  delete rhs;

  return 0;
}

