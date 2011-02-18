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

const int P_H = 0, P_V = 0;                       // Polynomial degrees of mesh elements in horizontal and vertical directions.
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

// Flux definition.
#include "fluxes.cpp"

// Boundary conditions.
BCType bc_types(int marker)
{
  return BC_NATURAL;
}

// Essential (Dirichlet) boundary condition values.
template<typename Real, typename Scalar>
Scalar g(int ess_bdy_marker, Real x, Real y)
{
  if (ess_bdy_marker == 1) return 1; else return 0;
}

// Weak forms.
#include "forms.cpp"

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
  for (int i=0; i<INIT_REF; i++) mesh.refine_all_elements();
  
  // Create an L2 space with default shapeset.
  L2Space space(&mesh, bc_types, NULL, Ord2(P_H, P_V));
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d", ndof);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form));
  wf.add_vector_form(callback(linear_form));
  wf.add_matrix_form_surf(callback(bilinear_form_boundary), H2D_DG_BOUNDARY_EDGE);
  wf.add_vector_form_surf(callback(linear_form_boundary), H2D_DG_BOUNDARY_EDGE);
  wf.add_matrix_form_surf(callback(bilinear_form_interface), H2D_DG_INNER_EDGE);

  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, &space, is_linear);
  
  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
  
  // Initialize the preconditioner in the case of SOLVER_AZTECOO.
  if (matrix_solver == SOLVER_AZTECOO) 
  {
    ((AztecOOSolver*) solver)->set_solver(iterative_method);
    ((AztecOOSolver*) solver)->set_precond(preconditioner);
    // Using default iteration parameters (see solver/aztecoo.h).
  }
    
  // Initialize the solution.
  Solution sln;
  
  // Assemble the stiffness matrix and right-hand side vector.
  info("Assembling the stiffness matrix and right-hand side vector.");
  dp.assemble(matrix, rhs);
  
  // Solve the linear system and if successful, obtain the solution.
  info("Solving the matrix problem.");
  if(solver->solve())
    Solution::vector_to_solution(solver->get_solution(), &space, &sln);
  else
    error ("Matrix solver failed.\n");
  
  // Time measurement.
  cpu_time.tick();
  
  // Clean up.
  delete solver;
  delete matrix;
  delete rhs;
  
  // Visualize the solution.
  ScalarView view1("Solution", new WinGeom(860, 0, 400, 350));
  view1.show(&sln);

  // Wait for all views to be closed.
  View::wait();
  
  return 0;
}

