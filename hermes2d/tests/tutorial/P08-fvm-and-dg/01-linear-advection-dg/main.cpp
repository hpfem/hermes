#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

// This test makes sure that example 60-linear-advection-dg works correctly.

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
  
  info("ndof = %d", ndof);
  info("Coordinate ( 0.1, 0.1) value = %lf", sln.get_pt_value(0.1, 0.1));
  info("Coordinate ( 0.3, 0.3) value = %lf", sln.get_pt_value(0.3, 0.3));
  info("Coordinate ( 0.5, 0.5) value = %lf", sln.get_pt_value(0.5, 0.5));
  info("Coordinate ( 0.7, 0.7) value = %lf", sln.get_pt_value(0.7, 0.7));

  double coor_x_y[4] = {0.1, 0.3, 0.5, 0.7};
  double value[4] = {0.999885, 0.844340, 0.000000, 0.000000};
  bool success = true;
  for (int i = 0; i < 4; i++)
  {
    if (abs(value[i] - sln.get_pt_value(coor_x_y[i], coor_x_y[i])) > 1E-6) success = false;
  }
  if (success) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}

