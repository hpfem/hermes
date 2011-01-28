#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "config.h"
//#include <getopt.h>
#include <hermes3d.h>

// Solving a simple heat equation to demonstrate how to use CUBIT with Hermes3D.
//
// Use mesh file 'cylinder2.e. Material IDs corresponds to elements markers, 
// sideset IDs correspond to face (BC) markers.
//
//  The following parameters can be changed:
const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.
const int P_INIT_X = 1,
          P_INIT_Y = 1,
          P_INIT_Z = 1;                           // Initial polynomial degree of all mesh elements.
bool solution_output = true;                      // Generate output files (if true).
const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "jacobi";            // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers). 
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  // preconditioner from IFPACK (see solver/aztecoo.h).
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Boundary condition types. 
BCType bc_types(int marker)
{
  if (marker == 1) return BC_ESSENTIAL;
  else return BC_NATURAL;
}

// Essential (Dirichlet) boundary condition values. 
scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
  return 10;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char **args)
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh. 
  Mesh mesh;
  ExodusIIReader mesh_loader;
  if (!mesh_loader.load("cylinder2.e", &mesh))
    error("Loading mesh file '%s' failed.\n", "cylinder2.e");

  // Perform initial mesh refinement.
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);
  
  // Create H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, Ord3(P_INIT_X, P_INIT_Y, P_INIT_Z));
  info("Number of DOF: %d.", Space::get_num_dofs(&space));

  // Initialize weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form1), HERMES_SYM, 1);
  wf.add_matrix_form(callback(bilinear_form2), HERMES_SYM, 2);
  wf.add_vector_form(callback(linear_form), HERMES_ANY);

  // Initialize discrete problem.
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

  // Assemble stiffness amtrix and load vector.
  info("Assembling the linear problem (ndof: %d).", Space::get_num_dofs(&space));
  dp.assemble(matrix, rhs);
	
  // Solve the linear system. If successful, obtain the solution.
  info("Solving the linear problem.");
  Solution sln(space.get_mesh());
  if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), &space, &sln);
  else error ("Matrix solver failed.\n");

  // Output solution and the boundary condition.
  if (solution_output) 
  {
    out_fn_vtk(&sln, "sln");
    out_bc_vtk(&mesh, "bc");
  }

  // Time measurement.
  cpu_time.tick();

  // Print timing information.
  info("Solution and the boundary condition saved. Total running time: %g s", cpu_time.accumulated());

  // Clean up.
  delete matrix;
  delete rhs;
  delete solver;

  return 0;
}
