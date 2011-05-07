#define HERMES_REPORT_ALL
#include <hermes3d.h>

// This example shows the Bridge model
// a hexahedral mesh in CTU format.
//
// NOTE: This is a temporary version in progress. After the 
// example works with the Poisson equation, we will extend it
// to linear elasticity. After that works, we will extend it 
// to an eigenproblem for linear elasticity. 
//
// PDE: Laplace equation -Laplace u = f, where f = CONST_F.

const int P_INIT_X = 1,
          P_INIT_Y = 1,
          P_INIT_Z = 1;                           // Initial polynomial degree of all mesh elements.
bool solution_output = true;                      // Generate output files (if true).

// Right-hand side.
const scalar CONST_F = -1.0;  

// Boundary markers.
int bdy_supported = 1;

MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.                                               

// Boundary condition types. 
BCType bc_types(int marker) 
{
  return (marker == bdy_supported) ? H3D_BC_ESSENTIAL : H3D_BC_NATURAL;
}

// Essential (Dirichlet) boundary condition values. 
scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
  return 0;
}

#include "definitions.cpp"

int main(int argc, char **args) 
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh. 
  Mesh mesh;

  //CTUReader mloader;
  H3DReader mloader;

  info("Loading mesh...");
  mloader.load("bridge.mesh3d", &mesh);

  // Create H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, Ord3(P_INIT_X, P_INIT_Y, P_INIT_Z));

  // Initialize weak formulation. 
  WeakForm wf;
  wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<Ord, Ord>, HERMES_SYM);
  wf.add_vector_form(linear_form<double, scalar>, linear_form<Ord, Ord>, HERMES_ANY_INT);

  // Initialize discrete problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, &space, is_linear);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Assemble stiffness matrix and load vector.
  info("Assembling the linear problem (ndof: %d).", Space::get_num_dofs(&space));
  dp.assemble(matrix, rhs);

  // Solve the linear system. If successful, obtain the solution.
  info("Solving the linear problem.");
  Solution sln(&mesh);
  if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), &space, &sln);
  else error ("Matrix solver failed.\n");

  // Output the solution for Paraview.
  if (solution_output) out_fn_vtk(&sln, "sln");

  // Time measurement.
  cpu_time.tick();

  // Print timing information.
  info("Solutions and mesh with polynomial orders saved. Total running time: %g s", 
       cpu_time.accumulated());

  // Clean up.
  delete matrix;
  delete rhs;
  delete solver;

  return 0;
}
