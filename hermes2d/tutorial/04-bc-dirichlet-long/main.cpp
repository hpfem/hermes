#define H2D_REPORT_INFO
#include "hermes2d.h"

// This is a long version of example 04-bc-dirichlet: function solve_linear() is not used.

int P_INIT = 2;                                   // Initial polynomial degree in all elements.
int INIT_REF_NUM = 2;                             // Number of initial uniform mesh refinements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, SOLVER_NOX, 
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_UMFPACK.

// Problem parameters.
double CONST_F = -4.0; 

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int marker, double x, double y)
{
  return (-CONST_F/4.0)*(x*x + y*y);
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);
  int ndof = get_num_dofs(&space);
  info("ndof = %d", ndof);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form));
  wf.add_vector_form(callback(linear_form));

  // Initialize the linear problem.
  bool is_linear = true;
  FeProblem fep(&wf, &space, is_linear);

  // Initialize matrix solver.
  Solver* solver;
  switch (matrix_solver) {
    case SOLVER_AMESOS: solver = new AmesosSolver("Amesos_Klu", &fep); info("Using Amesos"); break;
    case SOLVER_MUMPS: solver = new MumpsSolver(&fep); info("Using Mumps"); break;
    case SOLVER_NOX: solver = new NoxSolver(&fep); info("Using Nox"); break;
    case SOLVER_PARDISO: solver = new PardisoLinearSolver(&fep); info("Using Pardiso"); break;
    case SOLVER_PETSC: solver = new PetscLinearSolver(&fep); info("Using PETSc"); break;
    case SOLVER_UMFPACK: solver = new UMFPackLinearSolver(&fep); info("Using UMFPack"); break;
    default: error("Unknown matrix solver requested.");
  }

  // Solve the matrix problem.
  if (!solver->solve()) error ("Matrix solver failed.\n");

  // Extract solution vector.
  scalar* coeffs = solver->get_solution();

  // Convert coefficient vector into a Solution.
  Solution* sln = new Solution(&space, coeffs);

  // Destroy matrix solver.
  delete solver;

  // Visualize the solution.
  ScalarView view("Solution", new WinGeom(0, 0, 440, 350));
  view.show(sln);

  // Wait for the view to be closed.
  View::wait();
  return 0;
}
