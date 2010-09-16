#define H2D_REPORT_INFO
#include "hermes2d.h"

// This is a long version of example 05-bc-neumann: function solve_linear() is not used.

int P_INIT = 4;                                   // Initial polynomial degree in all elements.
int CORNER_REF_LEVEL = 12;                        // Number of mesh refinements towards the re-entrant corner.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, SOLVER_NOX, 
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_UMFPACK.

// Problem parameters.
double CONST_F = -1.0;                            // Right-hand side.
double CONST_GAMMA[3] = {-0.5, 1.0, -0.5};        // Outer normal derivative on Gamma_1,2,3.

// Boundary condition types.
// Note: "natural" boundary condition means that 
// the solution value is not prescribed.
BCType bc_types(int marker)
{
  return (marker == 4) ? BC_ESSENTIAL : BC_NATURAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return 0.0;
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
  mesh.refine_towards_vertex(3, CORNER_REF_LEVEL);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);
  int ndof = get_num_dofs(&space);
  info("ndof = %d", ndof);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form));
  wf.add_vector_form(callback(linear_form));
  wf.add_vector_form_surf(callback(linear_form_surf));

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

  // Compute and show gradient magnitude.
  // (Note that the gradient at the re-entrant
  // corner needs to be truncated for visualization purposes.)
  ScalarView gradview("Gradient", new WinGeom(450, 0, 400, 350));
  MagFilter grad(Tuple<MeshFunction *>(sln, sln), 
                 Tuple<int>(H2D_FN_DX, H2D_FN_DY));
  gradview.show(&grad);

  // Wait for the views to be closed.
  View::wait();
  return 0;
}
