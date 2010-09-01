#define H2D_REPORT_INFO
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

// This is a long version of example 08-system: function solve_linear() is not used.

const int P_INIT = 4;                                      // Initial polynomial degree of all elements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;           // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                           // SOLVER_MUMPS, and more are coming.

// Problem parameters.
const double E  = 200e9;                                   // Young modulus (steel).
const double nu = 0.3;                                     // Poisson ratio.
const double f_0  = 0;                                     // External force in x-direction.
const double f_1  = 1e4;                                   // External force in y-direction.
const double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));  // First Lame constant.
const double mu = E / (2*(1 + nu));                        // Second Lame constant.

// Boundary marker (external force).
const int GAMMA_3_BDY = 3;

// Boundary condition types.
BCType bc_types(int marker)
  { return (marker == 1) ? BC_ESSENTIAL : BC_NATURAL; }

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
  { return 0; }

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("sample.mesh", &mesh);

  // Perform uniform mesh refinement.
  mesh.refine_all_elements();

  // Create x- and y- displacement space using the default H1 shapeset.
  H1Space u_space(&mesh, bc_types, essential_bc_values, P_INIT);
  H1Space v_space(&mesh, bc_types, essential_bc_values, P_INIT);
  int ndof = get_num_dofs(Tuple<Space *>(&u_space, &v_space));
  info("ndof = %d.", ndof);

  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, callback(bilinear_form_0_0), H2D_SYM);  // Note that only one symmetric part is
  wf.add_matrix_form(0, 1, callback(bilinear_form_0_1), H2D_SYM);  // added in the case of symmetric bilinear
  wf.add_matrix_form(1, 1, callback(bilinear_form_1_1), H2D_SYM);  // forms.
  wf.add_vector_form_surf(0, callback(linear_form_surf_0), GAMMA_3_BDY);
  wf.add_vector_form_surf(1, callback(linear_form_surf_1), GAMMA_3_BDY);

  // Initialize the linear problem.
  LinearProblem lp(&wf, Tuple<Space *>(&u_space, &v_space));

  // Select matrix solver.
  Matrix* mat; Vector* rhs; CommonSolver* solver;
  init_matrix_solver(matrix_solver, ndof, mat, rhs, solver);

  // Assemble stiffness matrix and rhs.
  lp.assemble(mat, rhs);

  // Solve the matrix problem.
  if (!solver->solve(mat, rhs)) error ("Matrix solver failed.\n");

  // Convert coefficient vector into a Solution.
  Solution* u_sln = new Solution(&u_space, rhs);
  Solution* v_sln = new Solution(&v_space, rhs);

  // Visualize the solution.
  ScalarView view("Von Mises stress [Pa]", new WinGeom(0, 0, 800, 400));
  VonMisesFilter stress(Tuple<MeshFunction*>(u_sln, v_sln), lambda, mu);
  view.show_mesh(false);
  view.show(&stress, H2D_EPS_HIGH, H2D_FN_VAL_0, u_sln, v_sln, 1.5e5);

  // Wait for the view to be closed.
  View::wait();
  return 0;
}

