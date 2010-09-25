#define H2D_REPORT_INFO
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

// This example explains how to create multiple spaces over a mesh and use them
// to solve a simple problem of linear elasticity. Note how Tuples are used, 
// they replace variable-length argument lists. At the end, VonMises filter is 
// used to visualize the stress.
//
// PDE: Lame equations of linear elasticity.
//
// BC: du_1/dn = f_0 on Gamma_3 and du_1/dn = 0 on Gamma_2, Gamma_4, Gamma_5,
//     du_2/dn = f_1 on Gamma_3 and du_2/dn = 0 on Gamma_2, Gamma_4, Gamma_5,
//     u_1 = 0 and u_2 = 0 on Gamma_1.
//
// The following parameters can be changed:

const int P_INIT = 6;                                      // Initial polynomial degree of all elements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;           // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, SOLVER_NOX, 
                                                           // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_UMFPACK.

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
  info("ndof = %d.", get_num_dofs(Tuple<Space *>(&u_space, &v_space)));

  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, callback(bilinear_form_0_0), HERMES_SYM);  // Note that only one symmetric part is
  wf.add_matrix_form(0, 1, callback(bilinear_form_0_1), HERMES_SYM);  // added in the case of symmetric bilinear
  wf.add_matrix_form(1, 1, callback(bilinear_form_1_1), HERMES_SYM);  // forms.
  wf.add_vector_form_surf(0, callback(linear_form_surf_0), GAMMA_3_BDY);
  wf.add_vector_form_surf(1, callback(linear_form_surf_1), GAMMA_3_BDY);

  // Initialize the FE problem.
  bool is_linear = true;
  FeProblem fep(&wf, Tuple<Space *>(&u_space, &v_space), is_linear);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_solver(matrix_solver, matrix, rhs);

  // Initialize the solutions.
  Solution u_sln, v_sln;

  // Assemble the stiffness matrix and right-hand side vector.
  info("Assembling the stiffness matrix and right-hand side vector.");
  fep.assemble(matrix, rhs);

  // Solve the linear system and if successful, obtain the solutions.
  info("Solving the matrix problem.");
  if(solver->solve())
    vector_to_solutions(solver->get_solution(), Tuple<Space *>(&u_space, &v_space), Tuple<Solution *>(&u_sln, &v_sln));
  else
    error ("Matrix solver failed.\n");
  
  // Visualize the solution.
  ScalarView view("Von Mises stress [Pa]", new WinGeom(0, 0, 800, 400));
  VonMisesFilter stress(Tuple<MeshFunction *>(&u_sln, &v_sln), lambda, mu);
  view.show_mesh(false);
  view.show(&stress, HERMES_EPS_HIGH, H2D_FN_VAL_0, &u_sln, &v_sln, 1.5e5);

  // Wait for the view to be closed.
  View::wait();

  // Clean up.
  delete solver;
  delete matrix;
  delete rhs;

  return 0;
}

