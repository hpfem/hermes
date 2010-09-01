#define H2D_REPORT_INFO
#include "hermes2d.h"

// This example explains how to use Newton boundary conditions. Again,
// a Filter is used to visualize the solution gradient.
//
// PDE: Laplace equation -Laplace u = 0 (this equation describes, among
// many other things, also stationary heat transfer in a homogeneous linear
// material).
//
// BC: u = T1 ... fixed temperature on Gamma_3 (Dirichlet)
//     du/dn = 0 ... insulated wall on Gamma_2 and Gamma_4 (Neumann)
//     du/dn = H*(u - T0) ... heat flux on Gamma_1 (Newton)
//
// Note that the last BC can be written in the form  du/dn - H*u = -h*T0.
//
// The following parameters can be changed:

int P_INIT = 6;                                   // Uniform polynomial degree of all mesh elements.
int UNIFORM_REF_LEVEL = 2;                        // Number of initial uniform mesh refinements.
int CORNER_REF_LEVEL = 12;                        // Number of mesh refinements towards the re-entrant corner.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Problem parameters.
double T1 = 30.0;            // Prescribed temperature on Gamma_3.
double T0 = 20.0;            // Outer temperature on Gamma_1.
double H  = 0.05;            // Heat flux on Gamma_1.

// Boundary markers.
const int NEWTON_BDY = 1;

// Boundary condition types.
BCType bc_types(int marker)
  { return (marker == 3) ? BC_ESSENTIAL : BC_NATURAL; }

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
  { return T1; }

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i=0; i<UNIFORM_REF_LEVEL; i++) mesh.refine_all_elements();
  mesh.refine_towards_vertex(3, CORNER_REF_LEVEL);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);
  int ndof = get_num_dofs(&space);
  info("ndof = %d", ndof);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form));
  wf.add_matrix_form_surf(callback(bilinear_form_surf), NEWTON_BDY);
  wf.add_vector_form_surf(callback(linear_form_surf), NEWTON_BDY);

  // Solve the linear problem.
  Solution sln;
  solve_linear(&space, &wf, matrix_solver, &sln);

  // Visualize the solution.
  ScalarView view("Solution", new WinGeom(0, 0, 440, 350));
  view.show(&sln);

  // Compute and show gradient magnitude.
  // (Note that the gradient at the re-entrant
  // corner needs to be truncated for visualization purposes.)
  ScalarView gradview("Gradient", new WinGeom(450, 0, 400, 350));
  MagFilter grad(Tuple<MeshFunction *>(&sln, &sln), 
                 Tuple<int>(H2D_FN_DX, H2D_FN_DY));
  gradview.show(&grad);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
