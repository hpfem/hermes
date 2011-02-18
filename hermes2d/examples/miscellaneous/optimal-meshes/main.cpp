#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "function/function.h"

using namespace RefinementSelectors;


//  This example is part of a research project on the design of 
//  optimal meshes.
//
//  Solved PDE: - Laplace u + u = f
//
//  Domain: Square (-1, 1) x (-1, 1)
//  
//  Exact solution: u(x,y) = atan(K*x) 
//
//  BC: Neumann, given by exact solution.

int P_INIT = 2;                                   // Uniform polynomial degree of all mesh elements.
int UNIFORM_REF_LEVEL = 0;                        // Number of initial uniform mesh refinements.
const double THRESHOLD = 0.3;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                           // Adaptive strategy:
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO;          // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                                  // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                                  // See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 0.5;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.01;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // reference mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
double K = 10.0;

// Boundary markers.
const int BDY_LEFT_RIGHT = 1;
const int BDY_TOP_BOTTOM = 2;

// Exact solution.
#include "exact_solution.cpp"

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square_2_elem.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < UNIFORM_REF_LEVEL; i++) mesh.refine_all_elements();

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_neumann(Hermes::vector<int>(BDY_LEFT_RIGHT, BDY_TOP_BOTTOM));

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d", ndof);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form_vol));
  wf.add_vector_form(callback(linear_form_vol));
  wf.add_vector_form_surf(callback(linear_form_surf_right), 2);
  wf.add_vector_form_surf(callback(linear_form_surf_left), 2);

  Solution sln; 

  // NON-ADAPTIVE VERSION
  
  // Initialize the linear problem.
  bool is_linear = true;
  DiscreteProblem* dp = new DiscreteProblem(&wf, &space, is_linear);

  // Select matrix solver.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Assemble stiffness matrix and rhs.
  dp->assemble(matrix, rhs);
 
  // Solve the linear system of the reference problem. If successful, obtain the solutions.
  if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), &space, &sln);
  else error ("Matrix solver failed.\n");

  // Visualize the solution.
  ScalarView view("Solution", new WinGeom(0, 0, 440, 350));
  view.show(&sln);

  // Calculate error wrt. exact solution.
  Solution sln_exact;
  sln_exact.set_exact(&mesh, exact);
  double err = calc_abs_error(&sln, &sln_exact, HERMES_H1_NORM);
  printf("err = %g, err_squared = %g\n\n", err, err*err);
 
  // Wait for all views to be closed.
  View::wait();
  return 0;
}
