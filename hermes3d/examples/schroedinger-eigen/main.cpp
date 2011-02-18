#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "config.h"
#include "hermes3d.h"
#include <stdio.h>

//  This example solves the eigenproblem for the time-independent Schroedinger 
//  equation in a cube with zero boundary conditions. Python and Pysparse must
//  be installed. 
//
//  PDE: -Laplace u + V(x,y,z) u = lambda_k u,
//  where lambda_0, lambda_1, ... are the eigenvalues.
//
//  Domain: Cube (-pi/2, pi/2)^3.
//
//  BC:  Homogeneous Dirichlet.
//
//  The eigenvalues are of the form 
//  i^2+j^2+k^2, where i,j,k are natural numbers 
//  The following parameters can be changed:

using Teuchos::ptr;
using Teuchos::RCP;
using Teuchos::rcp;
using Hermes::EigenSolver;

int NUMBER_OF_EIGENVALUES = 1;                    // Desired number of eigenvalues.

int P_INIT_X = 4;                                 // Uniform polynomial degree of mesh elements.
int P_INIT_Y = 4;                                 // Uniform polynomial degree of mesh elements.
int P_INIT_Z = 4;                                 // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 3;                       // Number of initial mesh refinements.
double TARGET_VALUE = 3.0;                        // PySparse parameter: Eigenvalues in the vicinity of this number will be computed. 
double TOL = 1e-10;                               // Pysparse parameter: Error tolerance.
int MAX_ITER = 1000;                              // PySparse parameter: Maximum number of iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Boundary condition types.
// Note: "essential" means that solution value is prescribed.
BCType bc_types(int marker)
{
  if (marker > 0) return BC_ESSENTIAL;
    else
      return BC_NATURAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
  return 0;
}

// Potential.
double V( double x, double y, double z){
  return 0.0; 
  //return -1./sqrt(x*x + y*y + z*z);
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  info("Desired number of eigenvalues: %d", NUMBER_OF_EIGENVALUES);
  // Time measurement.
  TimePeriod cpu_time;
  // Load the mesh.
  info("Loading mesh...");
  Mesh mesh;
  H3DReader mloader;
  cpu_time.reset();
  mloader.load("hexahedron.mesh3d", &mesh);
  cpu_time.tick();
  info("Total running time for loading mesh: %g s", cpu_time.accumulated());

  // Perform initial mesh refinements (optional).
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, Ord3(P_INIT_X, P_INIT_Y, P_INIT_Z));
  int ndof = Space::get_num_dofs(&space);
  info("ndof: %d", ndof);
  // Initialize the weak formulation for the left hand side, i.e., H.
  info("Initializing weak form...");
  WeakForm wf_left, wf_right;
  wf_left.add_matrix_form(bilinear_form_left, bilinear_form_left_ord, HERMES_SYM, HERMES_ANY );
  wf_right.add_matrix_form(callback(bilinear_form_right), HERMES_SYM, HERMES_ANY );

  // Initialize matrices and matrix solver.
  RCP<SparseMatrix> matrix_left = rcp(new CSCMatrix());
  RCP<SparseMatrix> matrix_right = rcp(new CSCMatrix());
  Solver* solver = create_linear_solver(matrix_solver, matrix_left.get());

  // Assemble the matrices.
  cpu_time.reset();
  info("Assembling matrices...");
  bool is_linear = true;
  DiscreteProblem dp_left(&wf_left, &space, is_linear);
  dp_left.assemble(matrix_left.get());
  DiscreteProblem dp_right(&wf_right, &space, is_linear);
  dp_right.assemble(matrix_right.get());
  cpu_time.tick();
  info("Total running time for assembling matrices: %g s", cpu_time.accumulated());
  cpu_time.reset();

  // Initialize eigensolver.
  cpu_time.reset();
  EigenSolver es(matrix_left, matrix_right);
  cpu_time.tick();
  info("Total running time for initializing EigenSolver : %g s.", cpu_time.accumulated());


  // Calling Python eigensolver. Solution will be written to "eivecs.dat".
  cpu_time.reset();
  info("Using eigensolver...");
  es.solve(NUMBER_OF_EIGENVALUES, TARGET_VALUE, TOL, MAX_ITER);
  info("Total running time for solving generalized eigenvalue problem: %g s", cpu_time.accumulated());
  es.print_eigenvalues();
  // Initializing solution vector, solution and ScalarView.
  info("Initializing solution vector...");
  double* coeff_vec;
  Solution sln(space.get_mesh());


  // Reading solution vectors from file and visualizing.
  int neig = es.get_n_eigs();
  if (neig != NUMBER_OF_EIGENVALUES) error("Mismatched number of eigenvectors in the eigensolver output file.");  
  for (int ieig = 0; ieig < neig; ieig++) {
    int n;
    es.get_eigenvector(ieig, &coeff_vec, &n);
    // Convert coefficient vector into a Solution.
    Solution::vector_to_solution(coeff_vec, &space, &sln);

    out_fn_vtk(&sln, "sln", ieig );
  }  

  // Clean up.
  delete solver;
  
  return 0; 
};

