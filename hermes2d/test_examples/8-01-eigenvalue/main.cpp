#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"
#include <stdio.h>

using namespace RefinementSelectors;
using Teuchos::RCP;
using Teuchos::rcp;

//  This example solves a simple eigenproblem in a square. 
//  Python and Pysparse must be installed. 
//
//  PDE: -Laplace u + (x*x + y*y)u = lambda_k u,
//  where lambda_0, lambda_1, ... are the eigenvalues.
//
//  Domain: Square (0, pi)^2.
//
//  BC:  Homogeneous Dirichlet.
//
//  The following parameters can be changed:

const int NUMBER_OF_EIGENVALUES = 50;             // Desired number of eigenvalues.
const int P_INIT = 4;                             // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 3;                       // Number of initial mesh refinements.
const double TARGET_VALUE = 2.0;                  // PySparse parameter: Eigenvalues in the vicinity of 
                                                  // this number will be computed. 
const double TOL = 1e-10;                         // Pysparse parameter: Error tolerance.
const int MAX_ITER = 1000;                        // PySparse parameter: Maximum number of iterations.

int main(int argc, char* argv[])
{
  info("Desired number of eigenvalues: %d.", NUMBER_OF_EIGENVALUES);

  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinements (optional).
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize boundary conditions. 
  DefaultEssentialBCConst<double> bc_essential("Bdy", 0.0);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();
  info("ndof: %d.", ndof);

  // Initialize the weak formulation.
  WeakFormEigenLeft wf_left;
  WeakFormEigenRight wf_right;

  // Initialize matrices.
  RCP<SparseMatrix<double> > matrix_left = rcp(new CSCMatrix<double>());
  RCP<SparseMatrix<double> > matrix_right = rcp(new CSCMatrix<double>());

  // Assemble the matrices.
  DiscreteProblem<double> dp_left(&wf_left, &space);
  dp_left.assemble(matrix_left.get());
  DiscreteProblem<double> dp_right(&wf_right, &space);
  dp_right.assemble(matrix_right.get());

  EigenSolver<double> es(matrix_left, matrix_right);
  info("Calling Pysparse...");
  es.solve(NUMBER_OF_EIGENVALUES, TARGET_VALUE, TOL, MAX_ITER);
  info("Pysparse finished.");
  es.print_eigenvalues();

  // Initializing solution vector, solution and ScalarView.
  double* coeff_vec;
  Solution<double> sln;
  Views::ScalarView<double> view("Solution", new Views::WinGeom(0, 0, 440, 350));

  // Reading solution vectors and visualizing.
  double* eigenval = new double[NUMBER_OF_EIGENVALUES];
  int neig = es.get_n_eigs();
  if (neig != NUMBER_OF_EIGENVALUES) error("Mismatched number of eigenvectors in the eigensolver output file.");  
  for (int ieig = 0; ieig < neig; ieig++) {
    eigenval[ieig] = es.get_eigenvalue(ieig);
    int n;
    es.get_eigenvector(ieig, &coeff_vec, &n);
    // Convert coefficient vector into a Solution.
    Solution<double>::vector_to_solution(coeff_vec, &space, &sln);

    // Visualize the solution.
    char title[100];
    sprintf(title, "Solution %d, val = %g", ieig, eigenval[ieig]);
    view.set_title(title);
    view.show(&sln);

    // Wait for keypress.
    Views::View::wait(Views::HERMES_WAIT_KEYPRESS);
  }

  return 0; 
};

