#define HERMES_REPORT_ALL
#include "hermes2d.h"
#include <stdio.h>

// This test makes sure that example 50-eigenvalue works correctly.

int NUMBER_OF_EIGENVALUES = 1;                    // Desired number of eigenvalues.
int P_INIT = 4;                                   // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 0;                       // Number of initial mesh refinements.
double TARGET_VALUE = 2.0;                        // PySparse parameter: Eigenvalues in the vicinity of this number will be computed. 
double TOL = 1e-10;                               // Pysparse parameter: Error tolerance.
int MAX_ITER = 1000;                              // PySparse parameter: Maximum number of iterations.

using Teuchos::RCP;
using Teuchos::rcp;
using Hermes::EigenSolver;

// Boundary markers.
const int BDY_ALL = 1;

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  info("Desired number of eigenvalues: %d.", NUMBER_OF_EIGENVALUES);

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinements (optional).
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Enter boundary markers. 
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_ALL);

  // Enter Dirichlet boundary values.
  BCValues bc_values;
  bc_values.add_zero(BDY_ALL);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof: %d.", ndof);

  // Initialize the weak formulation for the left hand side, i.e., H.
  WeakForm wf_left, wf_right;
  wf_left.add_matrix_form(callback(bilinear_form_left));
  wf_right.add_matrix_form(callback(bilinear_form_right));

  // Initialize matrices.
  RCP<SparseMatrix> matrix_left = rcp(new CSCMatrix());
  RCP<SparseMatrix> matrix_right = rcp(new CSCMatrix());

  // Assemble the matrices.
  bool is_linear = true;
  DiscreteProblem dp_left(&wf_left, &space, is_linear);
  dp_left.assemble(matrix_left.get());
  DiscreteProblem dp_right(&wf_right, &space, is_linear);
  dp_right.assemble(matrix_right.get());

  EigenSolver es(matrix_left, matrix_right);
  info("Calling Pysparse...");
  es.solve(NUMBER_OF_EIGENVALUES, TARGET_VALUE, TOL, MAX_ITER);
  info("Pysparse finished.");
  es.print_eigenvalues();

  // Initializing solution vector, solution and ScalarView.
  double* coeff_vec;
  Solution sln;

  // Reading solution vectors and visualizing.
  double* eigenval = new double[NUMBER_OF_EIGENVALUES];
  int neig = es.get_n_eigs();
  if (neig != NUMBER_OF_EIGENVALUES) error("Mismatched number of eigenvectors in the eigensolver output file.");  
  for (int ieig = 0; ieig < neig; ieig++) {
    eigenval[ieig] = es.get_eigenvalue(ieig);
    int n;
    es.get_eigenvector(ieig, &coeff_vec, &n);
    // Convert coefficient vector into a Solution.
    Solution::vector_to_solution(coeff_vec, &space, &sln);
  }  

  info("ndof = %d", ndof);
  info("Coordinate ( 0.5, 0.5) value = %lf", sln.get_pt_value(0.5, 0.5));
  info("Coordinate ( 1.0, 0.5) value = %lf", sln.get_pt_value(1.0, 0.5));
  info("Coordinate ( 1.5, 0.5) value = %lf", sln.get_pt_value(1.5, 0.5));
  info("Coordinate ( 2.0, 0.5) value = %lf", sln.get_pt_value(2.0, 0.5));

  double coor_x[4] = {0.5, 1.0, 1.5, 2.0};
  double coor_y = 0.5;
  double t_value[4] = {0.146640, 0.257224, 0.304497, 0.277808};
  bool success = true;
  for (int i = 0; i < 4; i++)
  {
    double correct = t_value[i];
    double calculated = sln.get_pt_value(coor_x[i], coor_y);
    if (std::abs(correct - calculated) < 1E-6 ||
            std::abs(correct + calculated) < 1E-6) {
            // OK
    } else {
        info("Failed: i=%d, value[i]=%f, pt=%f, difference=%f", i,
                correct, calculated, std::abs(correct - calculated));
        success = false;
    }
  }
  if (success) {
    info("Success!");
    return ERR_SUCCESS;
  }
  else {
    info("Failure!");
    return ERR_FAILURE;
  }
}

