#define HERMES_REPORT_INFO
#include "hermes2d.h"
#include <stdio.h>

//  This example shows how one can perform adaptivity to a selected eigenfunction
//  without calling the eigensolver again in each adaptivity step. The eigensolver 
//  is only called once at the beginning. 
//
//  PDE: -Laplace u + V*u = lambda_k u,
//  where lambda_0, lambda_1, ... are the eigenvalues.
//
//  Domain: Square (-pi/2, pi/2)^2.
//
//  BC:  Homogeneous Dirichlet.
//
//  The following parameters can be changed:

int TARGET_EIGENFUNCTION = 1;                     // Desired eigenfunction: 1 for the first, 2 for the second, etc.

int P_INIT = 2;                                   // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial mesh refinements.
double TARGET_VALUE = 2.0;                        // PySparse parameter: Eigenvalues in the vicinity of this number will be computed. 
double TOL = 1e-10;                               // Pysparse parameter: Error tolerance.
int MAX_ITER = 1000;                              // PySparse parameter: Maximum number of iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
double V(double x, double y) {
  return 0;
  //double r = sqrt(x*x + y*y);
  //return -1./(0.001 + r*r);
}

// Boundary markers.
const int BDY_BOTTOM = 1;
const int BDY_RIGHT = 2;
const int BDY_TOP = 3;
const int BDY_LEFT = 4;

// Weak forms.
#include "forms.cpp"

// Write the matrix in Matrix Market format.
void write_matrix_mm(const char* filename, Matrix* mat) 
{
  // Get matrix size.
  int ndof = mat->get_size();
  FILE *out = fopen(filename, "w" );
  if (out == NULL) error("failed to open file for writing.");

  // Calculate the number of nonzeros.
  int nz = 0;
  for (int i = 0; i < ndof; i++) {
    for (int j = 0; j <= i; j++) { 
      double tmp = mat->get(i, j);
      if (std::abs(tmp) > 1e-15) nz++;
    }
  } 

  // Write the matrix in MatrixMarket format
  fprintf(out,"%%%%MatrixMarket matrix coordinate real symmetric\n");
  fprintf(out,"%d %d %d\n", ndof, ndof, nz);
  for (int i = 0; i < ndof; i++) {
    for (int j = 0; j <= i; j++) { 
      double tmp = mat->get(i, j);
      if (std::abs(tmp) > 1e-15) fprintf(out, "%d %d %24.15e\n", i + 1, j + 1, tmp);
    }
  } 

  fclose(out);
}

// Normalizes vector so that vec^T*mat*vec = 1. 
void normalize(UMFPackMatrix* mat, double* vec, int length) 
{
  double norm = 0;
  double* product = new double[length];
  mat->multiply(vec, product);
  for (int i=0; i<length; i++) norm += vec[i]*product[i];
  norm = sqrt(norm);
  if (fabs(norm) < 1e-7) error("normalize(): Vector norm too small.");
  for (int i = 0; i < length; i++) vec[i] /= norm;
  delete [] product;
}

int main(int argc, char* argv[])
{
  info("Desired eigenfunction to calculate: %d.", TARGET_EIGENFUNCTION);

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinements (optional).
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(Hermes::Tuple<int>(BDY_BOTTOM, BDY_RIGHT, BDY_TOP, BDY_LEFT));

  // Enter Dirichlet boundary values.
  BCValues bc_values;
  bc_values.add_zero(Hermes::Tuple<int>(BDY_BOTTOM, BDY_RIGHT, BDY_TOP, BDY_LEFT));

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof: %d.", ndof);

  // Initialize the weak formulation for the left hand side, i.e., H.
  WeakForm wf_left, wf_right;
  wf_left.add_matrix_form(bilinear_form_left, bilinear_form_left_ord);
  wf_right.add_matrix_form(callback(bilinear_form_right));

  // Initialize matrices and matrix solver.
  SparseMatrix* matrix_left = create_matrix(matrix_solver);
  SparseMatrix* matrix_right = create_matrix(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix_left);

  // Assemble matrices M and S.
  bool is_linear = true;
  DiscreteProblem dp_left(&wf_left, &space, is_linear);
  dp_left.assemble(matrix_left);
  DiscreteProblem dp_right(&wf_right, &space, is_linear);
  dp_right.assemble(matrix_right);

  // Write matrix_left in MatrixMarket format.
  write_matrix_mm("mat_left.mtx", matrix_left);

  // Write matrix_left in MatrixMarket format.
  write_matrix_mm("mat_right.mtx", matrix_right);

  // Calling Python eigensolver. Solution will be written to "eivecs.dat".
  info("Calling Pysparse...");
  char call_cmd[255];
  sprintf(call_cmd, "python solveGenEigenFromMtx.py mat_left.mtx mat_right.mtx %g %d %g %d", 
	  TARGET_VALUE, TARGET_EIGENFUNCTION, TOL, MAX_ITER);
  system(call_cmd);
  info("Pysparse finished.");

  // Initializing solution vector, solution and ScalarView.
  double* coeff_vec = new double[ndof];

  // Reading solution vectors from file and visualizing.
  double* eigenval =new double[TARGET_EIGENFUNCTION];
  FILE *file = fopen("eivecs.dat", "r");
  char line [64];                  // Maximum line size.
  fgets(line, sizeof line, file);  // ndof
  int n = atoi(line);            
  if (n != ndof) error("Mismatched ndof in the eigensolver output file.");  
  fgets(line, sizeof line, file);  // Number of eigenvectors in the file.
  int neig = atoi(line);
  if (neig != TARGET_EIGENFUNCTION) error("Mismatched number of eigenvectors in the eigensolver output file.");  
  for (int ieig = 0; ieig < neig; ieig++) {
    // Get next eigenvalue from the file
    fgets(line, sizeof line, file);  // eigenval
    eigenval[ieig] = atof(line);            
    // Get the corresponding eigenvector.
    for (int i = 0; i < ndof; i++) {  
      fgets(line, sizeof line, file);
      coeff_vec[i] = atof(line);
    }
    // Normalize the eigenvector.
    normalize((UMFPackMatrix*)matrix_right, coeff_vec, ndof);
  }  
  fclose(file);

  /*** Now we have the desired eigenfunction in terms of the coefficient vector ***/

  // Convert coefficient vector into a Solution.
  Solution sln;
  double lambda = eigenval[neig-1];
  info("Eigenvalue (coarse mesh): %g", lambda);
  Solution::vector_to_solution(coeff_vec, &space, &sln);

  // Visualize the solution.
  ScalarView sview("", new WinGeom(0, 0, 440, 350));
  char title[100];
  sprintf(title, "Eigenfunction no. %d, val = %g", neig, lambda);
  sview.set_title(title);
  sview.show(&sln);
  OrderView oview("Mesh", new WinGeom(450, 0, 410, 350));
  oview.show(&space);

  // Wait for keypress.
  View::wait(HERMES_WAIT_KEYPRESS);

  /*** Begin adaptivity ***/

  // Construct globally refined reference mesh and setup reference space.
  Space* ref_space = construct_refined_space(&space);
  int ndof_ref = Space::get_num_dofs(ref_space);
  info("ndof_ref = %d", ndof_ref);

  // Initialize matrices and matrix solver on reference mesh.
  SparseMatrix* matrix_left_ref = create_matrix(matrix_solver);
  SparseMatrix* matrix_right_ref = create_matrix(matrix_solver);

  // Assemble matrices M and S on reference mesh.
  is_linear = true;
  DiscreteProblem dp_left_ref(&wf_left, ref_space, is_linear);
  dp_left_ref.assemble(matrix_left_ref);
  DiscreteProblem dp_right_ref(&wf_right, ref_space, is_linear);
  dp_right_ref.assemble(matrix_right_ref);




  // Debug: dumping matrices.
  //FILE* f;
  //f = fopen("s.txt", "w");
  //matrix_left_ref->dump(f, "mat-left");
  //fclose(f);
  //f = fopen("m.txt", "w");
  //matrix_left_ref->dump(f, "mat-right");
  //fclose(f);




  // Project the coarse mesh eigenfunction to the reference mesh.
  info("Projecting reference solution on coarse mesh.");
  double* coeff_vec_ref = new double[ndof_ref];
  OGProjection::project_global(ref_space, &sln, coeff_vec_ref, matrix_solver); 
  // Normalize the eigenvector.
  normalize((UMFPackMatrix*)matrix_right_ref, coeff_vec_ref, ndof_ref);

  // Extracting the arrays Ap, Ai and Ax from the matrices M and S on reference mesh.
  int size = ((UMFPackMatrix*)matrix_left_ref)->get_matrix_size();
  int nnz = ((UMFPackMatrix*)matrix_left_ref)->get_nnz();
  int* ap = ((UMFPackMatrix*)matrix_left_ref)->get_Ap();
  int* ai = ((UMFPackMatrix*)matrix_left_ref)->get_Ai();
  double* ax_left = ((UMFPackMatrix*)matrix_left_ref)->get_Ax();
  double* ax_right = ((UMFPackMatrix*)matrix_right_ref)->get_Ax();

  // Calculating the vector MY.
  double* my_vec = new double[size+1];
  ((UMFPackMatrix*)matrix_right_ref)->multiply(coeff_vec_ref, my_vec);


  // Debug.
  //info("coeff_vec_ref:");
  //for (int i=0; i<ndof_ref; i++) printf("%g ", coeff_vec_ref[i]);
  //printf("\n");
  //info("my_vec:");
  //for (int i=0; i<ndof_ref; i++) printf("%g ", my_vec[i]);
  //printf("\n");



  // Constructing the augmented matrix for Newton's method.
  int new_size =  ndof_ref + 1;
  int new_nnz = nnz + 2*size;
  int* new_Ap = new int[new_size+1]; assert(new_Ap != NULL);
  int* new_Ai = new int[new_nnz];    assert(new_Ai != NULL);
  double* new_Ax = new double[new_nnz]; assert(new_Ax != NULL);
  // Filling the new Ap array.
  new_Ap[0] = ap[0];
  for (int i=1; i < size+1; i++) new_Ap[i] = ap[i] + i;
  new_Ap[size+1] = new_Ap[size] + size;


  // Debug.
  //info("old Ap:");
  //for (int i=0; i<size+1; i++) printf("%d ", ap[i]);
  //printf("\n");
  //info("new Ap:");
  //for (int i=0; i<new_size+1; i++) printf("%d ", new_Ap[i]);
  //printf("\n");




  // Filling the new Ai array.
  int count = 0;
  for (int j=0; j < size; j++) {                                // Index of a column.
    for (int i = ap[j]; i < ap[j + 1]; i++) {                   // Index of a row.
      new_Ai[count++] = ai[i];                                  // First size entries in each column are the same.
    }
    new_Ai[count++] = size;                                     // Accounting for last item in columns 0, 1, size-1.
  }
  for (int i=0; i < size; i++) new_Ai[count++] = i;             // Accounting for last column.


  // Debug.
  //info("old Ai:");
  //for (int i=0; i<nnz; i++) printf("%d ", ai[i]);
  //printf("\n");
  //info("new Ai:");
  //for (int i=0; i<new_nnz; i++) printf("%d ", new_Ai[i]);
  //printf("\n");




  // Filling the new Ax array.  
  count = 0;
  for (int j=0; j < size; j++) {                                // Index of a column.
    for (int i = ap[j]; i < ap[j + 1]; i++) {                   // Index of a row.
      new_Ax[count++] = ax_left[i] - lambda*ax_right[i];        // Block S minus lambda M
    }
    new_Ax[count++] = 2*my_vec[j];                              // 2*M*Y transposed is in the last row.
  }
  for (int i=0; i < size; i++) new_Ax[count++] = -my_vec[i];    // Minus M*Y is in the last column.
  // Creating the new matrix.
  UMFPackMatrix* new_matrix = new UMFPackMatrix();
  new_matrix->create(new_size, new_nnz, new_Ap, new_Ai, new_Ax);

  // Create the residual vector.
  // Multiply S times Y.
  double* residual_1 = new double[ndof_ref];
  ((UMFPackMatrix*)matrix_left_ref)->multiply(coeff_vec_ref, residual_1);
  // Multiply M times Y.
  double* residual_2 = new double[ndof_ref];
  ((UMFPackMatrix*)matrix_right_ref)->multiply(coeff_vec_ref, residual_2);
  // Calculate SY - lambda MY.
  double* residual = new double[ndof_ref+1];
  for (int i=0; i<ndof_ref; i++) residual[i] = residual_1[i] - lambda * residual_2[i];
  // Last component is Y^T M Y - 1.
  double last_val = 0;
  for (int i=0; i<ndof_ref; i++) last_val += coeff_vec_ref[i] * residual_2[i];
  residual[ndof_ref] = last_val - 1;
  // Residual is multiplied with -1.
  for (int i=0; i<ndof_ref+1; i++)  residual[i] *= -1.0;

  // Creating UMFPackVector for the residual.
  Vector* new_vector = new UMFPackVector(ndof_ref+1);
  for (int i=0; i<ndof_ref+1; i++) new_vector->set(i, residual[i]);

  // Solving the matrix problem.
  Solver* solver_new = create_linear_solver(matrix_solver, new_matrix, new_vector);
  if(!solver_new->solve()) error ("Matrix solver failed.\n");
  double* increment = solver_new->get_solution();
  // Update the eigenfunction and eigenvalue.
  for (int i=0; i<ndof_ref; i++) coeff_vec_ref[i] += increment[i];
  lambda += increment[ndof_ref];
  info("Eigenvalue (fine mesh): %g", lambda);

  // Show updated eigenfunction on the reference mesh.
  Solution* ref_sln = new Solution();
  Solution::vector_to_solution(coeff_vec_ref, ref_space, ref_sln);
  sprintf(title, "Eigenfunction no. %d, val = %g", neig, lambda);
  sview.set_title(title);
  sview.show(ref_sln);
  sprintf(title, "Reference mesh");
  oview.set_title(title);
  oview.show(ref_space);






  // Wait for keypress.
  View::wait(HERMES_WAIT_KEYPRESS);



  // Visualize the increment.
  DiffFilter increm(Hermes::Tuple<MeshFunction*>(&sln, ref_sln));
  sprintf(title, "Newton's increment");
  sview.set_title(title);
  sview.show(&increm);





  // Wait for keypress.
  View::wait(HERMES_WAIT_KEYPRESS);







  // Cleaning up.
  delete [] coeff_vec;
  delete [] coeff_vec_ref;

  return 0; 
};

