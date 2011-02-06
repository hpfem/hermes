#define HERMES_REPORT_ALL

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

// Calculate mass norm vec^T*mat*vec.
double calc_mass_product(UMFPackMatrix* mat, double* vec, int length)
{
  double result = 0;
  double* product = new double[length];
  mat->multiply_with_vector(vec, product);
  for (int i=0; i<length; i++) result += vec[i]*product[i];
  delete [] product;
  return result;
}

// Normalizes vector so that vec^T*mat*vec = 1. 
void normalize(UMFPackMatrix* mat, double* vec, int length) 
{
  double norm = sqrt(calc_mass_product(mat, vec, length));
  if (fabs(norm) < 1e-7) error("normalize(): Vector norm too small.");
  for (int i = 0; i < length; i++) vec[i] /= norm;
}

// Multiply two vectors.
double scalar_product(double* vec1, double* vec2, int length) 
{
  double val = 0;
  for (int i=0; i<length; i++) val += vec1[i] * vec2[i];
  return val;
}

void create_augmented_linear_system(SparseMatrix* matrix_S_ref, SparseMatrix* matrix_M_ref, 
                                    double* coeff_vec_ref, double lambda, UMFPackMatrix* new_matrix, 
                                    UMFPackVector* new_vector)
{
  // Extracting the arrays Ap, Ai and Ax from the matrices M and S on reference mesh.
  int size = ((UMFPackMatrix*)matrix_S_ref)->get_matrix_size();
  int ndof_ref = size;
  int nnz = ((UMFPackMatrix*)matrix_S_ref)->get_nnz();
  int* ap = ((UMFPackMatrix*)matrix_S_ref)->get_Ap();
  int* ai = ((UMFPackMatrix*)matrix_S_ref)->get_Ai();
  double* ax_S = ((UMFPackMatrix*)matrix_S_ref)->get_Ax();
  double* ax_M = ((UMFPackMatrix*)matrix_M_ref)->get_Ax();

  // Calculating the vector MY.
  double* my_vec = new double[size+1];
  ((UMFPackMatrix*)matrix_M_ref)->multiply_with_vector(coeff_vec_ref, my_vec);

  // Construct the augmented matrix for Newton's method.
  int new_size =  size + 1;
  int new_nnz = nnz + 2*size;
  int* new_Ap = new int[new_size+1]; assert(new_Ap != NULL);
  int* new_Ai = new int[new_nnz];    assert(new_Ai != NULL);
  double* new_Ax = new double[new_nnz]; assert(new_Ax != NULL);

  // Fill the new Ap array.
  new_Ap[0] = ap[0];
  for (int i=1; i < size+1; i++) new_Ap[i] = ap[i] + i;
  new_Ap[size+1] = new_Ap[size] + size;

  // Fill the new Ai array.
  int count = 0;
  for (int j=0; j < size; j++) {                                // Index of a column.
    for (int i = ap[j]; i < ap[j + 1]; i++) {                   // Index of a row.
      new_Ai[count++] = ai[i];                                  // First size entries in each column are the same.
    }
    new_Ai[count++] = size;                                     // Accounting for last item in columns 0, 1, size-1.
  }
  for (int i=0; i < size; i++) new_Ai[count++] = i;             // Accounting for last column.

  // Fill the new Ax array.  
  //double max = 0;
  count = 0;
  for (int j=0; j < size; j++) {                                // Index of a column.
    for (int i = ap[j]; i < ap[j + 1]; i++) {                   // Index of a row.
      new_Ax[count++] = ax_S[i] - lambda * ax_M[i];      // Block S minus lambda M
      //if (fabs(ax_S[i] - lambda * ax_M[i]) > 1e-10)
      //  warn("value[%d] = %.12f", i, ax_S[i] - lambda * ax_M[i]);
      //if (fabs(ax_S[i] - lambda * ax_M[i]) > max) max = fabs(ax_S[i] - lambda * ax_M[i]);
    }
    new_Ax[count++] = 2*my_vec[j];                              // 2*M*Y transposed is in the last row.
  }
  for (int i=0; i < size; i++) new_Ax[count++] = -my_vec[i];    // Minus M*Y is in the last column.
  //warn("max = %.12f", max);
  
  // Creating the new matrix.
  new_matrix->create(new_size, new_nnz, new_Ap, new_Ai, new_Ax);

  // Create the residual vector.
  // Multiply S times Y.
  double* vec_SY = new double[ndof_ref];
  ((UMFPackMatrix*)matrix_S_ref)->multiply_with_vector(coeff_vec_ref, vec_SY);
  // Multiply M times Y.
  double* vec_MY = new double[ndof_ref];
  ((UMFPackMatrix*)matrix_M_ref)->multiply_with_vector(coeff_vec_ref, vec_MY);
  // Calculate SY - lambda MY.
  double* residual = new double[ndof_ref+1];
  for (int i=0; i<ndof_ref; i++) residual[i] = vec_SY[i] - lambda * vec_MY[i];
  // Last component is Y^T M Y - 1.
  residual[ndof_ref] = calc_mass_product((UMFPackMatrix*)matrix_M_ref, coeff_vec_ref, ndof_ref) - 1;
  delete [] vec_SY;
  delete [] vec_MY;

  // Residual is multiplied with -1.
  for (int i=0; i<ndof_ref+1; i++) residual[i] *= -1.0;

  // Creating UMFPackVector for the residual.
  for (int i=0; i<ndof_ref+1; i++) new_vector->set(i, residual[i]);
}

bool solve_newton_eigen(Space* ref_space, UMFPackMatrix* matrix_S_ref, UMFPackMatrix* matrix_M_ref, 
                        double* coeff_vec_ref, double &lambda, MatrixSolverType matrix_solver,
                        double newton_tol, int newton_max_iter)
{
  ScalarView sview("", new WinGeom(0, 0, 440, 350));
  OrderView oview("", new WinGeom(450, 0, 410, 350));
  int ndof_ref = matrix_M_ref->get_size();
  double newton_err_rel;
  UMFPackMatrix* matrix_augm = new UMFPackMatrix();  
  UMFPackVector* vector_augm = new UMFPackVector(ndof_ref+1);
  Solver* solver_augm = create_linear_solver(matrix_solver, matrix_augm, vector_augm);
  Solution ref_sln_prev;
  Solution::vector_to_solution(coeff_vec_ref, ref_space, &ref_sln_prev);
  bool success = true;
  int it = 1;
  do {
    // Check the number of iterations.
    if (it >= newton_max_iter) {
      success = false;
      info("Newton's iteration not successful, returning false.");
      break;
    }

    // Normalize the eigenvector.
    //normalize((UMFPackMatrix*)matrix_M_ref, coeff_vec_ref, ndof_ref);
    //info("Eigenvector mass product: %g", calc_mass_product((UMFPackMatrix*)matrix_M_ref, 
    //                                     coeff_vec_ref, ndof_ref));

    // Create the augmented matrix and vector for Newton.
    //info("Creating augmented matrix and vector on reference mesh.");
    create_augmented_linear_system(matrix_S_ref, matrix_M_ref, coeff_vec_ref, lambda, 
                                   matrix_augm, vector_augm);

    // Solve the augmented matrix problem.
    //info("Solving augmented problem on reference mesh.");
    if(!solver_augm->solve()) {
      info("Matrix solver failed.\n");
      success = false;
      break;
    }
    double* increment = solver_augm->get_solution();

    // Update the eigenfunction and eigenvalue.
    //info("Updating reference solution.");
    for (int i=0; i<ndof_ref; i++) coeff_vec_ref[i] += increment[i];
    lambda += increment[ndof_ref];

    // Calculate relative error of the increment.
    Solution ref_sln_new;
    Solution::vector_to_solution(coeff_vec_ref, ref_space, &ref_sln_new);
    newton_err_rel = calc_rel_error(&ref_sln_prev, &ref_sln_new, HERMES_H1_NORM) * 100;

    // Updating reference solution.
    ref_sln_prev.copy(&ref_sln_new);

    // Debug.
    //char title1[100];
    //sprintf(title1, "Newton's iteration %d", it);
    //sview.set_title(title1);
    //sview.show(&ref_sln_prev);
    //View::wait(HERMES_WAIT_KEYPRESS);

    info("---- Newton iter %d, ndof %d, eigenvalue: %.12f, newton_err_rel %g%%", 
         it++, ndof_ref, lambda, newton_err_rel);
    //info("Eigenvalue increment: %.12f", increment[ndof_ref]);
    //info("Eigenvalue: %.12f", lambda);
  }
  while (newton_err_rel > newton_tol);

  // Clean up.
  delete matrix_augm;
  delete vector_augm;
  delete solver_augm;

  return success;
}

// This method always converges to the eigenvalue closest to the value of the argument lambda. 
// This is possible because the spectrum of the problem is shifted in such a way that the sought 
// eigenvalue comes to be very close to the origin where the method tends to converge.
bool solve_picard_eigen(Space* ref_space, UMFPackMatrix* matrix_S_ref, UMFPackMatrix* matrix_M_ref, 
                        double* coeff_vec_ref, double &lambda, MatrixSolverType matrix_solver,
                        double picard_tol, int picard_max_iter)
{
  int ndof_ref = matrix_M_ref->get_size();
  double picard_err_rel;
  UMFPackVector* vec_lambda_MY = new UMFPackVector(ndof_ref);
  Solver* solver = create_linear_solver(matrix_solver, matrix_S_ref, vec_lambda_MY);
  double* vec_MY = new double[ndof_ref]; 
  Solution ref_sln_prev;
  Solution::vector_to_solution(coeff_vec_ref, ref_space, &ref_sln_prev);
  bool success = true;
  double shift = lambda;
  // Construct shifted matrx.
  double *Sx = ((UMFPackMatrix*)matrix_S_ref)->get_Ax();
  double *Mx = ((UMFPackMatrix*)matrix_M_ref)->get_Ax();
  for (unsigned int i=0; i<((UMFPackMatrix*)matrix_S_ref)->get_nnz(); i++) Sx[i] = Sx[i] + shift * Mx[i];
  // Normalize the eigenvector.
  normalize((UMFPackMatrix*)matrix_M_ref, coeff_vec_ref, ndof_ref);
  // Init the eigenvalue for the shifted problem.
  lambda = calc_mass_product((UMFPackMatrix*)matrix_S_ref, coeff_vec_ref, ndof_ref)
             / calc_mass_product((UMFPackMatrix*)matrix_M_ref, coeff_vec_ref, ndof_ref);
  int it = 1;
  do {
    // Check the number of iterations.
    if (it >= picard_max_iter) {
      success = false;
      info("Picard's iteration not successful, returning false.");
      break;
    }
  
    matrix_M_ref->multiply_with_vector(coeff_vec_ref, vec_MY);
    for (int i=0; i<ndof_ref; i++) vec_lambda_MY->set(i, lambda*vec_MY[i]);

    // Solve the matrix problem.
    if(!solver->solve()) {
      info("Matrix solver failed.\n");
      success = false;
      break;
    }
    double* new_eigen_vec = solver->get_solution();

    // Copy the new eigen vector to coeff_vec_ref.
    for (int i=0; i<ndof_ref; i++) coeff_vec_ref[i] = new_eigen_vec[i];

    // Normalize the eigenvector.
    normalize((UMFPackMatrix*)matrix_M_ref, coeff_vec_ref, ndof_ref);

    // Update the eigenvalue.
    lambda = calc_mass_product((UMFPackMatrix*)matrix_S_ref, coeff_vec_ref, ndof_ref)
             / calc_mass_product((UMFPackMatrix*)matrix_M_ref, coeff_vec_ref, ndof_ref);

    // Calculate relative error of the increment.
    Solution ref_sln_new;
    Solution::vector_to_solution(coeff_vec_ref, ref_space, &ref_sln_new);
    picard_err_rel = calc_rel_error(&ref_sln_prev, &ref_sln_new, HERMES_H1_NORM) * 100;

    // Updating reference solution.
    ref_sln_prev.copy(&ref_sln_new);
    
    info("---- Picard iter %d, ndof %d, eigenvalue: %.12f, picard_err_rel %g%%", 
         it++, ndof_ref, lambda-shift, picard_err_rel);
  }
  while (picard_err_rel > picard_tol);
  
  // Unshift lambda
  lambda = lambda-shift;
  return success;

}
