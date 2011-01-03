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
  double norm = 0;
  double* product = new double[length];
  mat->multiply(vec, product);
  for (int i=0; i<length; i++) norm += vec[i]*product[i];
  delete [] product;
  return norm;
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
  ((UMFPackMatrix*)matrix_M_ref)->multiply(coeff_vec_ref, my_vec);


  // Debug.
  //info("coeff_vec_ref:");
  //for (int i=0; i<ndof_ref; i++) printf("%g ", coeff_vec_ref[i]);
  //printf("\n");
  //info("my_vec:");
  //for (int i=0; i<ndof_ref; i++) printf("%g ", my_vec[i]);
  //printf("\n");



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


  // Debug.
  //info("old Ap:");
  //for (int i=0; i<size+1; i++) printf("%d ", ap[i]);
  //printf("\n");
  //info("new Ap:");
  //for (int i=0; i<new_size+1; i++) printf("%d ", new_Ap[i]);
  //printf("\n");




  // Fill the new Ai array.
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
  ((UMFPackMatrix*)matrix_S_ref)->multiply(coeff_vec_ref, vec_SY);
  // Multiply M times Y.
  double* vec_MY = new double[ndof_ref];
  ((UMFPackMatrix*)matrix_M_ref)->multiply(coeff_vec_ref, vec_MY);
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



