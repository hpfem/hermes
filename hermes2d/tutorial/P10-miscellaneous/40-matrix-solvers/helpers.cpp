// Max row length in input file.
#define MAX_ROW_LEN	  1024

#ifndef INVALID_IDX
#define INVALID_IDX       ((unsigned int) -1)
#endif

// Matrix entry (coordinate format).
struct MatrixEntry {
  MatrixEntry() { }
  MatrixEntry(int m, int n, scalar value) {
  this->m = m;
  this->n = n;
  this->value = value;
  }
  void set(int m, int n, scalar value) {
    this->m = m;
    this->n = n;
    this->value = value;
    //info("Setting matrix element (%d, %d, %g).", m, n, value);
  }

  int m, n, id, used;
  scalar value;
};

// Vector entry (coordinate format).
struct VectorEntry {
  VectorEntry() { }
  VectorEntry(int m, scalar value) {
  this->m = m;
  this->value = value;
  }
  void set(int m, scalar value) {
    this->m = m;
    this->value = value;
    //info("Setting vector entry (%d, %g).", m, value);
  }

  int m, id, used;
  scalar value;
};

// Reads n numbers from the input file.
bool read_n_numbers(char *row, int n, double values[]) {
  int i = 0;
  char delims[] = " \t\n\r";
  char *token = strtok(row, delims);
  while (token != NULL && i < n) {
    int n;
    sscanf(token, "%d", &n);
    values[i++] = n;

    token = strtok(NULL, delims);
  }

  return (i == n);
}

// Processes the input file. Matrix and vector will be stored in 
// arrays of MatrixEntry and VectorEntry types, respectively.
bool read_matrix_and_rhs(char *file_name, int &n, 
                        Array<MatrixEntry> &mat, Array<VectorEntry> &rhs) 
{
  FILE *file = fopen(file_name, "r");
  if (file == NULL) return false;

  enum EState {
    STATE_N,
    STATE_MATRIX,
    STATE_RHS,
  } state = STATE_N;

  double buffer[3];
  char row[MAX_ROW_LEN];
  while (fgets(row, MAX_ROW_LEN, file) != NULL) {
    switch (state) {
      case STATE_N:
        if (read_n_numbers(row, 1, buffer)) {
          n = (int) buffer[0];
          state = STATE_MATRIX;
        } 
      break;

      case STATE_MATRIX:
        if (read_n_numbers(row, 3, buffer)) {
          MatrixEntry* me = mat.add();
          me->set((int) buffer[0], (int) buffer[1], buffer[2]);
        }
	else
        state = STATE_RHS;
      break;

      case STATE_RHS:
        if (read_n_numbers(row, 2, buffer)) {
          VectorEntry* ve = rhs.add();
          ve->set((int) buffer[0], (scalar) buffer[1]);
        }
      break;
    }
  }

  fclose(file);

  return true;
}

// Translates the matrix and vector from coordinate format 
// to a solver-specific format.
void build_matrix(int n, Array<MatrixEntry> &ar_mat, 
                  Array<VectorEntry> &ar_rhs, SparseMatrix *mat,
                  Vector *rhs) {
  mat->prealloc(n);
  for (int i = 0; i < ar_mat.get_size(); i++) {
    MatrixEntry &me = ar_mat.get(i);
    mat->pre_add_ij(me.m, me.n);
  }

  mat->alloc();
  for (int i = 0; i < ar_mat.get_size(); i++ ) {
    MatrixEntry &me = ar_mat.get(i);
    mat->add(me.m, me.n, me.value);
  }
  mat->finish();

  rhs->alloc(n);
  for (int i = 0; i < ar_rhs.get_size(); i++ ) {
    rhs->add(ar_rhs[i].m, ar_rhs[i].value);
  }
  rhs->finish();
}

// Block version of build_matrix().
void build_matrix_block(unsigned int n, Array<MatrixEntry> &ar_mat, Array<VectorEntry> &ar_rhs,
                        SparseMatrix *matrix, Vector *rhs) {
  matrix->prealloc(n);
  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = 0; j < n; j++)
      matrix->pre_add_ij(i, j);

  matrix->alloc();
  scalar **mat = new_matrix<scalar>(n, n);
  int *cols = new int[n];
  int *rows = new int[n];
  for (unsigned int i = 0; i < n; i++) {
    cols[i] = i;
    rows[i] = i;
  }
  for (int i = 0; i < ar_mat.get_size(); i++) {
    MatrixEntry &me = ar_mat.get(i);
    mat[me.m][me.n] = me.value;
  }
  matrix->add(n, n, mat, rows, cols);
  matrix->finish();

  rhs->alloc(n);
  scalar *rs = new scalar[n];
  for (int i = 0; i < ar_rhs.get_size(); i++) {
    VectorEntry &ve = ar_rhs.get(i);
    rs[ve.m] = ve.value;
  }

  unsigned int *u_rows = new unsigned int[n];
  for (unsigned int i = 0; i < n; i++) {
    u_rows[i] = i;
  }

  rhs->add(n, u_rows, rs);
  rhs->finish();
}

