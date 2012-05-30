#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE

#include "hermes_common.h"
#include <iostream>

using namespace Hermes::Algebra::DenseMatrixOperations;
using namespace Hermes::Solvers;

// Test of linear solvers.
// Read matrix and RHS from a file.

// Max row length in input file.
#define MAX_ROW_LEN	1024

class MatrixEntry
{
public:
  MatrixEntry() { }
  MatrixEntry(int m, int n, std::complex<double> value) {
    this->m = m;
    this->n = n;
    this->value = value;
  }

  int m, n;
  std::complex<double>  value;
};

void show_mat(const char *msg, std::map<unsigned int, MatrixEntry> mp)
{
  std::map<unsigned int, MatrixEntry>::iterator itr;

  std::cout << msg << std::endl;

  for(itr = mp.begin(); itr != mp.end(); ++itr)
    std::cout << " " << (int) itr->first << ": " <<
    (int) itr->second.m << " " <<
    (int) itr->second.n << " " <<
    (std::complex<double>) itr->second.value <<
    std::endl;

  std::cout << std::endl;
}

void show_rhs(const char *msg, std::map<unsigned int, std::complex<double> > mp) {
  std::map<unsigned int, std::complex<double> >::iterator itr;

  std::cout << msg << std::endl;

  for(itr = mp.begin(); itr != mp.end(); ++itr)
    std::cout << " " << (int) itr->first << ": " << (std::complex<double>) itr->second << std::endl;

  std::cout << std::endl;
}

bool testPrint(bool value, const char *msg, bool correct) {
  info("%s...", msg);
  if (value == correct) {
    info("OK.");
    return true;
  }
  else {
    info("failed.");
    return false;
  }
}

bool read_n_nums(char *row, int n, double values[]) {
  int i = 0;
  char delims[] = " \t\n\r";
  char *token = strtok(row, delims);
  while (token != NULL && i < n) {
    double entry_buffer;
    sscanf(token, "%lf", &entry_buffer);
    values[i++] = entry_buffer;

    token = strtok(NULL, delims);
  }

  return (i == n);
}

int read_matrix_and_rhs(char *file_name, int &n, int &nnz,
  std::map<unsigned int, MatrixEntry> &mat, std::map<unsigned int, std::complex<double> > &rhs, bool &cplx_2_real)
{
  FILE *file = fopen(file_name, "r");
  if (file == NULL) return TEST_FAILURE;

  enum EState {
    STATE_N,
    STATE_MATRIX,
    STATE_RHS,
    STATE_NNZ
  }
  state = STATE_N;

  // Variables needed to turn complex matrix into real.
  int k = 0;
  int l = 0;
  double* rhs_buffer = NULL;

  double buffer[4];
  char row[MAX_ROW_LEN];
  while (fgets(row, MAX_ROW_LEN, file) != NULL) {
    switch (state) {
    case STATE_N:
      if (read_n_nums(row, 1, buffer)) {
        if (cplx_2_real) {
          n = 2*((int) buffer[0]);
          rhs_buffer = new double[n];
          for (int i = 0; i < n; i++) {
            rhs_buffer[i] = 0.0;
          }
        }
        else
          n = (int) buffer[0];

        state = STATE_NNZ;
      }
      break;

    case STATE_NNZ:
      if (read_n_nums(row, 1, buffer))
        nnz = (int) buffer[0];

      state = STATE_MATRIX;
      break;

    case STATE_MATRIX:
      if (read_n_nums(row, 4, buffer)) {
        std::complex<double> cmplx_buffer(buffer[2], buffer[3]);
        mat[mat.size()] = (MatrixEntry((int) buffer[0], (int) buffer[1], (std::complex<double>) cmplx_buffer));
      }
      else
        state = STATE_RHS;
      break;

    case STATE_RHS:
      if (read_n_nums(row, 3, buffer)) {
        std::complex<double> cmplx_buffer(buffer[1], buffer[2]);
        rhs[(int) buffer[0]] = (std::complex<double>) cmplx_buffer;
      }
      break;
    }
  }

  fclose(file);

  // Free memory
  delete [] rhs_buffer;

  // Clear pointer.
  rhs_buffer = NULL;

  return TEST_SUCCESS;
}

void build_matrix(int n, std::map<unsigned int, MatrixEntry> &ar_mat, std::map<unsigned int, std::complex<double> > &ar_rhs,
  SparseMatrix<std::complex<double> > *mat, Vector<std::complex<double> > *rhs)
{
    mat->prealloc(n);
    for (std::map<unsigned int, MatrixEntry>::iterator it = ar_mat.begin(); it != ar_mat.end(); it++) {
      MatrixEntry &me = it->second;
      mat->pre_add_ij(me.m, me.n);
    }

    mat->alloc();
    for (std::map<unsigned int, MatrixEntry>::iterator it = ar_mat.begin(); it != ar_mat.end(); it++) {
      MatrixEntry &me = it->second;
      mat->add(me.m, me.n, me.value);
    }
    mat->finish();

    rhs->alloc(n);
    for (std::map<unsigned int, std::complex<double> >::iterator it = ar_rhs.begin(); it != ar_rhs.end(); it++) {
      rhs->add(it->first, it->second);
    }
    rhs->finish();
}

void build_matrix_block(int n, std::map<unsigned int, MatrixEntry> &ar_mat, std::map<unsigned int, std::complex<double> > &ar_rhs,
  SparseMatrix<std::complex<double> > *matrix, Vector<std::complex<double> > *rhs) {
    matrix->prealloc(n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        matrix->pre_add_ij(i, j);

    matrix->alloc();
    std::complex<double>  **mat = new_matrix<std::complex<double> >(n, n);
    int *cols = new int[n];
    int *rows = new int[n];
    for (int i = 0; i < n; i++) {
      cols[i] = i;
      rows[i] = i;
    }
    for (std::map<unsigned int, MatrixEntry>::iterator it = ar_mat.begin(); it != ar_mat.end(); it++) {
      MatrixEntry &me = it->second;
      mat[me.m][me.n] = me.value;
    }
    matrix->add(n, n, mat, rows, cols);
    matrix->finish();

    rhs->alloc(n);
    std::complex<double>  *rs = new std::complex<double>[n];
    for (std::map<unsigned int, std::complex<double> >::iterator it = ar_rhs.begin(); it != ar_rhs.end(); it++) {
      rs[it->first] = it->second;
    }
    unsigned int *u_rows = new unsigned int[n];
    for (int i = 0; i < n; i++)
      u_rows[i] = rows[i] >= 0 ? rows[i] : 0;
    rhs->add(n, u_rows, rs);
    rhs->finish();
}

// Test code.
void solve(LinearMatrixSolver<std::complex<double> > &solver, int n) {
  if (solver.solve()) {
    std::complex<double> *sln = solver.get_sln_vector();
    for (int i = 0; i < n; i++)
      if(sln[i].imag() < 0.0)
        std::cout << std::endl << sln[i].real() << sln[i].imag();
      else
        std::cout << std::endl << sln[i].real() << ' + ' << sln[i].imag();
  }
  else
    printf("Unable to solve.\n");
}

int main(int argc, char *argv[]) {
  int ret = TEST_SUCCESS;

  if (argc < 2) error("Not enough parameters.");

  int n;
  int nnz;
  bool cplx_2_real;

  std::map<unsigned int, MatrixEntry> ar_mat;
  std::map<unsigned int, std::complex<double> > ar_rhs;

  if (argc == 3 && strcasecmp(argv[2], "complex-matrix-to-real") == 0)
    cplx_2_real = true;
  else
    cplx_2_real = false;

  if (read_matrix_and_rhs((char*)"in/linsys-cplx-4", n, nnz, ar_mat, ar_rhs, cplx_2_real) != TEST_SUCCESS)
    error("Failed to read the matrix and rhs.");

  std::complex<double>* sln;

  if (strcasecmp(argv[1], "petsc") == 0) {
#ifdef WITH_PETSC
    PetscMatrix<std::complex<double> > mat;
    PetscVector<std::complex<double> > rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    PetscLinearMatrixSolver<std::complex<double> > solver(&mat, &rhs);
    solve(solver, n);
sln = solver.get_sln_vector();
#endif
  }
  else if (strcasecmp(argv[1], "petsc-block") == 0) {
#ifdef WITH_PETSC
    PetscMatrix<std::complex<double> > mat;
    PetscVector<std::complex<double> > rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    PetscLinearMatrixSolver<std::complex<double> > solver(&mat, &rhs);
    solve(solver, n);
sln = solver.get_sln_vector();
#endif
  }
  else if (strcasecmp(argv[1], "umfpack") == 0) {
#ifdef WITH_UMFPACK
    UMFPackMatrix<std::complex<double> > mat;
    UMFPackVector<std::complex<double> > rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    UMFPackLinearMatrixSolver<std::complex<double> > solver(&mat, &rhs);
    solve(solver, n);
sln = solver.get_sln_vector();
#endif
  }
  else if (strcasecmp(argv[1], "umfpack-block") == 0) {
#ifdef WITH_UMFPACK
    UMFPackMatrix<std::complex<double> > mat;
    UMFPackVector<std::complex<double> > rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    UMFPackLinearMatrixSolver<std::complex<double> > solver(&mat, &rhs);
    solve(solver, n);
sln = solver.get_sln_vector();
#endif
  }
  else if (strcasecmp(argv[1], "aztecoo") == 0) {
#ifdef WITH_TRILINOS
    EpetraMatrix<std::complex<double> > mat;
    EpetraVector<std::complex<double> > rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    AztecOOSolver<std::complex<double> > solver(&mat, &rhs);
    solve(solver, n);
sln = solver.get_sln_vector();
#endif
  }
  else if (strcasecmp(argv[1], "aztecoo-block") == 0) {
#ifdef WITH_TRILINOS
    EpetraMatrix<std::complex<double> > mat;
    EpetraVector<std::complex<double> > rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    AztecOOSolver<std::complex<double> > solver(&mat, &rhs);
    solve(solver, n);
sln = solver.get_sln_vector();
#endif
  }
  else if (strcasecmp(argv[1], "amesos") == 0) {
#ifdef WITH_TRILINOS
    EpetraMatrix<std::complex<double> > mat;
    EpetraVector<std::complex<double> > rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    if (AmesosSolver<std::complex<double> >::is_available("Klu")) {
      AmesosSolver<std::complex<double> > solver("Klu", &mat, &rhs);
      solve(solver, n);
sln = solver.get_sln_vector();
    }
#endif
  }
  else if (strcasecmp(argv[1], "amesos-block") == 0) {
#ifdef WITH_TRILINOS
    EpetraMatrix<std::complex<double> > mat;
    EpetraVector<std::complex<double> > rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    if (AmesosSolver<std::complex<double> >::is_available("Klu")) {
      AmesosSolver<std::complex<double> > solver("Klu", &mat, &rhs);
      solve(solver, n);
sln = solver.get_sln_vector();
    }
#endif
  }
  else if (strcasecmp(argv[1], "mumps") == 0) {
#ifdef WITH_MUMPS
    MumpsMatrix<std::complex<double> > mat;
    MumpsVector<std::complex<double> > rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    MumpsSolver<std::complex<double> > solver(&mat, &rhs);
    solve(solver, n);
sln = solver.get_sln_vector();
#endif
  }
  else if (strcasecmp(argv[1], "mumps-block") == 0) {
#ifdef WITH_MUMPS
    MumpsMatrix<std::complex<double> > mat;
    MumpsVector<std::complex<double> > rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    MumpsSolver<std::complex<double> > solver(&mat, &rhs);
    solve(solver, n);
sln = solver.get_sln_vector();
#endif
  }
  else
    ret = TEST_FAILURE;

std::cout << sln[0] << sln[1] << sln[2];

  if (std::abs(sln[0] - std::complex<double>(0.800000, -0.600000)) > 1E-6 || std::abs(sln[1] - std::complex<double>(0.470588, -0.882353)) > 1E-6 || std::abs(sln[2] - std::complex<double>(0.486486, -0.918919)) > 1E-6)
    ret = TEST_FAILURE;
  else
    ret = TEST_SUCCESS;

  // Test
  if (ret == TEST_FAILURE)
    printf("Failure!\n");
  else
    printf("Success!\n");
}
