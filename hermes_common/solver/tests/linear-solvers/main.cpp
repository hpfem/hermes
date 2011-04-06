#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE

//#include <getopt.h>
#include "common.h"
#include "config.h"

#include "solver/solver.h"
#include "solver/umfpack_solver.h"
#include "solver/superlu.h"
#include "solver/petsc.h"
#include "solver/epetra.h"
#include "solver/amesos.h"
#include "solver/aztecoo.h"
#include "solver/mumps.h"

#include <iostream>

// Test of linear solvers.
// Read matrix and RHS from a file.

// Max row length in input file.
#define MAX_ROW_LEN	1024

struct MatrixEntry {
  MatrixEntry() { }
  MatrixEntry(int m, int n, scalar value) {
  this->m = m;
  this->n = n;
  this->value = value;
  }

  int m, n;
  scalar value;
};

void show_mat(const char *msg, std::map<unsigned int, MatrixEntry> mp) {
  std::map<unsigned int, MatrixEntry>::iterator itr;

  std::cout << msg << std::endl;

  for(itr=mp.begin(); itr != mp.end(); ++itr)
    std::cout << " " << (int) itr->first << ": " << 
   (int) itr->second.m << " " << 
   (int) itr->second.n << " " << 
   (scalar) itr->second.value << 
   std::endl;

  std::cout << std::endl;
}

void show_rhs(const char *msg, std::map<unsigned int,scalar> mp) {
  std::map<unsigned int, scalar>::iterator itr;

  std::cout << msg << std::endl;

  for(itr=mp.begin(); itr != mp.end(); ++itr)
   std::cout << " " << (int) itr->first << ": " << (scalar) itr->second << std::endl;

  std::cout << std::endl;
}

// Helpers.

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
                        std::map<unsigned int, MatrixEntry> &mat, std::map<unsigned int, scalar> &rhs, bool &cplx_2_real) 
{
  FILE *file = fopen(file_name, "r");
  if (file == NULL) return ERR_FAILURE;

  enum EState {
    STATE_N,
    STATE_NNZ,
    STATE_MATRIX,
    STATE_RHS,
  } state = STATE_N;

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
          if (cplx_2_real){
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
        if (read_n_nums(row, 1, buffer)) {
             nnz = (int) buffer[0];
          } 

        state = STATE_MATRIX; 
      break;

#ifndef HERMES_COMMON_COMPLEX

      case STATE_MATRIX:
        if (cplx_2_real){
          if (read_n_nums(row, 4, buffer)) { 

            mat.insert(std::pair<unsigned int, MatrixEntry>(k,   MatrixEntry((int) buffer[0],     (int) buffer[1],     buffer[2])));

            mat.insert(std::pair<unsigned int, MatrixEntry>(k+1, MatrixEntry((int) buffer[0] + n/2, (int) buffer[1],     buffer[3])));

            mat.insert(std::pair<unsigned int, MatrixEntry>(k+2*nnz, MatrixEntry((int) buffer[0],     (int) buffer[1] + n/2, (-1)*buffer[3])));

            mat.insert(std::pair<unsigned int, MatrixEntry>(k+2*nnz+1, MatrixEntry((int) buffer[0] + n/2, (int) buffer[1] + n/2, buffer[2])));

            k=k+2;
          }
           else        
           state = STATE_RHS;
          
     
        }else{ // if cplx_2_real is false.
           if (read_n_nums(row, 3, buffer)) 
             mat[mat.size()] = (MatrixEntry((int) buffer[0], (int) buffer[1], buffer[2]));

           else        
           state = STATE_RHS;
         }
      break; //case STATE_MATRIX break.

      case STATE_RHS:
        if (cplx_2_real) {
          if (read_n_nums(row, 3, buffer)) {

              if (buffer[0] != (int) n/2-1) // Then this is not the last line in the input file
              {
                rhs[((int) buffer[0])] = (scalar) buffer[1];
                rhs_buffer[l] = (scalar) buffer[2];
                l=l+1;
              }
                         
              else // This is last line in the file.
              {
                // First read last line entry
                rhs[((int) buffer[0])] = (scalar) buffer[1];
                rhs_buffer[l] = (scalar) buffer[2];
                l=l+1;
                // Take imaginary parts you saved, 
                // and fill the rest of the rhs vector.
                for (int i=0; i < l; i++) 
                  {
                    rhs[rhs.size()] = rhs_buffer[i];                  
                  }
              }
           }
        }
        else { // if cplx_2_real is false.
          if (read_n_nums(row, 2, buffer)) 
            rhs[(int) buffer[0]] = (scalar) buffer[1];
        }                 
        break;
    }
  }

#else

      case STATE_MATRIX:
        if (read_n_nums(row, 4, buffer)) {
            std::complex<double> cmplx_buffer(buffer[2], buffer[3]);
          mat[mat.size()] = (MatrixEntry((int) buffer[0], (int) buffer[1], (scalar) cmplx_buffer));
        }
	else
        state = STATE_RHS;
      break;

      case STATE_RHS:
        if (read_n_nums(row, 3, buffer)) {
        std::complex<double> cmplx_buffer(buffer[1], buffer[2]);
          rhs[(int) buffer[0]] = (scalar) cmplx_buffer;
        }
      break;
    }
  }

#endif

  fclose(file);

  // Free memory
  delete [] rhs_buffer;
  
  //Clear pointers
  rhs_buffer = NULL;

  return ERR_SUCCESS;
}

void build_matrix(int n, std::map<unsigned int, MatrixEntry> &ar_mat, 
                  std::map<unsigned int, scalar> &ar_rhs, SparseMatrix *mat,
                  Vector *rhs) {
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
  for (std::map<unsigned int, scalar>::iterator it = ar_rhs.begin(); it != ar_rhs.end(); it++) {
    rhs->add(it->first, it->second);
  }
  rhs->finish();
}

void build_matrix_block(int n, std::map<unsigned int, MatrixEntry> &ar_mat, std::map<unsigned int, scalar> &ar_rhs,
                        SparseMatrix *matrix, Vector *rhs) {
  matrix->prealloc(n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      matrix->pre_add_ij(i, j);

  matrix->alloc();
  scalar **mat = new_matrix<scalar>(n, n);
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
  scalar *rs = new scalar[n];
  for (std::map<unsigned int, scalar>::iterator it = ar_rhs.begin(); it != ar_rhs.end(); it++) {
    rs[it->first] = it->second;
  }
  unsigned int *u_rows = new unsigned int[n];
  for (int i = 0; i < n; i++)
    u_rows[i] = rows[i] >= 0 ? rows[i] : 0;
  rhs->add(n, u_rows, rs);
  rhs->finish();
}

// Test code.
void solve(Solver &solver, int n) {
  if (solver.solve()) {
    scalar *sln = solver.get_solution();
    for (int i = 0; i < n; i++) {
      printf(SCALAR_FMT"\n", SCALAR(sln[i]));
    }
  }
  else {
    printf("Unable to solve.\n");
  }
}

int main(int argc, char *argv[]) {
  int ret = ERR_SUCCESS;

#ifndef HERMES_COMMON_COMPLEX
  if (argc < 3) error("Not enough parameters.");
#else
  if (argc < 3) error("Not enough parameters.");
#endif

  int n;
  int nnz;                            
  bool cplx_2_real;
  
  std::map<unsigned int, MatrixEntry> ar_mat;
  std::map<unsigned int, scalar> ar_rhs;

  if (argc == 4 && strcasecmp(argv[3],"complex-matrix-to-real") == 0)
     cplx_2_real = true;
  else
     cplx_2_real = false;

  if (read_matrix_and_rhs(argv[2], n, nnz, ar_mat, ar_rhs, cplx_2_real) != ERR_SUCCESS)
    error("Failed to read the matrix and rhs.");

  //show_mat("Here is the original ar_mat: ", ar_mat);

  //show_rhs("Here is the original ar_rhs: ", ar_rhs);

  if (strcasecmp(argv[1], "petsc") == 0) {
#ifdef WITH_PETSC
    PetscMatrix mat;
    PetscVector rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    PetscLinearSolver solver(&mat, &rhs);
    solve(solver, n);
#endif
  }
  else if (strcasecmp(argv[1], "petsc-block") == 0) {
#ifdef WITH_PETSC
    PetscMatrix mat;
    PetscVector rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    PetscLinearSolver solver(&mat, &rhs);
    solve(solver, n);
#endif
  }
  else if (strcasecmp(argv[1], "umfpack") == 0) {
#ifdef WITH_UMFPACK
    UMFPackMatrix mat;
    UMFPackVector rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    UMFPackLinearSolver solver(&mat, &rhs);
    solve(solver, n);
#endif
  }
  else if (strcasecmp(argv[1], "umfpack-block") == 0) {
#ifdef WITH_UMFPACK
    UMFPackMatrix mat;
    UMFPackVector rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    UMFPackLinearSolver solver(&mat, &rhs);
    solve(solver, n);
#endif
  }
  else if (strcasecmp(argv[1], "aztecoo") == 0) {
#ifdef WITH_TRILINOS
    EpetraMatrix mat;
    EpetraVector rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    AztecOOSolver solver(&mat, &rhs);
    solve(solver, n);
#endif
  }
  else if (strcasecmp(argv[1], "aztecoo-block") == 0) {
#ifdef WITH_TRILINOS
    EpetraMatrix mat;
    EpetraVector rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    AztecOOSolver solver(&mat, &rhs);
    solve(solver, n);
#endif
  }
  else if (strcasecmp(argv[1], "amesos") == 0) {
#ifdef WITH_TRILINOS
    EpetraMatrix mat;
    EpetraVector rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    if (AmesosSolver::is_available("Klu")) {
      AmesosSolver solver("Klu", &mat, &rhs);
      solve(solver, n);
    }
#endif
  }
  else if (strcasecmp(argv[1], "amesos-block") == 0) {
#ifdef WITH_TRILINOS
    EpetraMatrix mat;
    EpetraVector rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    if (AmesosSolver::is_available("Klu")) {
      AmesosSolver solver("Klu", &mat, &rhs);
      solve(solver, n);
    } 
#endif
  }
  else if (strcasecmp(argv[1], "mumps") == 0) {
#ifdef WITH_MUMPS
    MumpsMatrix mat;
    MumpsVector rhs;
    build_matrix(n, ar_mat, ar_rhs, &mat, &rhs);

    MumpsSolver solver(&mat, &rhs);
    solve(solver, n);
#endif
  }
  else if (strcasecmp(argv[1], "mumps-block") == 0) {
#ifdef WITH_MUMPS
    MumpsMatrix mat;
    MumpsVector rhs;
    build_matrix_block(n, ar_mat, ar_rhs, &mat, &rhs);

    MumpsSolver solver(&mat, &rhs);
    solve(solver, n);
#endif
  }  
  else
    ret = ERR_FAILURE;

  return ret;
}
