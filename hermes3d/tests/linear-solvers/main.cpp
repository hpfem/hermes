#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "config.h"
//#include <getopt.h>
#include <hermes3d.h>

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

int read_matrix_and_rhs(char *file_name, int &n, 
                        std::map<unsigned int, MatrixEntry> &mat, std::map<unsigned int, scalar> &rhs, bool &cplx_2_real) 
{
  FILE *file = fopen(file_name, "r");
  if (file == NULL) return ERR_FAILURE;

  enum EState {
    STATE_N,
    STATE_MATRIX,
    STATE_RHS,
  } state = STATE_N;

  // Variables needed to turn complex matrix into real.
  int k = 0; 
  int l = 0;
  int param = 0;
  int im_i_coord[n];
  int im_j_coord[n];
  int i_coord[n];
  double im_buffer[n];
  double rhs_buffer[n];
  bool fill_upper = true;
  bool fill_lower = false;
  fpos_t pos; //Object containing information to specify a position within a file.

  double buffer[4]; 
  char row[MAX_ROW_LEN];
  while (fgets(row, MAX_ROW_LEN, file) != NULL) {
    switch (state) {
      case STATE_N:
        if (read_n_nums(row, 1, buffer)) {
          if (cplx_2_real) {
             if (fill_upper==true) { //sumnjivo
             n = 2*((int) buffer[0]); 
             printf("%d\n",n);
             }//sumnjivo
             else 
               state = STATE_MATRIX;
          }   
          else{ 
             n = (int) buffer[0];
          } 
         //fgetpos (file, &pos); // Get current position in stream.                               
          state = STATE_MATRIX;
        } 
      break;

#ifndef H3D_COMPLEX

      case STATE_MATRIX:
        if (cplx_2_real){
          if (read_n_nums(row, 4, buffer)) {

              if (fill_upper == true) { //case FILL_UPPER:
           
                if (buffer[0] == k)
                {
                   mat[mat.size()] = (MatrixEntry((int) buffer[0], (int) buffer[1], buffer[2]));
                   std::cout << "Re upper " << mat.size() << std::endl;

                   im_buffer[l] = (-1)*buffer[3];
                   im_i_coord[l] = (int) buffer[0];
                   im_j_coord[l] = (int) buffer[1] + n; 
                   l=l+1;
                }
                else // we have read an element that belongs to next row. Now what?
                {
                   for (int i=0; i < l; i++) //OK, first finish the previous row with imaginaries.
                   {
                     mat[mat.size()] = (MatrixEntry(im_i_coord [i], im_j_coord [i], im_buffer[i]));
                     std::cout << "Fill the row with Im " << mat.size() << std::endl;
                   }

                   // Now begin the next row with real
                   mat[mat.size()] = (MatrixEntry((int) buffer[0], (int) buffer[1], buffer[2]));
                   std::cout <<"Next row with Re "<< mat.size() << std::endl;
                   // And remember imaginary part from this row
                   im_buffer[0] = (-1)*buffer[3];
                   im_i_coord[0] = (int) buffer[0];
                   im_j_coord[0] = (int) buffer[1] + n; 
                   l = 1;
 
                   while(k < buffer[0]) // Fast forward k, we need it for next read.
                     k++;
                }
              } // break;  // case FILL_UPPER break.

              if (fill_lower == true){ //case FILL_LOWER: // we fill in lower half of the matrix, by different rules
 
                if (buffer[0] == k)
                {
                   // Create next matrix entry
                   mat[mat.size()] = (MatrixEntry((int) buffer[0], (int) buffer[1], buffer[3]));
                   std::cout <<"Re lower " << mat.size() << std::endl;

                   // Save imaginary part for later
                   im_buffer[l] = buffer[2];
                   im_i_coord[l] = (int) buffer[0];
                   im_j_coord[l] = (int) buffer[1] + n; 
                   l=l+1;                 

                }
                else // we have read an element that belongs to next row. Now what?
                {
                   //OK, first finish the previous row with imaginaries.
                   for (int i=0; i < l; i++) 
                   {
                     mat[mat.size()] = (MatrixEntry(im_i_coord [i], im_j_coord [i], im_buffer[i]));
                     std::cout << " Finish the row with Im lower " << mat.size() << std::endl;
                   }

                   // Now begin the next row with real
                   mat[mat.size()] = (MatrixEntry((int) buffer[0], (int) buffer[1], buffer[3]));
                   std::cout <<" Next row lower " << mat.size() << std::endl;
                   // Save imaginary part for later
                   im_buffer[0] = buffer[2];
                   im_i_coord[0] = (int) buffer[0];
                   im_j_coord[0] = (int) buffer[1] + n; 
                   l = 1;

                   // Fast forward k, we need it for next read. 
                   while(k < buffer[0]) 
                     k++;
                }
                } //break; // case FILL_LOWER break.

 //----           } // switch state: FILL_UPPER, FILL_LOWER

          } // if read_n_nums block end.

 	  else {
           // We have reached the line in the input file, 
           // where the first rhs vector entries are defined. 
           // But before we go to STATE_RHS, 
           // let's finish filling the lower half of the new real matrix. 
           // Now we are sure we have all the elements we need.
              if (param == 0){
                rewind (file); //fsetpos (file, &pos); //sumnjivo
                state = STATE_N; // sumnjivo
                fill_lower = true;
                fill_upper = false;
                param = 1;
              }
              else {
                l = 0;
                state = STATE_RHS;
              }
          }      
     
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

              if (buffer[0] != n-1) // Then this is not the last line in the file
              {
                rhs[((int) buffer[0])] = buffer[1];
                rhs_buffer[l] = buffer[2];
                i_coord[l] = (int) buffer[0] + n;
                l++;
              }
                         
              else //if (buffer[0] == n-1). This is last line in the file.
              {
                // Take imaginary parts you saved, 
                // and fill the rest of the rhs vector.
                for (int i=0; i < l; i++) 
                  {
                    rhs[i_coord[i]] = rhs_buffer[i];                  
                  }
              }

//            rhs[((int) buffer[0])] = buffer[1];
//            rhs[((int) buffer[0]) + n] = buffer[2];
//             rhs.insert(pair<unsigned int, scalar>((int) buffer[0],     buffer[1]));
//             rhs.insert(pair<unsigned int, scalar>((int) buffer[0] + n, buffer[2]));
//              printf("Done both rhs!\n");
          }
        }
        else { // if cplx_2_real is false.
          if (read_n_nums(row, 2, buffer)) 
            rhs[(int) buffer[0]] = buffer[1];
        }                 
        break;
    }
  }

#else

      case STATE_MATRIX:
        if (read_n_nums(row, 4, buffer)) {
          complex<double> cmplx_buffer(buffer[2], buffer[3]);
          mat[mat.size()] = (MatrixEntry((int) buffer[0], (int) buffer[1], (scalar) cmplx_buffer));
        }
	else
        state = STATE_RHS;
      break;

      case STATE_RHS:
        if (read_n_nums(row, 3, buffer)) {
        complex<double> cmplx_buffer(buffer[1], buffer[2]);
          rhs[(int) buffer[0]] = (scalar) cmplx_buffer;
        }
      break;
    }
  }

#endif

  fclose(file);

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
  printf("I'm here\n");
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

#ifndef H3D_COMPLEX
  if (argc < 3) error("Not enough parameters.");
#else
  if (argc < 3) error("Not enough parameters.");
#endif

  int n;
  bool cplx_2_real;
  
  std::map<unsigned int, MatrixEntry> ar_mat;
  std::map<unsigned int, scalar> ar_rhs;

  if (argc == 4 && strcasecmp(argv[3],"complex-matrix-to-real") == 0)
     cplx_2_real = true;
  else
     cplx_2_real = false;

  if (read_matrix_and_rhs(argv[2], n, ar_mat, ar_rhs, cplx_2_real) != ERR_SUCCESS)
    error("Failed to read the matrix and rhs.");


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
