#define HERMES_REPORT_ALL

#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

bool solve_picard(WeakForm* wf, Space* space, Solution* sln_prev_iter,
                  MatrixSolverType matrix_solver, double picard_tol, 
                  int picard_max_iter, bool verbose) 
{
  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  Solution sln_new;

  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(wf, space, is_linear);

  int iter_count = 0;
  while (true) {
    // Assemble the stiffness matrix and right-hand side.
    dp.assemble(matrix, rhs);

    // Solve the linear system and if successful, obtain the solution.
    if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), space, &sln_new);
    else error ("Matrix solver failed.\n");

    double rel_error = calc_abs_error(sln_prev_iter, &sln_new, HERMES_H1_NORM) 
                       / calc_norm(&sln_new, HERMES_H1_NORM) * 100;
    if (verbose) info("---- Picard iter %d, ndof %d, rel. error %g%%", 
                 iter_count+1, Space::get_num_dofs(space), rel_error);

    // Stopping criterion.
    if (rel_error < picard_tol) {
      sln_prev_iter->copy(&sln_new);
      delete matrix;
      delete rhs;
      delete solver;
      return true;
    }
    
    if (iter_count >= picard_max_iter) {
      delete matrix;
      delete rhs;
      delete solver;
      return false;
    }
    
    // Saving solution for the next iteration;
    sln_prev_iter->copy(&sln_new);
   
    iter_count++;
  }
}

// Look for a file with precalculated constitutive relations. 
// If found, read it. If not, create it. 
bool get_constitutive_tables(const char* tables_filename, int method)
{
  bool file_found;
  ifstream inFile;

  // NOTE: The order of tables is important!
  inFile.open(tables_filename);
  if (inFile) {
    info("File with precalculated constitutive tables found, reading.");
    // Reading K.
    info("Reading K(h).");
    for (int j=0; j<4; j++) {
      for (int i=0; i< 1500000; i++) {
        inFile >> K_TABLE[j][i];
      }
    }
    // Reading dKdh.
    info("Reading dKdh(h).");
    for (int j=0; j<4; j++) {
      for (int i=0; i< 1500000; i++) {
        inFile >> dKdh_TABLE[j][i];
      }
    }
    // Reading C(h).
    info("Reading C(h).");
    for (int j=0; j<4; j++) {
      for (int i=0; i< 1500000; i++) {
        inFile >> C_TABLE[j][i];
      }
    }
    // If Picard, stop here.
    if (method != 1) {
      inFile.close();    
      return true;
    }
    // Reading ddKdhh.
    info("Reading ddKdhh(h).");
    for (int j=0; j<4; j++) {
      for (int i=0; i< 1500000; i++) {
        inFile >> ddKdhh_TABLE[j][i];
      }
    }
    // Reading dCdh(h).
    info("Reading dCdh(h).");
    for (int j=0; j<4; j++) {
      for (int i=0; i< 1500000; i++) {
        inFile >> dCdh_TABLE[j][i];
      }
    }
    // Closing input file and returning.
    inFile.close();
    return true;
  }

  /*** TABLES DATA FILE WAS NOT FOUND, CREATING ***/

  info("File %s with precalculated constitutive tables not found, creating.", tables_filename);
  info("This will take some time, but it is worthwhile.");
  FILE* f = fopen(tables_filename, "w");
  if (f == NULL) error("Could not open file %s for writing.", tables_filename);
  // Calculate and save K(h).
  info("Calculating and saving K(h).");
  for (int j=0; j<4; j++) {
    for (int i=0; i< 1500000; i++) {
      fprintf(f, "%15.10f ", K_TABLE[j][i] = K(-0.01*i, j));
    }
  }
  // Calculate and save dKdh(h).
  info("Calculating and saving dKdh(h).");
  for (int j=0; j<4; j++) {
    for (int i=0; i< 1500000; i++) {
      fprintf(f, "%15.10f ", dKdh_TABLE[j][i] = dKdh(-0.01*i, j));
    }
  }
  // Calculate and save C(h).
  info("Calculating and saving C(h).");
  for (int j=0; j<4; j++) {
    for (int i=0; i< 1500000; i++) {
      fprintf(f, "%15.10f ", C_TABLE[j][i] = C(-0.01*i, j));
    }
  }
  // Calculate and save ddKdhh(h).
  info("Calculating and saving ddKdhh(h).");
  for (int j=0; j<4; j++) {
    for (int i=0; i< 1500000; i++) {
      fprintf(f, "%15.10f ", ddKdhh_TABLE[j][i] = ddKdhh(-0.01*i, j));
    }
  }
  // Calculate and save dCdh(h).
  info("Calculating and saving dCdh(h).");
  for (int j=0; j<4; j++) {
    for (int i=0; i< 1500000; i++) {
      fprintf(f, "%15.10f ", dCdh_TABLE[j][i] = dCdh(-0.01*i, j));
    }
  }
      
  fclose(f);
  return true;
}

// Initialize polynomial approximation of constitutive relations close to full saturation.
// n - degree of polynomials
// low_limit - start point of the polynomial approximation
// points - array of points inside the interval bounded by <low_limit, 0> to improve the accuracy, at least one is recommended. 
// An approriate amount of points related to the polynomial degree should be selected.
// n_inside_point - number of inside points
// layer - material to be considered
int init_polynomials(int n, double low_limit, double *points, int n_inside_points, int layer){
  double** Aside;
  double* Bside;
  double* X;
  
  X = new double[7];
  
  if (POLYNOMIALS_ALLOCATED == false){

    POLYNOMIALS = new double**[MATERIAL_COUNT];
    
    for (int i=0; i<MATERIAL_COUNT; i++){
      POLYNOMIALS[i] = new double*[3] ;
    }
    
    for (int i=0; i<MATERIAL_COUNT; i++){
      for (int j=0; j<3; j++){
	POLYNOMIALS[i][j] = new double[n_inside_points+6] ;
      }
    }
    POLYNOMIALS_ALLOCATED = true;
  }

  Aside = new double*[n] ;
  Bside = new double[n] ;
  for (int i=0; i<n; i++){
    Aside[i] = new double[n] ;
  }

  // evaluate the first three rows of the matrix (zero, first and second derivative at point low_limit)
  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++){
      Aside[i][j] = 0.0;
    }
  }
  for (int i=0; i<n; i++){
    Aside[3][i] = pow(low_limit,i) ;
    Aside[4][i] = i*pow(low_limit, i-1) ;
    Aside[5][i] = i*(i-1)*pow(low_limit, i-2) ;
  }
  Bside[3] = K(low_limit, layer) ;
  Bside[4] = dKdh(low_limit, layer) ;
  Bside[5] = ddKdhh(low_limit, layer) ; 
  
  //evaluate the second three rows of the matrix (zero, first and second derivative at point zero)
  Aside[0][0] = 1.0 ;

  //for the both first and second derivative it does not really matter what value is placed there.
  Aside[1][1] = 1.0 ;
  Aside[2][2] = 2.0 ;
 
  Bside[0] = K(0.0, layer) ;
  Bside[1] = 0.0;
  Bside[2] = 0.0;
 
  for (int i=6; i<(6+n_inside_points); i++){
    for (int j=0; j<n; j++) {
      Aside[i][j]=pow(points[i-6],j) ;
    }
    printf("poradi, %i %lf %lf \n", i, K(points[i-6], layer), points[i-6]);

    Bside[i] = K(points[i-6], layer) ;
    printf("layer, %i \n", layer);
  }

  gem_full(Aside, Bside, POLYNOMIALS[layer][0], (n_inside_points+6));
  
  for (int i=1; i<3; i++){
    for (int j=0; j< (n_inside_points+5); j++){
      POLYNOMIALS[layer][i][j] = (j+1)*POLYNOMIALS[layer][i-1][j+1];
    }
    POLYNOMIALS[layer][i][n_inside_points+6-i] = 0.0;
  }

  delete [] Aside;
  delete [] Bside;
  delete [] X;
  
  return 0;
}








 





