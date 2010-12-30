#define HERMES_REPORT_ALL

bool solve_picard(WeakForm* wf, Space* space, Solution* sln_prev_iter,
                 MatrixSolverType matrix_solver, double PICARD_TOL, int MAX_PICARD_ITER_NUM, bool verbose) 
{
  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  Solution sln_new;

  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(wf, space, is_linear);

  int i_err;
  int iter_count = 0;
  while (true) {
    // Assemble the stiffness matrix and right-hand side.
    info("Assembling the stiffness matrix and right-hand side vector.");
    dp.assemble(matrix, rhs);

    // Solve the linear system and if successful, obtain the solution.
    info("Solving the matrix problem.");
    if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), space, &sln_new);
    else error ("Matrix solver failed.\n");

    double rel_error = calc_abs_error(sln_prev_iter, &sln_new, HERMES_H1_NORM) 
                       / calc_norm(&sln_new, HERMES_H1_NORM) * 100;
    info("Relative error: %g%%", rel_error);

    // Stopping criterion.
    if (rel_error < PICARD_TOL) {
      sln_prev_iter->copy(&sln_new);
      delete matrix;
      delete rhs;
      delete solver;
      return true;
    }
    
    if (iter_count >= MAX_PICARD_ITER_NUM) {
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

// this procedure should be moved into constitutive_genuchten.cpp
double precalculate_constitutive_values()
{
  info("Precalculating constitutive relations.");
  info("This will take some time, but it is worthwhile.");
  
  info("Precalculating K(h).");
  for (int j=0; j<4; j++) {
    for (int i=0; i< 1500000; i++) {
      K_TABLE[j][i] = K(-0.01*i, j);
    }
  }
  
  info("Precalculating dKdh(h).");
  for (int j=0; j<4; j++) {
    for (int i=0; i< 1500000; i++) {
      dKdh_TABLE[j][i] = dKdh(-0.01*i, j);
    }
  }
  
  info("Precalculating ddKdhh(h).");
  for (int j=0; j<4; j++) {
    for (int i=0; i< 1500000; i++) {
      ddKdhh_TABLE[j][i] = ddKdhh(-0.01*i, j);
    }
  }
  
  info("Precalculating function C(h).");
  for (int j=0; j<4; j++) {
    for (int i=0; i< 1500000; i++) {
      C_TABLE[j][i] = C(-0.01*i, j);
    }
  }
   
  info("Precalculating function dCdh(h).");
  for (int j=0; j<4; j++) {
    for (int i=0; i< 1500000; i++) {
      dCdh_TABLE[j][i] = dCdh(-0.01*i, j);
    }
  }
   
  CONSTITUTIVE_TABLES_READY = true;
  return 0;
}

//initialize polynomial approximation of the strange LOCO functions close to a full saturation
// n - degree of polynomials
// low_limit - start point of the polynomial approximation
// points - array of points inside the interval bounded by <low_limit, 0> to improve the accuracy, at least one is recommended. 
//An approriate amount of points related to the polynomial degree should be selected.
//n_inside_point - number of inside points
// layer - material to be considered
int init_polynomials(int n, double low_limit, double *points, int n_inside_points, int layer){
  double** Aside;
  double* Bside;
  double* X;
  
  X = new double[7];
  
  if (POLYNOMIALS_ALLOCATED < 0){

    POLYNOMIALS = new double**[MATERIAL_COUNT];
    
    for (int i=0; i<MATERIAL_COUNT; i++){
      POLYNOMIALS[i] = new double*[3] ;
    }
    
    for (int i=0; i<MATERIAL_COUNT; i++){
      for (int j=0; j<3; j++){
	POLYNOMIALS[i][j] = new double[n_inside_points+6] ;
      }
    }
    POLYNOMIALS_ALLOCATED = 1;
  }
  

  
  Aside = new double*[n] ;
  Bside = new double[n] ;
  for (int i=0; i<n; i++){
    Aside[i] = new double[n] ;
  }

  //evaluate the first three rows of the matrix (zero, first and second derivation at point low_limit)
  
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
  
  //evaluate the second three rows of the matrix (zero, first and second derivation at point zero)


  
 Aside[0][0] = 1.0 ;
 //for the both first and second derivation it does not really matter what value is placed there.
 Aside[1][1] = 1.0 ;
 Aside[2][2] = 2.0 ;
 
 Bside[0] = K(0.0, layer) ;
 Bside[1] = 0.0;
 Bside[2] = 0.0;
 

 
 for (int i=6; i<(6+n_inside_points); i++){
   for (int j=0; j<n; j++){
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


 
  delete [] Aside ;
  delete [] Bside ;
  delete [] X;
  
  return 0;
}








 





