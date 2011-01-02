#define HERMES_REPORT_ALL

#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

//creates a table of precalculated constitutive functions 
bool get_constitutive_tables(const char* tables_filename, int method)
{
  
  info("Creating tables of constitutive functions (complicated real exponent relations) ");
  info("This will take some time, but it is worthwhile. \n");

  // table values dimension
  int bound = int(-TABLE_LIMIT/TABLE_PRECISION)+1 ;
  
  //allocating arrays 
  K_TABLE = new double*[MATERIAL_COUNT] ;
  for (int i=0; i<MATERIAL_COUNT; i++) {
    K_TABLE[i] = new double[bound];
  }
  
  dKdh_TABLE = new double*[MATERIAL_COUNT] ;
  for (int i=0; i<MATERIAL_COUNT; i++) {
    dKdh_TABLE[i] = new double[bound];
  }
  
  dKdh_TABLE = new double*[MATERIAL_COUNT] ;
  for (int i=0; i<MATERIAL_COUNT; i++) {
    dKdh_TABLE[i] = new double[bound];
  }
  
  ddKdhh_TABLE = new double*[MATERIAL_COUNT] ;
  for (int i=0; i<MATERIAL_COUNT; i++) {
    ddKdhh_TABLE[i] = new double[bound];
  }

  C_TABLE = new double*[MATERIAL_COUNT] ;
  for (int i=0; i<MATERIAL_COUNT; i++) {
    C_TABLE[i] = new double[bound];
  }
  
  
  dCdh_TABLE = new double*[MATERIAL_COUNT] ;
  for (int i=0; i<MATERIAL_COUNT; i++) {
    dCdh_TABLE[i] = new double[bound];
  }
  
  
  // Calculate and save K(h).
  info("Calculating and saving K(h).");
  for (int j=0; j<MATERIAL_COUNT; j++) {
    for (int i=0; i< bound; i++) {
      K_TABLE[j][i] = K(-TABLE_PRECISION*i, j);
    }
  }
  // Calculate and save dKdh(h).
  info("Calculating and saving dKdh(h).");
  for (int j=0; j<MATERIAL_COUNT; j++) {
    for (int i=0; i< bound; i++) {
      dKdh_TABLE[j][i] = dKdh(-TABLE_PRECISION*i, j);
    }
  }
  // Calculate and save C(h).
  info("Calculating and saving C(h).");
  for (int j=0; j<MATERIAL_COUNT; j++) {
    for (int i=0; i< bound; i++) {
      C_TABLE[j][i] = C(-TABLE_PRECISION*i, j);
    }
  }
  // Calculate and save ddKdhh(h).
  info("Calculating and saving ddKdhh(h).");
  for (int j=0; j<MATERIAL_COUNT; j++) {
    for (int i=0; i< bound; i++) {
      ddKdhh_TABLE[j][i] = ddKdhh(-TABLE_PRECISION*i, j);
    }
  }
  // Calculate and save dCdh(h).
  info("Calculating and saving dCdh(h).");
  for (int j=0; j<MATERIAL_COUNT; j++) {
    for (int i=0; i< bound; i++) {
      dCdh_TABLE[j][i] = dCdh(-TABLE_PRECISION*i, j);
    }
  }
      

  return true;
}


//simple Gaussian elimnation for full matrices called from init_polynomials procedure
bool gem_full(double** A, double* b, double* X, int n){
  int i,j,k;
  
  double** aa;
  double dotproduct, tmp;
  aa = new double*[n];
  
  for (i=0; i<n; i++){
    aa[i] = new double[n+1];
  }
  
  for (i=0;i<n; i++){
    for (j=0; j<n; j++){
      aa[i][j] = A[i][j] ;
    }
    aa[i][n] = b[i];
  }

  for (j=0; j<(n-1); j++){
    for (i=j+1; i<n; i++){
    tmp = aa[i][j]/aa[j][j];

      for (k=0; k<(n+1); k++){
	aa[i][k] = aa[i][k] - tmp*aa[j][k] ;
      }
    }
  }
  

  for (i=n-1; i>-1; i--){
    dotproduct=0.0;
    for (j=i+1; j<n; j++){
      dotproduct = dotproduct + aa[i][j]*X[j] ;
    }
    X[i] = (aa[i][n]-dotproduct)/aa[i][i] ;
  }
    
  
  delete []aa;
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








 





