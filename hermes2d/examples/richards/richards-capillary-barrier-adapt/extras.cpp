#define HERMES_REPORT_ALL

#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;


//Debugging matrix printer.
bool printmatrix(double** A, int n, int m){
  for (int i=0; i<n; i++){
    for (int j=0; j<m; j++){
      printf(" %lf ", A[i][j]) ;
    }
    printf(" \n");
  }
  printf("----------------------------------\n");
  return true;
}


//Debugging vector printer.
bool printvector(double* vect, int n){
  for (int i=0; i<n; i++){
    printf(" %lf ", vect[i]);
  }
  printf("\n");
  printf("----------------------------------\n");
  return true;
}

// Creates a table of precalculated constitutive functions.
bool get_constitutive_tables(int method)
{
  info("Creating tables of constitutive functions (complicated real exponent relations).");

  // Table values dimension.
  int bound = int(-TABLE_LIMIT/TABLE_PRECISION)+1;
  
  // Allocating arrays. 
  K_TABLE = new double*[MATERIAL_COUNT];
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
 

  C_TABLE = new double*[MATERIAL_COUNT] ;
  for (int i=0; i<MATERIAL_COUNT; i++) {
    C_TABLE[i] = new double[bound];
  }
  
  //If Newton method (method==1) selected constitutive function derivations are required.
  if (method==1){
    dCdh_TABLE = new double*[MATERIAL_COUNT] ;
    for (int i=0; i<MATERIAL_COUNT; i++) {
      dCdh_TABLE[i] = new double[bound];
    }
    
    ddKdhh_TABLE = new double*[MATERIAL_COUNT] ;
    for (int i=0; i<MATERIAL_COUNT; i++) {
      ddKdhh_TABLE[i] = new double[bound];
    }
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
  
  
  //If Newton method (method==1) selected constitutive function derivations are required.
  if (method==1){
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
  }	
      
  return true;
}

// Simple Gaussian elimination for full matrices called from init_polynomials().
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

// Initialize polynomial approximation of constitutive relations close to full saturation for CONSTITUTIVE_TABLE_METHOD=1.
// For CONSTITUTIVE_TABLE_METHOD=2 all constitutive functions are approximated by polynomials, K(h) function by quintic spline, C(h) 
// function by cubic splines. Discretization is managed by variable int NUM_OF_INTERVALS and double* INTERVALS_4_APPROX.
// ------------------------------------------------------------------------------
// For CONSTITUTIVE_TABLE_METHOD=1 this function requires folowing arguments:
// n - degree of polynomials
// low_limit - start point of the polynomial approximation
// points - array of points inside the interval bounded by <low_limit, 0> to improve the accuracy, at least one is recommended. 
// An approriate amount of points related to the polynomial degree should be selected.
// n_inside_point - number of inside points
// layer - material to be considered.
//------------------------------------------------------------------------------
// For CONSTITUTIVE_TABLE_METHOD=2, all parameters are obtained from global definitions.
bool init_polynomials(int n, double low_limit, double *points, int n_inside_points, int layer){
  double** Aside;
  double* Bside;
  double* X;
  switch (CONSTITUTIVE_TABLE_METHOD) 
    { 
      // no approximation 
      case 0 :
	break ;
      // polynomial approximation only for the the K(h) function surroundings close to zero
      case 1 :
	
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

	Aside = new double*[n+n_inside_points] ;
	Bside = new double[n+n_inside_points] ;
	for (int i=0; i<n; i++){
	  Aside[i] = new double[n+n_inside_points] ;
	}

	// Evaluate the first three rows of the matrix (zero, first and second derivative at point low_limit).
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
	
	// Evaluate the second three rows of the matrix (zero, first and second derivative at point zero).
	Aside[0][0] = 1.0 ;

	// For the both first and second derivative it does not really matter what value is placed there.
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
	break ;
      // polynomial approximation for all functions at interval (TABLE_LIMIT, 0)
      case 2 :
	
	int pts = 0;
	
	if (POLYNOMIALS_ALLOCATED == false) {
	  // K(h) function is approximated by quintic spline.
	  K_POLS = new double***[NUM_OF_INTERVALS];
	  //C(h) function is approximated by cubic spline.
	  C_POLS = new double***[NUM_OF_INTERVALS];
	  
	  for (int i=0; i<NUM_OF_INTERVALS; i++){
	    K_POLS[i] = new double**[MATERIAL_COUNT];
	    C_POLS[i] = new double**[MATERIAL_COUNT];
	    for (int j=0; j<MATERIAL_COUNT; j++){
	      K_POLS[i][j] = new double*[3];
	      C_POLS[i][j] = new double*[2];
	      for (int k=0; k<3; k++) {
		K_POLS[i][j][k] = new double[6 + pts];
		if (k<2)
		  C_POLS[i][j][k] = new double[4];
	      }
	    }
	  }
	  
// 	  //allocate POL_SEARCH_HELP array -- an index array with locations for particular pressure head functions
	  POL_SEARCH_HELP = new int[int(-TABLE_LIMIT)+1];
	  
	  for (int i=0; i<int(-TABLE_LIMIT); i++){
	    for (int j=0; j<NUM_OF_INTERVALS; j++){
	      if (j < 1) {
		if (-i > INTERVALS_4_APPROX[j]) {
		  POL_SEARCH_HELP[i] = j ;
		  break ;
		}
	      }
	      else {
		if (-i > INTERVALS_4_APPROX[j] && -i <= INTERVALS_4_APPROX[j-1]) {
		  POL_SEARCH_HELP[i] = j ;
		  break ;
		 }
	       }
	     }
	   }
	  
	  POLYNOMIALS_ALLOCATED = true;
	}
	//create matrix
	Aside = new double*[6 + pts];
	for (int i=0; i< (6 + pts); i++){
	  Aside[i] = new double[6+pts];
	}
	
	Bside = new double[6 + pts];


        for (int i=0; i<NUM_OF_INTERVALS; i++) {

	  if (i < 1){
	    for (int j=0; j<3; j++){
	      for (int k=0; k<(6 + pts); k++){
		Aside[j][k] = 0.0;
	      }
	    }
	   // Evaluate the second three rows of the matrix (zero, first and second derivative at point zero).
	    Aside[0][0] = 1.0 ;

	    // For the both first and second derivative it does not really matter what value is placed there.
	    Aside[1][1] = 1.0 ;
	    Aside[2][2] = 2.0 ;
	  
	  }
	  else {
	   for (int j=0; j<(6+pts); j++){
	     Aside[0][j] = pow(INTERVALS_4_APPROX[i-1],j);
	     Aside[1][j] = j*pow(INTERVALS_4_APPROX[i-1], j-1);
	     Aside[2][j] = j*(j-1)*pow(INTERVALS_4_APPROX[i-1], j-2);
	   }
	  }
	  for (int j=0; j<(6+pts); j++){	   
	    Aside[3][j] = pow(INTERVALS_4_APPROX[i],j);
	    Aside[4][j] = j*pow(INTERVALS_4_APPROX[i], j-1);
	    Aside[5][j] =  j*(j-1)*pow(INTERVALS_4_APPROX[i], j-2);
	    switch (pts) {
	      case 0 :
		break;
	      case 1 : 
		if (i > 0) {
		  Aside[6][j] =  pow((INTERVALS_4_APPROX[i]+INTERVALS_4_APPROX[i-1])/2, j);
		}
		else {
		  Aside[6][j] =  pow((INTERVALS_4_APPROX[i])/2, j) ;
		}
		break;
	      default :
		printf("too many of inside points in polynomial approximation; not implemented!!! (check extras.cpp) \n");
		exit(1);
	    }
	  }
	  //Evaluate K(h) function.
	  if (i<1){
	    Bside[0] = K(0.0, layer);
	    Bside[1] = 0.0;
	    Bside[2] = 0.0;
	  }
	  else {
	    Bside[0] = K(INTERVALS_4_APPROX[i-1], layer);
	    Bside[1] = dKdh(INTERVALS_4_APPROX[i-1], layer);
	    Bside[2] = ddKdhh(INTERVALS_4_APPROX[i-1], layer);
	  }
	 
	  Bside[3] = K(INTERVALS_4_APPROX[i], layer);
	  Bside[4] = dKdh(INTERVALS_4_APPROX[i], layer);
	  Bside[5] = ddKdhh(INTERVALS_4_APPROX[i], layer);  
	  
	  switch (pts) {
	      case 0 :
		break;
	      case 1 :
		if (i > 0) {
		  Bside[6] = K((INTERVALS_4_APPROX[i]+INTERVALS_4_APPROX[i-1])/2, layer);
		}
		else {
		  Bside[6] = K((INTERVALS_4_APPROX[i])/2, layer);
		}
		break;
	  }



	  gem_full(Aside, Bside, K_POLS[i][layer][0], (6+pts));

	  for (int j=1; j<3; j++){
	    for (int k=0; k<5; k++){
	      K_POLS[i][layer][j][k] = (k+1)*K_POLS[i][layer][j-1][k+1];
	    }
	    K_POLS[i][layer][j][6-j] = 0.0;
	  }
	  
	  //Evaluate C(h) functions.
	  if (i<1){
	    Bside[0] = C(0.0, layer);
	    Bside[1] = dCdh(0.0, layer);
	  }
	  else {
	    Bside[0] = C(INTERVALS_4_APPROX[i-1], layer);
	    Bside[1] = dCdh(INTERVALS_4_APPROX[i-1], layer);
	  }
	 
	  //The first two lines of the matrix Aside stays the same.
	  for (int j=0; j<4; j++){	   
	    Aside[2][j] = pow(INTERVALS_4_APPROX[i],j);
	    Aside[3][j] = j*pow(INTERVALS_4_APPROX[i], j-1);
	  }

	  
	  Bside[2] = C(INTERVALS_4_APPROX[i], layer);
	  Bside[3] = dCdh(INTERVALS_4_APPROX[i], layer);
	 
	  
	  gem_full(Aside, Bside, C_POLS[i][layer][0], 4);
	  
	  for (int k=0; k<5; k++){
	    C_POLS[i][layer][1][k] = (k+1)*C_POLS[i][layer][0][k+1];
	  }
	  
	  C_POLS[i][layer][1][5] = 0.0;
	
    }
    delete [] Aside;
    delete [] Bside;
    break ;
  }
      
 
  return true;
}







 





