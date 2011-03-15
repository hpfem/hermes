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