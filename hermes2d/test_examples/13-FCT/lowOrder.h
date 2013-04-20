#ifndef __LOWORDER_H
#define __LOWORDER_H

#include "hermes2d.h"


class Low_Order{

public:
	//Constructor
	Low_Order(double theta);

	//Destructor
	~Low_Order();

	void assemble_Low_Order(UMFPackMatrix<double> * conv_matrix,UMFPackMatrix<double>* diffusion, UMFPackMatrix<double> * lumped_matrix);
	double* solve_Low_Order(UMFPackMatrix<double> * lumped_matrix, double* u_n, double time_step);
	double* explicit_Correction(double* flux_correction);


protected:

	UMFPackMatrix<double> * low_matrix;  
	UMFPackMatrix<double> * lowmat_rhs;
	double* u_L;
	double* rhs;
	double* u_new;	
	double theta;


};
























#endif

