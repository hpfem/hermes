#ifndef __HIGHORDER_H
#define __HIGHORDER_H

#include "hermes2d.h"

class High_Order{
public:
	//Constructor
	High_Order(double theta);

	//Destructor
	~High_Order();

	void assemble_High_Order(UMFPackMatrix<double> * conv_matrix, UMFPackMatrix<double> * mass_matrix);
	double* solve_High_Order(double* u_n);


protected:

	UMFPackMatrix<double> * high_matrix;  
	UMFPackMatrix<double> * highmat_rhs;
	double* u_H;
	double theta;

};






#endif

