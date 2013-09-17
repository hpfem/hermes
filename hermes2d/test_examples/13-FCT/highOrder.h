#ifndef __HIGHORDER_H
#define __HIGHORDER_H

#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Solvers;
using namespace Hermes::Algebra;
using namespace Hermes::Hermes2D;

class High_Order{
public:
	//Constructor
	High_Order(double theta);

	//Destructor
	~High_Order();

	void assemble_High_Order(CSCMatrix<double> * conv_matrix, CSCMatrix<double> * mass_matrix);
	double* solve_High_Order(double* u_n);


protected:

	CSCMatrix<double> * high_matrix;  
	CSCMatrix<double> * highmat_rhs;
	double* u_H;
	double theta;

};






#endif

