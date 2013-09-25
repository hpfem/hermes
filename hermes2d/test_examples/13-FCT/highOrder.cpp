#include "highOrder.h"

using namespace Hermes;

High_Order::High_Order(double theta): theta(theta)
{
  highmat_rhs =new CSCMatrix<double> ;;
  high_matrix =new CSCMatrix<double> ;;
  u_H =NULL;
}

High_Order::~High_Order()
{
  if(high_matrix!=NULL) delete high_matrix;
  if(highmat_rhs!=NULL) delete highmat_rhs;
  if(u_H!=NULL) delete [] u_H;
}

//high-order scheme:  (M\tau -theta K) u_H   = (M\tau + (1-theta) K) u_n 
// high_matrix = (M\tau -theta K); highmat_rhs = (M\tau + (1-theta) K); M/tau = mass_matrix; K = conv_matrix
void High_Order::assemble_High_Order(CSCMatrix<double> * conv_matrix, CSCMatrix<double> * mass_matrix)
{
  if(high_matrix!=NULL) 
    high_matrix->free();
  if(highmat_rhs!=NULL) 
    highmat_rhs->free();
  high_matrix->create(conv_matrix->get_size(),conv_matrix->get_nnz(), conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
  high_matrix->multiply_with_Scalar(-theta);
  high_matrix->add_sparse_matrix(mass_matrix);  
  highmat_rhs->create(conv_matrix->get_size(),conv_matrix->get_nnz(), conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
  if(theta!=0.) 
    highmat_rhs->multiply_with_Scalar((1.0-theta));
  highmat_rhs->add_sparse_matrix(mass_matrix); 
}

// (M\tau -theta K) u_H   = (M\tau + (1-theta) K) u_n 
double* High_Order::solve_High_Order(double* u_n)
{
  int ndof = high_matrix->get_size();
  if((high_matrix==NULL)||(highmat_rhs==NULL) )
    throw Exceptions::Exception("matrices have to be calculated first, call High_Order::assemble_High_Order().");
  if(u_H!=NULL) delete [] u_H;
  u_H = new double[ndof];			
  double* coeff_vec_2 = new double[ndof];
  SimpleVector<double> * vec_rhs = new SimpleVector<double> (ndof);	

  //----------vec_rhs = (M\tau + (1-theta) K) u_n  --------
  highmat_rhs->multiply_with_vector(u_n, coeff_vec_2, true); 
  vec_rhs->add_vector(coeff_vec_2);

  //// Solve the linear system and if successful, obtain the solution.   (M\tau -theta K) u_H   = vec_rhs
  UMFPackLinearMatrixSolver<double> * highOrd = new UMFPackLinearMatrixSolver<double> (high_matrix,vec_rhs);	
  try
  {
    highOrd->solve();
  }
  catch(Hermes::Exceptions::Exception e)
  {
    e.print_msg();
  }	

  for(int i = 0; i<ndof;i++) 	u_H[i] = highOrd->get_sln_vector()[i];  

  delete highOrd;
  delete [] coeff_vec_2;
  delete vec_rhs;

  return u_H;


}


