#include "lowOrder.h"

using namespace Hermes;

Low_Order::Low_Order(double theta): theta(theta)
{
	lowmat_rhs = new CSCMatrix<double> ; 
	low_matrix = new CSCMatrix<double> ; 
	u_L =NULL;
	rhs = NULL;
	u_new = NULL;
}


Low_Order::~Low_Order()
{
	if(low_matrix!=NULL) delete low_matrix;
	if(lowmat_rhs!=NULL) delete lowmat_rhs;
	if(u_L!=NULL) delete [] u_L;
	if(rhs!=NULL) delete [] rhs;
	if(u_new!=NULL) delete [] u_new;
}


//low-order: (M_L/tau - theta (K+D)) u_n+1 = (M_L/tau - theta (K+D)) u_n
//lowmat_rhs = (M_L/tau - theta (K+D)) u_n; low_matrix: (M_L/tau - theta (K+D)) ; M_L/tau = lumped_matrix; K = conv_matrix, D = diffusion
void Low_Order::assemble_Low_Order(CSCMatrix<double> * conv_matrix,CSCMatrix<double>* diffusion, CSCMatrix<double> * lumped_matrix)
{
	if(low_matrix!=NULL) low_matrix->free();	 		
	if(lowmat_rhs!=NULL) lowmat_rhs->free();

	lowmat_rhs->create(conv_matrix->get_size(),conv_matrix->get_nnz(), conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
	lowmat_rhs->add_sparse_matrix(diffusion); 
	low_matrix->create(lowmat_rhs->get_size(),lowmat_rhs->get_nnz(), lowmat_rhs->get_Ap(), lowmat_rhs->get_Ai(),lowmat_rhs->get_Ax());

	//(-theta)(K+D)
	if(theta==0) low_matrix->zero();
	else	low_matrix->multiply_with_Scalar(-theta);
	//(1-theta)(K+D)
	if(theta ==1) lowmat_rhs->zero();
	else lowmat_rhs->multiply_with_Scalar((1.0-theta));
	//M_L/tau - theta(D+K)
	low_matrix->add_sparse_matrix(lumped_matrix);  
	//M_L/tau+(1-theta)(K+D)
	lowmat_rhs->add_sparse_matrix(lumped_matrix);	
}

// solve: (M_L/tau) u_L = (M_L/tau - theta (K+D)) u_n
//lumped_matrix = M_L !!!!!
double* Low_Order::solve_Low_Order(CSCMatrix<double> * lumped_matrix, double* u_n, double time_step  )
{
			if((low_matrix==NULL)||(lowmat_rhs==NULL) )
				throw Exceptions::Exception("matrices have to be calculated first, call  Low_Order::assemble_Low_Order.");
			int ndof = low_matrix->get_size();
			if(u_L!=NULL) delete [] u_L;
			u_L = new double[ndof];
			if(rhs!=NULL) delete [] rhs;
			rhs = new double[ndof];
			SimpleVector<double>* vec_rhs = new SimpleVector<double> (ndof); 
			
						lumped_matrix->multiply_with_Scalar(1./time_step); //M_L/tau
				//----------vec_rhs = (M_L/tau - theta (K+D)) u_n  --------
			lowmat_rhs->multiply_with_vector(u_n, rhs, true); 
			 vec_rhs->add_vector(rhs);
			 
// Solve the linear system and if successful, obtain the solution. M_L/tau u^L=  M_L/tau+ (1-theta)(K+D) u^n
			UMFPackLinearMatrixSolver<double> * lowOrd = new UMFPackLinearMatrixSolver<double> (lumped_matrix,vec_rhs);	
			try{
				lowOrd->solve();
			}
			catch(Hermes::Exceptions::Exception e)
			{
				 e.print_msg();
			}	
				for(int i = 0; i<ndof;i++) u_L[i] = lowOrd->get_sln_vector()[i];  
			lumped_matrix->multiply_with_Scalar(time_step);  // M_L
			
		delete lowOrd;
		delete vec_rhs;
		
		return u_L;

}

//(M_L/tau - theta (K+D)) u_new = (M_L/tau - theta (K+D)) u_n + flux_correction
double* Low_Order::explicit_Correction(double* flux_correction )
{
			if((low_matrix==NULL)||(lowmat_rhs==NULL) )
				throw Exceptions::Exception("matrices have to be calculated first, call  Low_Order::assemble_Low_Order.");
			if(rhs ==NULL) 
				throw Exceptions::Exception("rhs have to be calculated first, call Low_Order::solve_Low_Order() ");
			if(u_new!=NULL) delete [] u_new;
			int ndof = low_matrix->get_size();
			u_new = new double[ndof];
			SimpleVector<double>* vec_rhs = new SimpleVector<double> (ndof); 
			//rhs = (M_L/tau - theta (K+D)) u_n 			
			vec_rhs->add_vector(rhs);
			vec_rhs->add_vector(flux_correction);
			
	// Solve the linear system and if successful, obtain the solution. (M_L/tau - theta (K+D)) u_new = rhs+flux_correction	
			UMFPackLinearMatrixSolver<double> * newSol = new UMFPackLinearMatrixSolver<double> (low_matrix,vec_rhs);
			try{
				newSol->solve();
			}
			catch(Hermes::Exceptions::Exception e)
			{
				 e.print_msg();
			}	

			for(int i = 0; i<ndof;i++) u_new[i] = newSol->get_sln_vector()[i];  

			delete newSol;
			delete vec_rhs;
				
			return u_new;
}

