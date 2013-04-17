#include "lumped_projection.h"


void Lumped_Projection::project_internal(SpaceSharedPtr<double> space, WeakForm<double>* wf, double* target_vec, UMFPackMatrix<double>*  mat)
{

  if(space == NULL) printf("this->space == NULL in Lumped_Projection::project_internal().");
  // Get dimension of the space.
  int ndof = space->get_num_dofs();

  if(mat!=NULL) if(mat->get_size()!=ndof) printf("matrix=%i, ndof=%i", mat->get_size(),ndof);

  // Initialize DiscreteProblem.
  DiscreteProblem<double>* dp = new DiscreteProblem<double>(wf, space);
  UMFPackMatrix<double>* matrix = new UMFPackMatrix<double>;	
  UMFPackVector<double>* rhs = new UMFPackVector<double>(ndof);
  double* coeff_vec =NULL; 
  if(mat==NULL) 		//=> masslumping	
  {
    UMFPackMatrix<double>* lumped_matrix = new UMFPackMatrix<double>;   //M_L 
    dp->assemble(matrix, rhs);  		 
    int size = matrix->get_size();
    double* diag = new double[size];
    int nnz = matrix->get_nnz();
    int* row = new int[size]; 
    int* col = new int[size+1];
    for(int i = 0; i<size; i++)
    {    
      diag[i] = 0;
      row[i]= i;
      col[i]=i;
    }
    col[size]=size;

    for(int i = 0; i<nnz; i++)    
      diag[matrix->get_Ai()[i]] += matrix->get_Ax()[i]; 
    lumped_matrix->create(size, size, col, row, diag);  
    UMFPackLinearMatrixSolver<double>* solver = new UMFPackLinearMatrixSolver<double>(lumped_matrix,rhs);		
    try
    {
      solver->solve();
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.print_msg();
    }
    coeff_vec = solver->get_sln_vector();	

    if (target_vec != NULL)
      for (int i=0; i < ndof; i++) 
        target_vec[i] = coeff_vec[i];
    delete solver;
    delete lumped_matrix;
    delete [] diag;
    delete [] row;
    delete [] col;
  }
  else
  { 
    dp->assemble(rhs);
    UMFPackLinearMatrixSolver<double>* solver = new UMFPackLinearMatrixSolver<double>(mat,rhs);	

    try
    {
      solver->solve();
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.print_msg();
    }

    coeff_vec = solver->get_sln_vector();			

    if (target_vec != NULL)
      for (int i=0; i < ndof; i++) target_vec[i] = coeff_vec[i];
    delete solver;
  }


  delete matrix;
  delete rhs;
  delete dp;


}

void Lumped_Projection::project_lumped( const  SpaceSharedPtr<double> space, MeshFunctionSharedPtr<double> source_meshfn,
  double* target_vec, UMFPackMatrix<double>*  mat )
{
  // Sanity checks.
  if (target_vec == NULL) throw Exceptions::NullException(3);
  // Define temporary projection weak form.
  WeakForm<double>* proj_wf = new WeakForm<double>(1);
  proj_wf->set_ext(source_meshfn);

  ProjectionLumpedMatrixFormVol* matrix_form =	new ProjectionLumpedMatrixFormVol(0, 0);
  ProjectionLumpedVectorFormVol* vector_form = new ProjectionLumpedVectorFormVol(0);  
  //ProjectionLumpedVectorFormVol* vector_form = new ProjectionLumpedVectorFormVol(0, source_meshfn);  
  proj_wf->add_matrix_form(matrix_form);
  proj_wf->add_vector_form(vector_form);
  // Call main function.
  project_internal(space, proj_wf, target_vec, mat);

  // Clean up.
  delete proj_wf;

}