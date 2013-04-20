#include "fct.h"

Flux_Correction::Flux_Correction(double theta)
{
  this->theta = theta;
  fct = NULL; 
  space =NULL;			
  al = new AsmList<double>;	 
  P_plus= NULL;  P_minus= NULL;  Q_plus= NULL;  Q_minus= NULL;   R_plus= NULL;  R_minus= NULL;
};
Flux_Correction::~Flux_Correction()
{
  free();
  delete al;

};

void Flux_Correction::free()
{
  if(P_plus!=NULL)  delete [] P_plus;
  if(P_minus!=NULL) delete [] P_minus;
  if(Q_plus!=NULL)  delete [] Q_plus;
  if(Q_minus!=NULL) delete [] Q_minus;
  if(R_plus!=NULL)  delete [] R_plus;
  if(R_minus!=NULL) delete [] R_minus;
  if(fct!=NULL) delete [] fct; 

};



void Flux_Correction::init(SpaceSharedPtr<double> new_space)
{
  free();
  space = new_space;					
  int ndof = space->get_num_dofs();
  fct = new bool[ndof];	
  P_plus = new double[ndof]; P_minus = new double[ndof];
  Q_plus = new double[ndof]; Q_minus = new double[ndof];	
  R_plus = new double[ndof]; R_minus = new double[ndof];
  for(int i=0; i<ndof;i++)
    fct[i]=false;	
  Element* e =NULL;
  Element* elem_neigh=NULL;
  bool more = false;
  bool  p2_neighbor =false;
  int elem_id,id;

  //Determine which DOFs can be used for FCT/masslumping/artificial Diffusion
  for_all_active_elements(e, space->get_mesh())
  {  
    if((space->get_element_order(e->id)== H2D_MAKE_QUAD_ORDER(1, 1))||(space->get_element_order(e->id)==1))
    {
      elem_id= e->id;
      for (unsigned int iv = 0; iv < e->get_nvert(); iv++)
      {  
        // Class for finding all neighbors of the element e on the mesh space->get_mesh().
        NeighborSearch<double> ns(e, space->get_mesh());
        // Set this to ignore errors of the edge iv being a boundary edge (which has no neighbors).
        ns.set_ignore_errors(true);
        // Look for the neighbors over the edge iv.
        ns.set_active_edge(iv);
        // Iterate through the found neighbors.
        for(unsigned int i = 0; i < ns.get_neighbors()->size(); i++)
        {
          id = ns.get_neighbors()->at(i)->id;
          if((space->get_element_order(id)== H2D_MAKE_QUAD_ORDER(2, 2))||(space->get_element_order(id)==2))
          {
            p2_neighbor =true; 
            break;
          }
        }
        if(e->vn[iv]->is_constrained_vertex() ==true)	
        {		
          p2_neighbor =true; 
          break;
        }
      }

      if(p2_neighbor==false)
      {
        space->get_element_assembly_list(e, al);
        for (unsigned int iv = 0; iv < e->get_nvert(); iv++)
        {   		
          int index = space->get_shapeset()->get_vertex_index(iv,HERMES_MODE_QUAD);
          Node* vn = e->vn[iv];
          if (space->get_element_order(elem_id) == 0) break;
          if (!vn->is_constrained_vertex())
          {  
            for(unsigned int j = 0; j < al->get_cnt(); j ++)
            {			 
              if((al->get_idx()[j]==index)&&(al->get_dof()[j]!=-1.0))
              { 
                if(fct[al->get_dof()[j]]==false)
                  fct[al->get_dof()[j]]=true;															
              }
            }
          }
        }					
      } 
      else 
        p2_neighbor =false;							
    }
  }

}



UMFPackMatrix<double>* Flux_Correction::artificialDiffusion(UMFPackMatrix<double>* conv_matrix)
{
  if(fct==NULL) 
    throw Exceptions::Exception("fct-list=NULL");
  int size = conv_matrix->get_size();
  int nnz = conv_matrix->get_nnz();
  double a,b;
  UMFPackMatrix<double>* diffusion = new UMFPackMatrix<double>;  
  diffusion->create(size, nnz, conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
  diffusion->zero();  

  double* Ax = conv_matrix->get_Ax();
  int* Ai = conv_matrix->get_Ai();
  int* Ap = conv_matrix->get_Ap();

  for(int j = 0; j<size; j++)
  { //iterate through columns
    if(fct[j]== false) continue;
    for(int indx = Ap[j]; indx<Ap[j+1];indx++)
    {	
      int i = Ai[indx];	
      if(fct[i]== false) continue;
      if(j<i)
      {
        a= -Ax[indx];
        b= -conv_matrix->get(j,i);     
        if((a>=b)&&(a>0.0))
        {
          diffusion->add(j,i,a);
          diffusion->add(i,j,a);	
          diffusion->add(j,j,-a);
          diffusion->add(i,i,-a);	
        }
        else if((b>a)&&(b>0.0))
        {
          diffusion->add(i,j,b);
          diffusion->add(j,i,b);
          diffusion->add(j,j,-b);
          diffusion->add(i,i,-b);
        }
      }
    }
  }


  return diffusion;

}


UMFPackMatrix<double>* Flux_Correction::massLumping(UMFPackMatrix<double>* mass_matrix)
{  
  if(fct==NULL)
    throw Exceptions::Exception("fct-list=NULL");
  UMFPackMatrix<double>* lumped_matrix = new UMFPackMatrix<double>;   //M_L
  int size = mass_matrix->get_size();
  int nnz = mass_matrix->get_nnz();
  lumped_matrix->create(size, nnz, mass_matrix->get_Ap(), mass_matrix->get_Ai(),mass_matrix->get_Ax());
  double* Ax = lumped_matrix->get_Ax();
  int* Ai = lumped_matrix->get_Ai();
  int* Ap = lumped_matrix->get_Ap();
  double a =0.0;

  for(int j = 0; j<size; j++)
  { //iterate through columns
    if(fct[j]== false) continue;
    for(int indx = Ap[j]; indx<Ap[j+1];indx++)
    {	
      int i = Ai[indx];	
      if(fct[i]== false) continue;
      if(j<i)
      {
        a = Ax[indx];
        if(a!=0.0)
        {
          lumped_matrix->add(i,i,a);    
          lumped_matrix->add(j,j,a);    
          lumped_matrix->add(i,j,-a);	  
          lumped_matrix->add(j,i,-a);	
        }	
      }
    }
  }
  return lumped_matrix;
}



//Assemble antidiffusive fluxes & limit these
void Flux_Correction::antidiffusiveFlux(UMFPackMatrix<double>* mass_matrix,UMFPackMatrix<double>* lumped_matrix,UMFPackMatrix<double>* conv_matrix,UMFPackMatrix<double>* diffusion,double* u_high, double* u_L, double* u_old,double* flux_scalar,double time_step,Regularity_Estimator* regEst)
{ 
  if(fct==NULL) throw Exceptions::Exception("fct-list=NULL");

  int* smooth_dof =NULL;
  if(regEst!=NULL) smooth_dof = regEst->get_smooth_dofs(space,u_L,mass_matrix);		

  int ndof = conv_matrix->get_size();
  double alpha,f, plus, minus,mass, diff;
  double* Ax_mass = mass_matrix->get_Ax();
  int* Ai_mass = mass_matrix->get_Ai();
  int* Ap_mass = mass_matrix->get_Ap();

  for(int i=0; i<ndof;i++){ P_plus[i]=0.0; P_minus[i]=0.0; Q_plus[i]=0.; Q_minus[i]=0.; flux_scalar[i]=0.0;}

  //Compute P&Q
  for(int j = 0; j<ndof; j++)
  { 
    if(fct[j]== false) continue;
    for(int indx = Ap_mass[j]; indx<Ap_mass[j+1];indx++)
    {	
      int i = Ai_mass[indx];	
      if(fct[i]== false) continue;
      if(((mass=Ax_mass[indx])!=0.)&&(j<i))
      {
        diff = diffusion->get(i,j);
        f = (mass/time_step+ diff*theta)*(u_high[i]- u_high[j])
          -(mass/time_step - diff*(1.-theta)) *(u_old[i]- u_old[j]);	
        if( (f*(u_L[j]- u_L[i])) > 0.0) f = 0.0; //prelimiting step										
        if(f>0.0)	
        { 
          P_plus[i]+=f;
          P_minus[j]-=f;
        }
        else if (f<0.0)
        {
          P_minus[i]+=f;
          P_plus[j]-=f;
        }										
        f = lumped_matrix->get(i,i)*(u_L[j]-u_L[i])/time_step; 
        if(f>Q_plus[i]) Q_plus[i] = f;				
        if(f<Q_minus[i]) Q_minus[i] = f;			
        f= lumped_matrix->get(j,j)*(u_L[i]-u_L[j])/time_step;
        if(f>Q_plus[j]) Q_plus[j] = f;	
        if(f<Q_minus[j]) Q_minus[j] = f;
      }							
    }			
  }

  //Compute R	
  for(int i=0; i<ndof;i++)
  {
    if(fct[i]== false) continue;
    plus = 1.0; minus = 1.0;		
    if(P_plus[i]!=0.0)  plus = Q_plus[i]/P_plus[i];		
    if(P_minus[i]!=0.0) minus = Q_minus[i]/P_minus[i];			
    if(plus>=1.0) R_plus[i]= 1.0;
    else 	     R_plus[i]= plus;
    if(minus>=1.0) R_minus[i]= 1.0;
    else 	     R_minus[i]= minus;
    if(smooth_dof!=NULL)
    {
      if(smooth_dof[i]==1){ R_plus[i]= 1.0;R_minus[i]= 1.0;}
    }	
  }

  //Compute alpha & f_i
  alpha = 1.0;
  for(int j = 0; j<ndof; j++)
  { 
    if(fct[j]== false) continue;
    for(int indx = Ap_mass[j]; indx<Ap_mass[j+1];indx++)
    {	
      int i = Ai_mass[indx];	
      if(fct[i]== false) continue;
      if(((mass=Ax_mass[indx])!=0.)&&(j<i))
      {
        diff = diffusion->get(i,j);
        f = (mass/time_step+ diff*theta)*(u_high[i]- u_high[j])
          -(mass/time_step - diff*(1.-theta)) *(u_old[i]- u_old[j]);	
        if( (f*(u_L[j]- u_L[i])) > 0.0) f = 0.0; //prelimiting step										
        if(f>0.)
        {					
          if(R_plus[i]>R_minus[j]) alpha = R_minus[j];
          else 	alpha = R_plus[i];
        }
        else if(f<0.)
        {
          if(R_minus[i]>R_plus[j]) alpha = R_plus[j];
          else 	alpha = R_minus[i]; 
        }									
        flux_scalar[i] += alpha*f;
        flux_scalar[j] -= alpha*f;		
      }
    }			
  }


}

//FCT for projection 
void Flux_Correction::lumped_flux_limiter(UMFPackMatrix<double>* mass_matrix,UMFPackMatrix<double>* lumped_matrix, double* u_L, double* u_H, double time_step, int* smooth_dof)
{	
  if(fct==NULL) 
    throw Exceptions::Exception("fct-list=NULL");
  int ndof = mass_matrix->get_size();
  double* rhs = new double[ndof];
  lumped_matrix->multiply_with_vector(u_L, rhs); 

  double alpha,f, plus, minus,mass;

  double* Ax_mass = mass_matrix->get_Ax();
  int* Ai_mass = mass_matrix->get_Ai();
  int* Ap_mass = mass_matrix->get_Ap();


  for(int i=0; i<ndof;i++){ P_plus[i]=0.0;P_minus[i]=0.0;Q_plus[i]=0.;Q_minus[i]=0.;}
  //Compute P&Q
  for(int j = 0; j<ndof; j++)
  { 
    if(fct[j]== false) continue;
    for(int indx = Ap_mass[j]; indx<Ap_mass[j+1];indx++)
    {	
      int i = Ai_mass[indx];
      if(fct[i]== false) continue;	
      if(((mass=Ax_mass[indx])!=0.)&&(j<i))
      {
        f = mass*(u_H[i]- u_H[j]);
        if( (f*(u_L[j]- u_L[i])) > 0.0) f = 0.0; //prelimiting step										
        if(f>0.0)	
        { 
          P_plus[i]+=f;
          P_minus[j]-=f;
        }
        else if(f<0.0)
        {
          P_minus[i]+=f;
          P_plus[j]-=f;
        }										
        f = lumped_matrix->get(i,i)*(u_L[j]-u_L[i]); 
        if(f>Q_plus[i]) Q_plus[i] = f;				
        if(f<Q_minus[i]) Q_minus[i] = f;			
        f= lumped_matrix->get(j,j)*(u_L[i]-u_L[j]);
        if(f>Q_plus[j]) Q_plus[j] = f;	
        if(f<Q_minus[j]) Q_minus[j] = f;	
      }
    }			
  }


  //Compute R
  for(int i=0; i<ndof;i++)
  {
    if(fct[i]== false) continue;	
    plus = 1.0; minus = 1.0;		
    if(P_plus[i]!=0.0)  plus = Q_plus[i]/P_plus[i];		
    if(P_minus[i]!=0.0) minus = Q_minus[i]/P_minus[i];			
    if(plus>=1.0) R_plus[i]= 1.0;
    else 	     R_plus[i]= plus;
    if(minus>=1.0) R_minus[i]= 1.0;
    else 	     R_minus[i]= minus;	
    if(smooth_dof!=NULL)
    {
      if(smooth_dof[i]==1){ R_plus[i]= 1.0;R_minus[i]= 1.0;}
    }		
  }	

  //Compute alpha & f_i
  alpha = 1.0;
  for(int j = 0; j<ndof; j++)
  { 
    if(fct[j]== false) continue;
    for(int indx = Ap_mass[j]; indx<Ap_mass[j+1];indx++)
    {	
      int i = Ai_mass[indx];	
      if(fct[i]== false) continue;
      if(((mass=Ax_mass[indx])!=0.)&&(j<i))
      {
        f = mass*(u_H[i]- u_H[j]);
        if( (f*(u_L[j]- u_L[i])) > 0.0) f = 0.0; //prelimiting step										
        if(f>0.)
        {					
          if(R_plus[i]>R_minus[j]) alpha = R_minus[j];
          else 	alpha = R_plus[i];
        }
        else if(f<0.)
        {
          if(R_minus[i]>R_plus[j]) alpha = R_plus[j];
          else 	alpha = R_minus[i]; 
        }									
        rhs[i] += alpha*f;
        rhs[j] -= alpha*f;									

      }
    }			
  }



  UMFPackVector<double>* vec_rhs = new UMFPackVector<double>(ndof);	
  vec_rhs->zero(); vec_rhs->add_vector(rhs);
  double* sol =NULL;
  UMFPackLinearMatrixSolver<double>* lowOrd = new UMFPackLinearMatrixSolver<double>(lumped_matrix,vec_rhs);	
  try
  {
    lowOrd->solve();
  }
  catch(Hermes::Exceptions::Exception e)
  {
    e.print_msg();
  }		

  for(int i=0; i<ndof;i++) u_L[i] =lowOrd->get_sln_vector()[i];

  delete lowOrd;	
  delete vec_rhs;
  delete [] rhs;
}

void Flux_Correction::project_FCT(MeshFunctionSharedPtr<double> sln,double* coeff_vec, double* coeff_vec_2,UMFPackMatrix<double>* mass_matrix,UMFPackMatrix<double>* lumped_matrix, double time_step,OGProjection<double>* ogProjection,	Lumped_Projection* lumpedProjection, Regularity_Estimator* regEst)
{
  if(sln==NULL)
    throw Exceptions::Exception("project_FCT: sln=NULL");
  int* smooth_dof =NULL;
  //low-order projection
  lumpedProjection->project_lumped(space, sln, coeff_vec, lumped_matrix);			
  if(regEst!=NULL)
    smooth_dof = regEst->get_smooth_dofs(space,coeff_vec,mass_matrix);			//=>in smooth elements no FCT needs to be applied
  //high-order projection		
  ogProjection->project_global(space, sln, coeff_vec_2, HERMES_L2_NORM);	
  //calculation of fluxes, limiting & explicit correction				
  lumped_flux_limiter(mass_matrix, lumped_matrix, coeff_vec, coeff_vec_2,time_step,smooth_dof); 
}


