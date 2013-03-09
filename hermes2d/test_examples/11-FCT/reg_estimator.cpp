#include "reg_estimator.h"

Regularity_Estimator::Regularity_Estimator(double epsilon): epsilon(epsilon), R_h_1(MeshFunctionSharedPtr<double>(new Solution<double>)), 
  R_h_2(MeshFunctionSharedPtr<double>(new Solution<double>)), sln(MeshFunctionSharedPtr<double>(new Solution<double>))
{
  space =NULL;
  
  smooth_elem_patch= NULL;
  smooth_dof = NULL;
  al = new 	AsmList<double>;
  rhs_1 =NULL;
  rhs_2=NULL;
  grad_1= new GradientReconstruction_1(sln);
  grad_2= new GradientReconstruction_2(sln);
};

Regularity_Estimator::~Regularity_Estimator()
{
  free();
  delete al;	
  delete grad_1; delete grad_2;	
}

void Regularity_Estimator::free()
{
  if(smooth_dof!=NULL) delete [] smooth_dof;
  if(smooth_elem_patch!=NULL) delete [] smooth_elem_patch;
  if(rhs_1 !=NULL) delete rhs_1;
  if(rhs_2 !=NULL) delete rhs_2;
  smooth_elem_patch= NULL;
  smooth_dof = NULL;
  rhs_1 =NULL;
  rhs_2=NULL;
}

void Regularity_Estimator::set_space(SpaceSharedPtr<double> new_space)
{
  if(new_space == NULL)
    throw Hermes::Exceptions::Exception("smoothness_indicator: set_space");
  free();
  space = new_space;
  int ndof = space->get_num_dofs(); 
  smooth_dof = new int[ndof];  
  smooth_elem_patch = new int[new_space->get_mesh()->get_max_element_id()];
  rhs_1 = new UMFPackVector<double>(ndof);
  rhs_2 = new UMFPackVector<double>(ndof);
}

int* Regularity_Estimator::get_smooth_elems(SpaceSharedPtr<double> new_space,double* coeff_vec,UMFPackMatrix<double> * mass_matrix)
{
  if(new_space == NULL) throw Hermes::Exceptions::Exception("smoothness_indicator: space =NULL");
  if(coeff_vec==NULL) throw Hermes::Exceptions::Exception("smoothness_indicator: coeff_vec=NULL");
  if(space!=new_space) 
    set_space(new_space);
  Solution<double>::vector_to_solution(coeff_vec, space, sln);
  smoothness_indicator(mass_matrix);
  return smooth_elem_patch;
}

int* Regularity_Estimator::get_smooth_dofs(SpaceSharedPtr<double> new_space,double* coeff_vec,UMFPackMatrix<double> * mass_matrix)
{
  if(new_space == NULL) throw Hermes::Exceptions::Exception("smoothness_indicator: space =NULL");
  if(coeff_vec==NULL) throw Hermes::Exceptions::Exception("smoothness_indicator: coeff_vec=NULL");
  if(space!=new_space) 
    set_space(new_space);
  Solution<double>::vector_to_solution(coeff_vec, space, sln);		
  smoothness_indicator(mass_matrix);
  return smooth_dof;
}



//linear approximation of sln in the neighborhood of element e   (at (x_i,y_i))
// u_h_hat = u_h (p_c) + R_h(p_c) *(p - p_c)
// R_h = (R_h_1, R_h_2), p = (x,y)
double Regularity_Estimator::linear_approx(Element* e, double x_i, double y_i,double x_c, double y_c, MeshFunctionSharedPtr<double> sln)
{
  //center of reference element
  double x_c_ref = 0.;
  double y_c_ref = 0.; 

  double u_h_x_c = (dynamic_cast<Solution<double>*>(sln.get()))->get_ref_value(e, x_c_ref, y_c_ref, 0, 0);
  double u_h_hat = u_h_x_c + (dynamic_cast<Solution<double>*>(R_h_1.get()))->get_ref_value(e, x_c_ref, y_c_ref, 0, 0)*(x_i-x_c)+(dynamic_cast<Solution<double>*>(R_h_2.get()))->get_ref_value(e, x_c_ref, y_c_ref, 0, 0)*(y_i-y_c);

  return u_h_hat;
}

//linear approximation of the first component of the gradient in the neighborhood of element e   (at (x_i,y_i))
// u_h_hat_1 = u_h_dx (p_c) + grad(R_h_1)(p_c) *(p - p_c)
double Regularity_Estimator::linear_approx_dx(Element* e, double x_i, double y_i,double x_c, double y_c,MeshFunctionSharedPtr<double> sln)
{
  //center of reference element	
  double x_c_ref = 0.;
  double y_c_ref = 0.; 

  double d_u_h_x_c = (dynamic_cast<Solution<double>*>(sln.get()))->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 1);

  double u_h_hat = d_u_h_x_c + (dynamic_cast<Solution<double>*>(R_h_1.get()))->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 1)*(x_i-x_c)
    + (dynamic_cast<Solution<double>*>(R_h_1.get()))->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 2)*(y_i-y_c);
  return u_h_hat;
}

//linear approximation of the second component of the gradient in the neighborhood of element e   (at (x_i,y_i))
// u_h_hat_2 = u_h_dy (p_c) + grad(R_h_2)(p_c) *(p - p_c)
double Regularity_Estimator::linear_approx_dy(Element* e, double x_i, double y_i,double x_c, double y_c,MeshFunctionSharedPtr<double> sln)
{

  //center of reference element	
  double x_c_ref = 0.;
  double y_c_ref = 0.; 

  double d_u_h_x_c = (dynamic_cast<Solution<double>*>(sln.get()))->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 2);

  double u_h_hat = d_u_h_x_c + (dynamic_cast<Solution<double>*>(R_h_2.get()))->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 1)*(x_i-x_c)
    +(dynamic_cast<Solution<double>*>(R_h_2.get()))->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 2)*(y_i-y_c);
  return u_h_hat;

}




void Regularity_Estimator::smoothness_indicator(UMFPackMatrix<double> * mass_matrix)
{
  int ndof = space->get_num_dofs(); 
  if((mass_matrix!=NULL)&&(mass_matrix->get_matrix_size() != ndof)) 
    throw Hermes::Exceptions::Exception("smoothness_indicator: mass_matrix size unequal ndof");


  for(int i=0; i<space->get_mesh()->get_max_element_id(); i++) smooth_elem_patch[i] = 4;		

  DiscreteProblem<double> * dp_1 = new DiscreteProblem<double> (grad_1, space);
  DiscreteProblem<double> * dp_2 = new DiscreteProblem<double> (grad_2, space);
  bool own_mass_matrix = false;
  if(mass_matrix ==NULL) 
  {
    mass_matrix = new UMFPackMatrix<double> ; 	
    dp_1->assemble(mass_matrix,rhs_1);
    own_mass_matrix = true; 
  }
  else 
    dp_1->assemble(rhs_1);

  dp_2->assemble(rhs_2);

  UMFPackLinearMatrixSolver<double> * solver_1 = new UMFPackLinearMatrixSolver<double> (mass_matrix,rhs_1);
  UMFPackLinearMatrixSolver<double> * solver_2 = new UMFPackLinearMatrixSolver<double> (mass_matrix,rhs_2);	


  try{
    solver_1->solve();
  }
  catch(Hermes::Exceptions::Exception e)
  {
    e.print_msg();
  }

  try{
    solver_2->solve();
  }
  catch(Hermes::Exceptions::Exception e)
  {
    e.print_msg();
  }	


  Solution<double> ::vector_to_solution(solver_1->get_sln_vector() , space, R_h_1);		//L2-proj of the first comp. of the gradient		
  Solution<double> ::vector_to_solution(solver_2->get_sln_vector() , space, R_h_2);		//L2-proj of the second comp. of the gradient

  std::list<int>* dof_elem_list = new std::list<int>[ndof];	

  double* u_c = new double[space->get_mesh()->get_max_element_id()];
  double* d_u_c_dx = new double[space->get_mesh()->get_max_element_id()];
  double* d_u_c_dy = new double[space->get_mesh()->get_max_element_id()];
  //center of reference element
  double x_c_ref = 0.;
  double y_c_ref = 0.; 
  Element* e =NULL; 
  int* index = new int[H2D_MAX_NUMBER_VERTICES];
  for(int i =0;i<H2D_MAX_NUMBER_VERTICES;i++)
    index[i] =  space->get_shapeset()->get_vertex_index(i,HERMES_MODE_QUAD);

  //determine value of the solution/gradient at the center of each element
  // for each dof determine list of elements
  for_all_active_elements(e, space->get_mesh())
  {
    u_c[e->id]= (dynamic_cast<Solution<double>*>(sln.get()))->get_ref_value(e, x_c_ref, y_c_ref, 0, 0);
    d_u_c_dx[e->id]= (dynamic_cast<Solution<double>*>(sln.get()))->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 1);
    d_u_c_dy[e->id]= (dynamic_cast<Solution<double>*>(sln.get()))->get_ref_value_transformed(e, x_c_ref, y_c_ref, 0, 2);
    space->get_element_assembly_list(e, al);
    for (unsigned int iv = 0; iv < e->get_nvert(); iv++)
    { 
      for(unsigned int j = 0; j < al->get_cnt(); j ++)
      {			 
        if((al->get_idx()[j]==index[iv])&&(al->get_dof()[j]!=-1.0))    							
          dof_elem_list[al->get_dof()[j]].push_back(e->id); 	

      }
    }	

  }

  //------------------------------------------------

  Node* vn=NULL;
  double x_c =0.; double y_c =0.; double u_h_x_c =0.;
  double* x = new double[H2D_MAX_NUMBER_VERTICES]; double* y = new double[H2D_MAX_NUMBER_VERTICES];
  double u_i; double u_dx ,u_dy;
  double u_min, u_max, u_min_dx, u_min_dy, u_max_dx, u_max_dy;
  std::list<int>::iterator elem_id; 
  int smooth_fct_in_elem,smooth_dx_in_elem,smooth_dy_in_elem;


  bool non_smooth = false; bool non_smooth_dx = false;bool non_smooth_dy = false;

  for_all_active_elements(e, space->get_mesh())
  {
    non_smooth = false; non_smooth_dx = false; non_smooth_dy = false;
    smooth_fct_in_elem =1; //assume as smooth	
    smooth_dx_in_elem =1; 
    smooth_dy_in_elem =1; 
    smooth_elem_patch[e->id]=0;
    space->get_element_assembly_list(e, al);		
    x_c = 0.; y_c = 0.;
    for (unsigned int iv = 0; iv < e->get_nvert(); iv++)
    { 		 
      vn = e->vn[iv];
      x[iv]=vn->x; y[iv]=vn->y;  //x/y-coord of nodes
      x_c+= (vn->x); y_c+= (vn->y);    // center
    }
    x_c/=e->get_nvert(); y_c/=e->get_nvert();
    for (unsigned int iv = 0; iv < e->get_nvert(); iv++)
    { 			
      for(unsigned int j = 0; j < al->get_cnt(); j ++)
      {			 
        if((al->get_idx()[j]==index[iv])&&(al->get_dof()[j]!=-1.0))
        {
          u_i = linear_approx(e, x[iv], y[iv],x_c, y_c,sln);	
          u_dx =	linear_approx_dx(e, x[iv], y[iv],x_c, y_c,sln);
          u_dy =	linear_approx_dy(e, x[iv], y[iv],x_c, y_c,sln);	
          u_min = u_c[e->id]; u_max =u_c[e->id];	
          u_min_dx= d_u_c_dx[e->id];	u_max_dx= d_u_c_dx[e->id];	
          u_min_dy= d_u_c_dy[e->id];	u_max_dy= d_u_c_dy[e->id];								
          for(elem_id=dof_elem_list[al->get_dof()[j]].begin();elem_id!=dof_elem_list[al->get_dof()[j]].end();elem_id++)
          {	
            if(u_c[*elem_id]>u_max) u_max =	u_c[*elem_id];
            else if(u_c[*elem_id]<u_min) u_min =	u_c[*elem_id];	
            if(d_u_c_dx[*elem_id]>u_max_dx) u_max_dx =	d_u_c_dx[*elem_id];
            else if(d_u_c_dx[*elem_id]<u_min_dx) u_min_dx =	d_u_c_dx[*elem_id];	
            if(d_u_c_dy[*elem_id]>u_max_dy) u_max_dy =	d_u_c_dy[*elem_id];
            else if(d_u_c_dy[*elem_id]<u_min_dy) u_min_dy =	d_u_c_dy[*elem_id];								
          }
          if((u_i>u_min+epsilon)&&(u_i<u_max-epsilon)) non_smooth = false;
          else non_smooth = true;
          if((u_dx>u_min_dx+epsilon)&&(u_dx<u_max_dx-epsilon)) non_smooth_dx = false;
          else non_smooth_dx = true;
          if((u_dy>u_min_dy+epsilon)&&(u_dy<u_max_dy-epsilon)) non_smooth_dy = false;
          else non_smooth_dy = true;			
          break;
        }
      }
      if(non_smooth == true) smooth_fct_in_elem=0; 
      if(non_smooth_dx == true) smooth_dx_in_elem=0; 
      if(non_smooth_dy == true) smooth_dy_in_elem=0; 
      if((non_smooth == true)&&(non_smooth_dx == true)&&(non_smooth_dy == true)) break;
    }

    if(std::max(smooth_fct_in_elem, std::min(smooth_dx_in_elem,smooth_dy_in_elem))==1)
      smooth_elem_patch[e->id]=1;


  }



  for(int i =0; i<ndof;i++)
  {
    non_smooth = false; 
    smooth_dof[i]=0;
    for(elem_id=dof_elem_list[i].begin();elem_id!=dof_elem_list[i].end();elem_id++)
    {
      if(smooth_elem_patch[*elem_id]==0){ non_smooth = true; break;}
    }
    if(non_smooth==true) 
    {
      for(elem_id=dof_elem_list[i].begin();elem_id!=dof_elem_list[i].end();elem_id++){
        if(smooth_elem_patch[*elem_id]==1)	smooth_elem_patch[*elem_id]=2;
      }			
    }
    else
      smooth_dof[i]=1;

  }




  //Clean-up
  delete solver_1;
  delete solver_2;
  delete [] u_c;
  delete [] d_u_c_dx ;
  delete [] d_u_c_dy ;
  delete [] index;
  delete [] x ;
  delete [] y;
  delete dp_1;
  delete dp_2;
  if(own_mass_matrix == true) delete mass_matrix; 
  delete [] dof_elem_list;	
}




//---------Gradient Reconstruction---------------
//R_H^1 
GradientReconstruction_1::GradientReconstruction_1( MeshFunctionSharedPtr<double> sln) : WeakForm<double>(1) 
{  this->set_ext(sln);
add_matrix_form(new GradientReconstructionMatForm_1(0, 0));
GradientReconstructionVectorForm_1* vector_form = new GradientReconstructionVectorForm_1(0);
add_vector_form(vector_form);
};



template<typename Real, typename Scalar>
Scalar GradientReconstructionMatForm_1 ::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
  Func<Real> *v, Geom<Real> *e, Func<Scalar>  **ext) const 
{

  Scalar result = Scalar(0); 
  for (int i = 0; i < n; i++) 
    result += wt[i] * (u->val[i] * v->val[i] );	
  return result;

};

double GradientReconstructionMatForm_1 ::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
  Func<double> *v, Geom<double> *e, Func<double>  **ext) const 
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
};

Ord GradientReconstructionMatForm_1 ::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
  Geom<Ord> *e, Func<Ord>  **ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
};

MatrixFormVol<double>* GradientReconstructionMatForm_1::clone() const
{
  return new GradientReconstructionMatForm_1(*this);
}

template<typename Real, typename Scalar>
Scalar GradientReconstructionVectorForm_1::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, Func<Scalar>  **ext) const {
  Scalar result = Scalar(0);
  Func<Scalar>* u_h = ext[0];
  for (int i = 0; i < n; i++) 	
    result += wt[i]*(v->val[i]*u_h->dx[i]);	
  return result;

};

double GradientReconstructionVectorForm_1::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>  **ext) const {
  return vector_form<double, double>(n, wt, u_ext, v, e, ext);
};

Ord GradientReconstructionVectorForm_1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord>  **ext) const 
{
  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
};
VectorFormVol<double>* GradientReconstructionVectorForm_1::clone() const
{
  return new GradientReconstructionVectorForm_1(*this);
}

//R_H^2
GradientReconstruction_2::GradientReconstruction_2( MeshFunctionSharedPtr<double> sln) : WeakForm<double>(1) 
{
  this->set_ext(sln);
  add_matrix_form(new GradientReconstructionMatForm_2(0, 0));
  GradientReconstructionVectorForm_2* vector_form = new GradientReconstructionVectorForm_2(0);
  add_vector_form(vector_form);
};



template<typename Real, typename Scalar>
Scalar GradientReconstructionMatForm_2 ::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
  Func<Real> *v, Geom<Real> *e, Func<Scalar>  **ext) const 
{
  Scalar result = Scalar(0); 
  for (int i = 0; i < n; i++) 
    result += wt[i] * (u->val[i] * v->val[i] );	
  return result;
};

double GradientReconstructionMatForm_2 ::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
  Func<double> *v, Geom<double> *e, Func<double>  **ext) const 
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
};
Ord GradientReconstructionMatForm_2 ::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
  Geom<Ord> *e, Func<Ord>  **ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
};

MatrixFormVol<double>* GradientReconstructionMatForm_2::clone() const
{
  return new GradientReconstructionMatForm_2(*this);
}

template<typename Real, typename Scalar>
Scalar GradientReconstructionVectorForm_2::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, Func<Scalar>  **ext) const {
  Scalar result = Scalar(0);
  Func<Scalar>* u_h = ext[0];
  for (int i = 0; i < n; i++) 	
    result += wt[i]*(v->val[i]*u_h->dy[i]);	
  return result;

};

double GradientReconstructionVectorForm_2::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>  **ext) const {
  return vector_form<double, double>(n, wt, u_ext, v, e, ext);
};

Ord GradientReconstructionVectorForm_2::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord>  **ext) const 
{
  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
};

VectorFormVol<double>* GradientReconstructionVectorForm_2::clone() const
{
  return new GradientReconstructionVectorForm_2(*this);
}


