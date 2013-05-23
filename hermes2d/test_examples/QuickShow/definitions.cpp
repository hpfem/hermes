#include "definitions.h"

CustomWeakFormSteadyState::CustomWeakFormSteadyState(Hermes1DFunction<double>* thermal_conductivity,
                                                     Hermes2DFunction<double>* heat_source)
                                                     : WeakForm<double>(1)
{
#ifdef LINEAR_NONLINEAR_SWITCH
  
  // Matrix.
  this->add_matrix_form(new DefaultMatrixFormDiffusion<double>(0, 0, HERMES_ANY_MARKER, thermal_conductivity));

  // Rhs
  add_vector_form(new DefaultVectorFormVol<double>(0, HERMES_ANY_MARKER, heat_source));

#else

  // Jacobian.
  this->add_matrix_form(new DefaultJacobianDiffusion<double>(0, 0, HERMES_ANY_MARKER, thermal_conductivity));

  // Residual.
  add_vector_form(new DefaultResidualDiffusion<double>(0, HERMES_ANY_MARKER, thermal_conductivity));
  heat_source->const_value *= -1;
  add_vector_form(new DefaultVectorFormVol<double>(0, HERMES_ANY_MARKER, heat_source));

#endif
};



CustomWeakFormTimeDependent::CustomWeakFormTimeDependent(Hermes1DFunction<double>* thermal_conductivity,
                                                         Hermes2DFunction<double>* heat_source,
                                                         MeshFunctionSharedPtr<double> prev_sln)
                                                         : CustomWeakFormSteadyState(thermal_conductivity, heat_source)
{
  // Jacobian.
  this->add_matrix_form(new CustomMatrixFormVol(0, 0));

  // Residual.
  add_vector_form(new CustomResidualFormVol(0));
  
  CustomVectorFormVol* vec_form_vol = new CustomVectorFormVol(0);
  vec_form_vol->set_ext(prev_sln);
  add_vector_form(vec_form_vol);
};



CustomNonlinearity::CustomNonlinearity(double alpha): Hermes1DFunction<double>()
{
  this->is_const = false;
  this->alpha = alpha;
}

double CustomNonlinearity::value(double u) const
{
  return 1 + std::pow(u, alpha);
}

Ord CustomNonlinearity::value(Ord u) const
{ 
  return Ord(10);
}

double CustomNonlinearity::derivative(double u) const
{
  return alpha * std::pow(u, alpha - 1.0);
}

Ord CustomNonlinearity::derivative(Ord u) const
{
  return Ord(10);
}

CustomEssentialBCNonConst::CustomEssentialBCNonConst(std::string marker) : EssentialBoundaryCondition<double>(Hermes::vector<std::string>())
{
  markers.push_back(marker);
}

EssentialBoundaryCondition<double>::EssentialBCValueType CustomEssentialBCNonConst::get_value_type() const 
{ 
  return EssentialBoundaryCondition<double>::BC_FUNCTION; 
}

double CustomEssentialBCNonConst::value(double x, double y, double n_x, double n_y, double t_x, double t_y) const
{
  return (x + 10) * (y + 10) / 100.;
}

void CustomInitialCondition::derivatives (double x, double y, double& dx, double& dy) const 
{
  dx = (y + 10)/100.;
  dy = (x + 10)/100.;
};

double CustomInitialCondition::value (double x, double y) const 
{
  return (x + 10) * (y + 10) / 100.;
}

Ord CustomInitialCondition::ord(Ord x, Ord y) const 
{
  return x*y;
}
  
MeshFunction<double>* CustomInitialCondition::clone() const
{
  if(this->get_type() == HERMES_SLN)
    return Solution<double>::clone();
  else return new CustomInitialCondition(mesh);
}

CustomInitialCondition::~CustomInitialCondition()
{
}

CustomWeakFormTimeDependent::CustomMatrixFormVol::CustomMatrixFormVol(int i, int j) : DefaultMatrixFormVol<double>(i, j, HERMES_ANY_MARKER) 
{
}

double CustomWeakFormTimeDependent::CustomMatrixFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double> **ext) const
{
  return DefaultMatrixFormVol<double>::value(n, wt, u_ext, u, v, e, ext) / this->wf->get_current_time_step();
}

MatrixFormVol<double>* CustomWeakFormTimeDependent::CustomMatrixFormVol::clone() const
{
  return new CustomWeakFormTimeDependent::CustomMatrixFormVol(*this);
}

CustomWeakFormTimeDependent::CustomVectorFormVol::CustomVectorFormVol(int i) : DefaultResidualVol<double>(i, HERMES_ANY_MARKER)
{
}

double CustomWeakFormTimeDependent::CustomVectorFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const
{
  double result = 0;
  
  for (int i = 0; i < n; i++)
    result += wt[i] * ext[0]->val[i] * v->val[i];
         
  return -result / this->wf->get_current_time_step();
}

VectorFormVol<double>* CustomWeakFormTimeDependent::CustomVectorFormVol::clone() const
{
  return new CustomWeakFormTimeDependent::CustomVectorFormVol(*this);
}

CustomWeakFormTimeDependent::CustomResidualFormVol::CustomResidualFormVol(int i) : DefaultResidualVol<double>(i, HERMES_ANY_MARKER)
{
}

double CustomWeakFormTimeDependent::CustomResidualFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double> **ext) const
{
  double result = 0;
  
  for (int i = 0; i < n; i++)
    result += wt[i] * u_ext[0]->val[i] * v->val[i];
         
  return result / this->wf->get_current_time_step();
}

VectorFormVol<double>* CustomWeakFormTimeDependent::CustomResidualFormVol::clone() const
{
  return new CustomWeakFormTimeDependent::CustomResidualFormVol(*this);
}