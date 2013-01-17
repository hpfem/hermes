#include "definitions.h"

CustomNonlinearity::CustomNonlinearity(double alpha): Hermes1DFunction<double>()
{
  this->is_const = false;
  this->alpha = alpha;
}

double CustomNonlinearity::value(double u) const
{
  return -1 - std::pow(u, alpha);
}

Ord CustomNonlinearity::value(Ord u) const
{ 
  return Ord(10);
}

double CustomNonlinearity::derivative(double u) const
{
  return -alpha * std::pow(u, alpha - 1.0);
}

Ord CustomNonlinearity::derivative(Ord u) const
{
  return Ord(10);
}

EssentialBCNonConst::EssentialBCNonConst(std::string marker) : EssentialBoundaryCondition<double>(Hermes::vector<std::string>())
{
  markers.push_back(marker);
}

EssentialBoundaryCondition<double>::EssentialBCValueType EssentialBCNonConst::get_value_type() const 
{ 
  return EssentialBoundaryCondition<double>::BC_FUNCTION; 
}

double EssentialBCNonConst::value(double x, double y, double n_x, double n_y, double t_x, double t_y) const
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
  return (x + 10) * (y + 10) / 100.;
}
  
MeshFunction<double>* CustomInitialCondition::clone() const
{
  Mesh* new_mesh = new Mesh();
  new_mesh->copy(this->mesh);
  return new CustomInitialCondition(new_mesh);
}

CustomInitialCondition::~CustomInitialCondition()
{
  delete mesh;
}

CustomWeakFormPoisson::CustomWeakFormPoisson(Hermes1DFunction<double>* coeff, Hermes2DFunction<double>* f)
{
  if(1)
  {
    // Jacobian.
    NonlinearElemwiseDiffusionDerivativeForm<double>* form_derivative = new NonlinearElemwiseDiffusionDerivativeForm<double>(0, 0);
    NonlinearElemwiseDiffusionValueForm<double>* form_value = new NonlinearElemwiseDiffusionValueForm<double>(0, 0);

    form_derivative->set_elemwise_parameter(new ElemwiseParameterNonlinearHermesFunc<double>(coeff, ElemwiseParameterNonlinearDerivative));
    form_value->set_elemwise_parameter(new ElemwiseParameterNonlinearHermesFunc<double>(coeff, ElemwiseParameterNonlinearValue));

    this->add_matrix_form(form_derivative);
    this->add_matrix_form(form_value);
  }

  if(0)
  {
    this->add_matrix_form(new DefaultJacobianDiffusion<double>(0, 0, HERMES_ANY, coeff));
  }

  // Residual.
  this->add_vector_form(new DefaultResidualDiffusion<double>(0, HERMES_ANY, coeff));
  this->add_vector_form(new DefaultVectorFormVol<double>(0, HERMES_ANY, f));
}