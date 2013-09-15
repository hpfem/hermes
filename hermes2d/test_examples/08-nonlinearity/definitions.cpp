#include "definitions.h"

CustomNonlinearity::CustomNonlinearity(double alpha): Hermes1DFunction<double>()
{
  this->is_const = false;
  this->alpha = alpha;
}

double CustomNonlinearity::value(double u) const
{
  return 1 + Hermes::pow(u, alpha);
}

Ord CustomNonlinearity::value(Ord u) const
{
  return Ord(10);
}

double CustomNonlinearity::derivative(double u) const
{
  return alpha * Hermes::pow(u, alpha - 1.0);
}

Ord CustomNonlinearity::derivative(Ord u) const
{
  // Same comment as above applies.
  return Ord(10);
}

double CustomInitialCondition::value(double x, double y) const 
{
  return (x+10) * (y+10) / 100. + 2;
}

void CustomInitialCondition::derivatives(double x, double y, double& dx, double& dy) const 
{   
  dx = (y+10) / 100.;
  dy = (x+10) / 100.;
}

Ord CustomInitialCondition::ord(double x, double y) const 
{
  return Hermes::Ord(x*y);
}
  
MeshFunction<double>* CustomInitialCondition::clone() const
{
  return new CustomInitialCondition(this->mesh);
}

EssentialBoundaryCondition<double>::EssentialBCValueType CustomEssentialBCNonConst::get_value_type() const 
{ 
  return EssentialBoundaryCondition<double>::BC_FUNCTION; 
}

double CustomEssentialBCNonConst::value(double x, double y, double n_x, double n_y, 
                                        double t_x, double t_y) const
{
  return (x+10) * (y+10) / 100.;
}
