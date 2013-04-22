#include "definitions.h"

static double custom_fn(int n, double *wt, Func<double> *u, Func<double> *v)
{
  double result = double(0);
  for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * v->val[i];
  return result;
}

static double custom_fn(int n, double *wt, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v)
{
  double result = double(0);
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] - u->val_neighbor[i]) * (v->val[i] - v->val_neighbor[i]);
  return result;
}

CustomNormFormVol::CustomNormFormVol(int i, int j) : NormFormVol<double>(i, j)
{
  this->set_area(HERMES_ANY);
}

CustomNormFormSurf::CustomNormFormSurf(int i, int j) : NormFormSurf<double>(i, j)
{
  this->set_area(HERMES_ANY);
}

CustomNormFormDG::CustomNormFormDG(int i, int j) : NormFormDG<double>(i, j)
{
}

double CustomNormFormVol::value(int n, double *wt, Func<double> *u, Func<double> *v, Geom<double> *e) const
{
  return custom_fn(n, wt, u, v);
}

double CustomNormFormSurf::value(int n, double *wt, Func<double> *u, Func<double> *v, Geom<double> *e) const
{
  return custom_fn(n, wt, u, v);
}

double CustomNormFormDG::value(int n, double *wt, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v, Geom<double> *e) const
{
  return custom_fn(n, wt, u, v);
}


CustomExactSolutionScalar::CustomExactSolutionScalar(MeshSharedPtr mesh)
      : ExactSolutionScalar<double>(mesh)
{
}

double CustomExactSolutionScalar::value(double x, double y) const 
{
  if(this->element_id == 1 || this->element_id == 3)
    return 0.0;
  else
    return 2.0;
}

void CustomExactSolutionScalar::derivatives(double x, double y, double& dx, double& dy) const 
{
  dx = 0.0;
  dy = 0.0;
}

Ord CustomExactSolutionScalar::ord(Ord x, Ord y) const 
{
  return Ord(1);
}

CustomExactSolutionScalar::~CustomExactSolutionScalar() 
{
}

MeshFunction<double>* CustomExactSolutionScalar::clone() const
{
  return new CustomExactSolutionScalar(this->mesh);
}

void CustomExactSolutionScalar::set_active_element(Element* e)
{
  Solution<double>::set_active_element(e);
  this->element_id = e->id;
}