#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::WeakFormsH1;

CustomWeakFormPoisson::CustomWeakFormPoisson() : Hermes::Hermes2D::WeakForm<double>(3)
{
  this->add_matrix_form(new DefaultJacobianDiffusion<double>(0, 0));
  this->add_matrix_form(new DefaultJacobianDiffusion<double>(1, 1));
  this->add_matrix_form(new DefaultJacobianDiffusion<double>(2, 2));

  // Residual.
  this->add_vector_form(new DefaultResidualDiffusion<double>(0));
  this->add_vector_form(new DefaultVectorFormVol<double>(0));
  this->add_vector_form(new DefaultResidualDiffusion<double>(1));
  this->add_vector_form(new DefaultVectorFormVol<double>(1));
  this->add_vector_form(new CustomVectorFormVol(2));
  this->add_vector_form(new DefaultResidualDiffusion<double>(2));
  this->add_vector_form(new DefaultVectorFormVol<double>(2));
};

CustomVectorFormVol::CustomVectorFormVol(int i)
  : VectorFormVol<double>(i)
{}

double CustomVectorFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
  Geom<double> *e, ExtData<double> *ext) const 
{
  double result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * u_ext[0]->val[i] * v->val[i];
  return result;
}

Ord CustomVectorFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
  Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
    result += wt[i] * u_ext[0]->val[i] * v->val[i];
  return result;
}