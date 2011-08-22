#include "definitions.h"

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
  this->add_vector_form(new DefaultResidualDiffusion<double>(2));
  this->add_vector_form(new DefaultVectorFormVol<double>(2));
};
