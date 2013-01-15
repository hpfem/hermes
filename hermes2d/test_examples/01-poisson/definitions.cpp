#include "definitions.h"

/* Weak forms */

CustomWeakFormPoisson::CustomWeakFormPoisson() : Hermes::Hermes2D::WeakForm<double>(1)
{
  // Jacobian forms.
  add_matrix_form(new Hermes::Hermes2D::ConstantWeakFormsH1::ConstantMatrixFormDx<double>(0, 0));
  add_matrix_form(new Hermes::Hermes2D::ConstantWeakFormsH1::ConstantMatrixFormDy<double>(0, 0));

  // Residual forms.
  add_vector_form(new Hermes::Hermes2D::ConstantWeakFormsH1::ConstantVectorFormVol<double>(0, Hermes::HERMES_ANY));
};

CustomWeakFormPoissonCacheCalculation::CustomWeakFormPoissonCacheCalculation() : Hermes::Hermes2D::WeakForm<double>(1)
{
  // Jacobian forms.
  add_matrix_form(new Hermes::Hermes2D::WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0));

  // Residual forms.
  add_vector_form(new Hermes::Hermes2D::WeakFormsH1::DefaultVectorFormVol<double>(0, Hermes::HERMES_ANY));
};