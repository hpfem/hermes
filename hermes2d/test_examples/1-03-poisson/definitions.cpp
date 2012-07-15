#include "definitions.h"

/* Weak forms */

CustomWeakFormPoisson::CustomWeakFormPoisson(Hermes::Hermes1DFunction<double>* lambda_al, std::string mat_al, 
                                             Hermes::Hermes1DFunction<double>* lambda_cu, std::string mat_cu, 
                                             Hermes::Hermes2DFunction<double>* src_term) : Hermes::Hermes2D::WeakForm<double>(1)
{
  // Jacobian forms.
  add_matrix_form(new Hermes::Hermes2D::WeakFormsH1::DefaultJacobianDiffusion<double>(0, 0, lambda_al, mat_al));
  add_matrix_form(new Hermes::Hermes2D::WeakFormsH1::DefaultJacobianDiffusion<double>(0, 0, lambda_cu, mat_cu));

  // Residual forms.
  add_vector_form(new Hermes::Hermes2D::WeakFormsH1::DefaultResidualDiffusion<double>(0, lambda_al, mat_al));
  add_vector_form(new Hermes::Hermes2D::WeakFormsH1::DefaultResidualDiffusion<double>(0, lambda_cu, mat_cu));
  add_vector_form(new Hermes::Hermes2D::WeakFormsH1::DefaultVectorFormVol<double>(0, src_term));
};
