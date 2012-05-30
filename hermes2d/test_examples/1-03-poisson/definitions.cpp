#include "definitions.h"

/* Weak forms */

CustomWeakFormPoisson::CustomWeakFormPoisson(std::string mat_al, Hermes::Hermes1DFunction<double>* lambda_al,
                                             std::string mat_cu, Hermes::Hermes1DFunction<double>* lambda_cu,
                                             Hermes::Hermes2DFunction<double>* src_term) : Hermes::Hermes2D::WeakForm<double>(1)
{
  // Jacobian forms.
  add_matrix_form(new Hermes::Hermes2D::WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0, mat_al, lambda_al));
  add_matrix_form(new Hermes::Hermes2D::WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0, mat_cu, lambda_cu));

  // Residual forms.
  add_vector_form(new Hermes::Hermes2D::WeakFormsH1::DefaultVectorFormVol<double>(0, Hermes::HERMES_ANY, src_term));
};
