#include "definitions.h"

/* Weak forms */

CustomWeakFormPoisson::CustomWeakFormPoisson(std::string mat_al, HermesFunction<double>* lambda_al,
                                             std::string mat_cu, HermesFunction<double>* lambda_cu,
                                             HermesFunction<double>* src_term) : WeakForm<double>(1)
{
  // Jacobian forms.
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(0, 0, mat_al, lambda_al));
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(0, 0, mat_cu, lambda_cu));

  // Residual forms.
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(0, mat_al, lambda_al));
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(0, mat_cu, lambda_cu));
  add_vector_form(new WeakFormsH1::DefaultVectorFormVol<double>(0, HERMES_ANY, src_term));
};
