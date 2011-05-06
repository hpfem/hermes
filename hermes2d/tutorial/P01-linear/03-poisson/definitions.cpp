#include "definitions.h"

/* Weak forms */

CustomWeakFormPoisson::CustomWeakFormPoisson(std::string mat_al, double lambda_al,
                                             std::string mat_cu, double lambda_cu,
                                             double vol_heat_src) : WeakForm(1)
{
  // Jacobian forms - volumetric.
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(0, 0, mat_al, lambda_al));
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(0, 0, mat_cu, lambda_cu));

  // Residual forms - volumetric.
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion(0, mat_al, lambda_al));
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion(0, mat_cu, lambda_cu));
  add_vector_form(new WeakFormsH1::DefaultVectorFormVol(0, HERMES_ANY, -vol_heat_src));
};
