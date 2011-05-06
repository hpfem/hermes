#include "definitions.h"

/* Weak forms */

CustomWeakFormPoissonDirichlet::CustomWeakFormPoissonDirichlet(std::string mat_al, double lambda_al,
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

/* Custom non-constant Dirichlet condition */

CustomDirichletCondition::CustomDirichletCondition(Hermes::vector<std::string> markers, 
                                                   double A, double B, double C)
  : EssentialBoundaryCondition(markers), A(A), B(B), C(C) { }

EssentialBoundaryCondition::EssentialBCValueType CustomDirichletCondition::get_value_type() const
{ 
  return EssentialBoundaryCondition::BC_FUNCTION; 
}

scalar CustomDirichletCondition::value(double x, double y, double n_x, double n_y, 
                                       double t_x, double t_y) const 
{
  return A*x + B*y + C;
}
