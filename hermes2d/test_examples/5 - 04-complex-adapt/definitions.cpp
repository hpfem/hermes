#include "definitions.h"
  
CustomWeakForm::CustomWeakForm(std::string mat_air,  double mu_air,
                               std::string mat_iron, double mu_iron, double gamma_iron,
                               std::string mat_wire, double mu_wire, scalar j_ext, double omega) : WeakForm(1) 
{
  scalar ii =  cplx(0.0, 1.0);

  // Jacobian.
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(0, 0, mat_air,  new HermesFunction(1.0/mu_air)));
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(0, 0, mat_iron, new HermesFunction(1.0/mu_iron)));
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(0, 0, mat_wire, new HermesFunction(1.0/mu_wire)));
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol(0, 0, mat_iron, new HermesFunction(ii * omega * gamma_iron)));

  // Residual.
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion(0, mat_air, new HermesFunction(1.0/mu_air)));
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion(0, mat_iron, new HermesFunction(1.0/mu_iron)));
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion(0, mat_wire, new HermesFunction(1.0/mu_wire)));
  add_vector_form(new WeakFormsH1::DefaultVectorFormVol(0, mat_wire, new HermesFunction(-j_ext)));
  add_vector_form(new WeakFormsH1::DefaultResidualVol(0, mat_iron, new HermesFunction(ii * omega * gamma_iron)));
}
