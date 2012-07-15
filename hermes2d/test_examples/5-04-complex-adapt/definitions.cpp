#include "definitions.h"

CustomWeakForm::CustomWeakForm(double mu_air, std::string mat_air,  
                               double mu_iron, double gamma_iron, std::string mat_iron, 
                               double mu_wire, std::string mat_wire, std::complex<double> j_ext, double omega) : Hermes::Hermes2D::WeakForm<std::complex<double> >(1)
{
  std::complex<double> ii =  std::complex<double>(0.0, 1.0);

  // Jacobian.
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<std::complex<double> >(0, 0, new Hermes1DFunction<std::complex<double> >(1.0/mu_air), mat_air));
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<std::complex<double> >(0, 0, new Hermes1DFunction<std::complex<double> >(1.0/mu_iron), mat_iron));
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<std::complex<double> >(0, 0, new Hermes1DFunction<std::complex<double> >(1.0/mu_wire), mat_wire));
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol<std::complex<double> >(0, 0, new Hermes2DFunction<std::complex<double> >(ii * omega * gamma_iron), mat_iron));

  // Residual.
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<std::complex<double> >(0, new Hermes1DFunction<std::complex<double> >(1.0/mu_air), mat_air));
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<std::complex<double> >(0, new Hermes1DFunction<std::complex<double> >(1.0/mu_iron), mat_iron));
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<std::complex<double> >(0, new Hermes1DFunction<std::complex<double> >(1.0/mu_wire), mat_wire));
  add_vector_form(new WeakFormsH1::DefaultVectorFormVol<std::complex<double> >(0, new Hermes2DFunction<std::complex<double> >(-j_ext), mat_wire));
  add_vector_form(new WeakFormsH1::DefaultResidualVol<std::complex<double> >(0, new Hermes2DFunction<std::complex<double> >(ii * omega * gamma_iron), mat_iron));
}
