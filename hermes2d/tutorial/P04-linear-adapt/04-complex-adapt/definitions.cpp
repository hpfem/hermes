#include "weakform/weakform.h"
#include "weakform_library/h1.h"

using namespace WeakFormsH1;
  
/* Weak forms */
class CustomWeakFormMagnetics : public WeakForm
{ 
public:
  CustomWeakFormMagnetics(std::string mat_air,  double mu_air,
                          std::string mat_iron, double mu_iron, double gamma_iron,
                          std::string mat_wire, double mu_wire, scalar j_ext, double omega)
  : WeakForm(1) {
    scalar ii =  cplx(0.0, 1.0);

    // Jacobian.
    add_matrix_form(new DefaultJacobianDiffusion(0, 0, mat_air,  1.0/mu_air));
    add_matrix_form(new DefaultJacobianDiffusion(0, 0, mat_iron, 1.0/mu_iron));
    add_matrix_form(new DefaultJacobianDiffusion(0, 0, mat_wire, 1.0/mu_wire));
    add_matrix_form(new DefaultMatrixFormVol(0, 0, mat_iron, ii * omega * gamma_iron));

    // Residual.
    add_vector_form(new DefaultResidualDiffusion(0, mat_air, 1.0/mu_air));
    add_vector_form(new DefaultResidualDiffusion(0, mat_iron, 1.0/mu_iron));
    add_vector_form(new DefaultResidualDiffusion(0, mat_wire, 1.0/mu_wire));
    add_vector_form(new DefaultVectorFormVol(0, mat_wire, -j_ext));
    add_vector_form(new DefaultResidualVol(0, mat_iron, ii * omega * gamma_iron));
  };
};
