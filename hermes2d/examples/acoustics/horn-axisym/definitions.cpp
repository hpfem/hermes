#include "weakform/weakform.h"
#include "weakform_library/h1.h"

using namespace WeakFormsH1;

/* Weak forms */

class CustomWeakFormAcoustics : public WeakForm
{
public:
  CustomWeakFormAcoustics(std::string bdy_newton, double rho,
                          double sound_speed, double omega)
  : WeakForm(1) {
    scalar ii =  cplx(0.0, 1.0);

    // Jacobian.
    add_matrix_form(new DefaultJacobianDiffusion(0, 0, HERMES_ANY, 1.0/rho, HERMES_DEFAULT_FUNCTION, HERMES_SYM));
    add_matrix_form(new DefaultMatrixFormVol(0, 0, HERMES_ANY, - sqr(omega) / rho / sqr(sound_speed), HERMES_DEFAULT_FUNCTION, HERMES_SYM));
    add_matrix_form_surf(new DefaultMatrixFormSurf(0, 0, bdy_newton, -ii * omega / rho / sound_speed));

    // Residual.
    add_vector_form(new DefaultResidualDiffusion(0, HERMES_ANY, 1.0/rho));
    add_vector_form(new DefaultResidualVol(0, HERMES_ANY, - sqr(omega) / rho / sqr(sound_speed)));
    add_vector_form_surf(new DefaultResidualSurf(0, bdy_newton, -ii * omega / rho / sound_speed));
  };
};
