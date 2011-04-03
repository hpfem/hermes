#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/h1.h"

using namespace WeakFormsH1;

/* Weak forms */

class CustomWeakFormPoissonNewton : public WeakForm
{
public:
  CustomWeakFormPoissonNewton(double h, double T0, double lambda, std::string bdy_newton) : WeakForm(1)
  {
    add_matrix_form(new VolumetricMatrixForms::DefaultLinearDiffusion(0, 0, lambda, HERMES_SYM, HERMES_AXISYM_Y));
    add_matrix_form_surf(new SurfaceMatrixForms::DefaultMatrixFormSurf(0, 0, bdy_newton, lambda * h, HERMES_AXISYM_Y));
    add_vector_form_surf(new SurfaceVectorForms::DefaultVectorFormSurf(0, bdy_newton, lambda * h * T0, HERMES_AXISYM_Y));
  };
};

