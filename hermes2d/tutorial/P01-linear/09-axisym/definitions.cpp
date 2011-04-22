#include "weakform/weakform.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/h1.h"

using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::VolumetricVectorForms;
using namespace WeakFormsH1::SurfaceMatrixForms;
using namespace WeakFormsH1::SurfaceVectorForms;

/* Weak forms */

class CustomWeakFormPoissonNewton : public WeakForm
{
public:
  CustomWeakFormPoissonNewton(double lambda, double alpha, double T0, std::string bdy_heat_flux)
        : WeakForm(1)
  {
    // Jacobian form - volumetric.
    add_matrix_form(new DefaultLinearDiffusion(0, 0, HERMES_ANY, lambda, HERMES_SYM, HERMES_AXISYM_Y));

    // Jacobian form - surface.
    add_matrix_form_surf(new DefaultMatrixFormSurf(0, 0, bdy_heat_flux, alpha, HERMES_AXISYM_Y));

    // Residual forms - volumetric.
    add_vector_form(new DefaultResidualLinearDiffusion(0, HERMES_ANY, lambda, HERMES_AXISYM_Y));

    // Residual form - surface.
    add_vector_form_surf(new DefaultResidualSurfConst(0, bdy_heat_flux, alpha, HERMES_AXISYM_Y));
    add_vector_form_surf(new DefaultVectorFormSurf(0, bdy_heat_flux, -alpha * T0, HERMES_AXISYM_Y));
  };
};

