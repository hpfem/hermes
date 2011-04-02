#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/h1.h"

using namespace WeakFormsH1;

/* Weak forms */

class CustomWeakFormPoissonNewton : public WeakForm
{
public:
  CustomWeakFormPoissonNewton(double h, double T0, std::string natural_bc_bnd_part) : WeakForm(1)
  {
    add_matrix_form(new VolumetricMatrixForms::DefaultLinearDiffusion(0, 0, 1.0, HERMES_SYM, HERMES_AXISYM_Y));
    add_matrix_form_surf(new SurfaceMatrixForms::DefaultMatrixFormSurf(0, 0, natural_bc_bnd_part, h, HERMES_AXISYM_Y));
    add_vector_form_surf(new SurfaceVectorForms::DefaultVectorFormSurf(0, natural_bc_bnd_part, h * T0, HERMES_AXISYM_Y));
  };
};

