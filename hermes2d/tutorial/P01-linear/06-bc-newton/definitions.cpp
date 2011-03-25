#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/laplace.h"

class CustomWeakFormPoissonNewton : public WeakForm
{
public:
  // Problem parameters.
  CustomWeakFormPoissonNewton(double h, double T0, std::string natural_bc_bnd_part) : WeakForm(1)
  {
    add_matrix_form(new DefaultMatrixFormStiffness(0, 0));
    add_matrix_form_surf(new DefaultMatrixFormSurf(0, 0, natural_bc_bnd_part, h));
    add_vector_form_surf(new DefaultVectorFormSurf(0, natural_bc_bnd_part, h * T0));
  };
};

