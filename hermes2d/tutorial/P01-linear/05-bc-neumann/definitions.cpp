#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/laplace.h"

class CustomWeakFormPoissonNeumann : public WeakForm
{
public:
  CustomWeakFormPoissonNeumann(double const_f,
                               double const_gamma_bottom,
                               double const_gamma_outer,
                               double const_gamma_left) : WeakForm(1)
  {
    add_matrix_form(new DefaultMatrixFormVolConst(0, 0));
    add_vector_form(new DefaultVectorFormVolConst(0, const_f));

    add_vector_form_surf(new DefaultVectorFormSurfConst(0, BDY_BOTTOM, const_gamma_bottom));
    add_vector_form_surf(new DefaultVectorFormSurfConst(0, BDY_OUTER, const_gamma_outer));
    add_vector_form_surf(new DefaultVectorFormSurfConst(0, BDY_LEFT, const_gamma_left));
  };
};
