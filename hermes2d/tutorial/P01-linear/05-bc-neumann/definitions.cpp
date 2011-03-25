#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/laplace.h"

class CustomWeakFormPoissonNeumann : public WeakForm
{
public:
  CustomWeakFormPoissonNeumann(double const_f,
                               std::string bdy_bottom,
                               double const_gamma_bottom,
                               std::string bdy_outer,
                               double const_gamma_outer,
                               std::string bdy_left,
                               double const_gamma_left) : WeakForm(1)
  {
    add_matrix_form(new DefaultMatrixFormStiffness(0, 0));
    add_vector_form(new DefaultVectorFormVolConst(0, const_f));

    add_vector_form_surf(new DefaultVectorFormSurf(0, bdy_bottom, const_gamma_bottom));
    add_vector_form_surf(new DefaultVectorFormSurf(0, bdy_outer, const_gamma_outer));
    add_vector_form_surf(new DefaultVectorFormSurf(0, bdy_left, const_gamma_left));
  };
};
