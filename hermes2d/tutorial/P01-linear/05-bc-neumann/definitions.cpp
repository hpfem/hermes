#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform/sample_weak_forms.h"

class WeakFormNeumann : public WeakFormPoisson
{
public:
  WeakFormNeumann(double const_f,
                  double const_gamma_bottom,
                  double const_gamma_outer,
                  double const_gamma_left) : WeakFormPoisson(1)
  {
    add_matrix_form(new WeakFormPoisson::MatrixFormVol(0, 0));
    add_vector_form(new WeakFormPoisson::VectorFormVol(0, const_f));

    add_vector_form_surf(new WeakFormPoisson::VectorFormSurfNeumann(0, BDY_BOTTOM, const_gamma_bottom));
    add_vector_form_surf(new WeakFormPoisson::VectorFormSurfNeumann(0, BDY_OUTER, const_gamma_outer));
    add_vector_form_surf(new WeakFormPoisson::VectorFormSurfNeumann(0, BDY_LEFT, const_gamma_left));
  };
};
