#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform/sample_weak_forms.h"

class CustomWeakFormPoisson : public WeakFormPoisson
{
public:
  CustomWeakFormPoisson(double const_f) : WeakFormPoisson(1)
  {
    add_matrix_form(new WeakFormPoisson::MatrixFormVol(0, 0));
    add_vector_form(new WeakFormPoisson::VectorFormVol(0, const_f));
  };
};

