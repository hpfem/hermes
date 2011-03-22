#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform/sample_weak_forms.h"

class MyWeakFormPoisson : public WeakFormLaplace
{
public:
  MyWeakFormPoisson(double const_f) : WeakFormLaplace() {
    add_vector_form(new VectorFormConstant(0, const_f));
  };
private:
  double const_f;
};

