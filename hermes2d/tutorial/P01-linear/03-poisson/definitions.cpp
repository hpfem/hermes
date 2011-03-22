#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform/sample_weak_forms.h"

class MyWeakFormPoisson : public WeakFormLaplace
{
public:
  MyWeakFormPoisson(double const_f) : WeakFormLaplace() {
    VectorFormConstant* my_vector_form = new VectorFormConstant(0);
    my_vector_form->param.push_back(const_f);
    add_vector_form(my_vector_form);
  };
};

