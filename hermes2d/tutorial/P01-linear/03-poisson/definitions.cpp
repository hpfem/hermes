#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/laplace.h"

using namespace Laplace::VolumetricMatrixForms;
using namespace Laplace::VolumetricVectorForms;

class CustomWeakFormPoisson : public WeakForm
{
public:
  CustomWeakFormPoisson(double const_f) : WeakForm(1)
  {
    add_matrix_form(new DefaultMatrixFormStiffness(0, 0));
    add_vector_form(new DefaultVectorFormConst(0, const_f));
  };
};
