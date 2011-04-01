#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/h1.h"

using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::VolumetricVectorForms;

class CustomWeakFormPoisson : public WeakForm
{
public:
  CustomWeakFormPoisson(double const_f) : WeakForm(1)
  {
    add_matrix_form(new DefaultLinearDiffusion(0, 0));
    add_vector_form(new DefaultVectorFormConst(0, const_f));
  };
};
