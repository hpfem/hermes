#include "weakform/weakform.h"
#include "weakform_library/h1.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"

using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::VolumetricVectorForms;

/* Weak forms */

class WeakFormPoisson : public WeakForm
{
public:
  WeakFormPoisson(double const_f) : WeakForm(1) {
    add_matrix_form(new DefaultLinearDiffusion(0, 0));
    add_vector_form(new DefaultVectorFormConst(0, const_f));
  };
};

