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
    // Jacobian.
    add_matrix_form(new DefaultJacobianDiffusion(0, 0));

    // Residual.
    add_vector_form(new DefaultResidualDiffusion(0));
    add_vector_form(new DefaultVectorFormVol(0, HERMES_ANY, -const_f));
  };
};

