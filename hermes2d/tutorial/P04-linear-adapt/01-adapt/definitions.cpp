#include "weakform/weakform.h"
#include "weakform_library/h1.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::VolumetricVectorForms;

/* Weak forms */

class CustomWeakFormPoisson : public WeakForm
{
public:
  CustomWeakFormPoisson(std::string mat_motor, double eps_motor, 
                        std::string mat_air, double eps_air) : WeakForm(1)
  {
    add_matrix_form(new DefaultLinearDiffusion(0, 0, mat_motor, eps_motor));
    add_matrix_form(new DefaultLinearDiffusion(0, 0, mat_air, eps_air));
    add_vector_form(new DefaultResidualLinearDiffusion(0, mat_motor, eps_motor));
    add_vector_form(new DefaultResidualLinearDiffusion(0, mat_air, eps_air));
  };
};
