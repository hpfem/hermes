#include "weakform/weakform.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/h1.h"

using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::VolumetricVectorForms;

/* Weak forms */

class CustomWeakFormPoisson : public WeakForm
{
public:
  CustomWeakFormPoisson(std::string mat_al, double lambda_al, 
                        std::string mat_cu, double lambda_cu, 
                        double vol_heat_src) : WeakForm(1)
  {
    // Jacobian forms - volumetric.
    add_matrix_form(new DefaultLinearDiffusion(0, 0, mat_al, lambda_al));
    add_matrix_form(new DefaultLinearDiffusion(0, 0, mat_cu, lambda_cu));

    // Residual forms - volumetric.
    add_vector_form(new DefaultResidualLinearDiffusion(0, mat_al, lambda_al));
    add_vector_form(new DefaultResidualLinearDiffusion(0, mat_cu, lambda_cu));
    add_vector_form(new DefaultVectorFormConst(0, -vol_heat_src));
  };
};
