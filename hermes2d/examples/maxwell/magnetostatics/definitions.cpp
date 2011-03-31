#include "weakform/weakform.h"
#include "weakform_library/h1.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

using namespace WeakFormsH1;
using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::VolumetricVectorForms;

/* Weak forms */

class CustomWeakFormMagnetostatics : public WeakForm
{
public:
  CustomWeakFormMagnetostatics(CubicSpline* mu_iron, std::string material_air, double mu_air, 
                               std::string material_copper, double mu_copper, double current_density) 
    : WeakForm(1) {
    // Jacobian.
    add_matrix_form(new DefaultLinearMagnetostatics(0, 0, material_air, mu_air));
    add_matrix_form(new DefaultLinearMagnetostatics(0, 0, material_copper, mu_copper));
    add_matrix_form(new DefaultJacobianNonlinearMagnetostatics(0, 0, std::string("1"), mu_iron));
    add_matrix_form(new DefaultJacobianNonlinearMagnetostatics(0, 0, std::string("2"), mu_iron));           //FIXME
    // Residual.
    add_vector_form(new DefaultResidualLinearMagnetostatics(0, material_air, mu_air));
    add_vector_form(new DefaultResidualLinearMagnetostatics(0, material_copper, mu_copper));
    add_vector_form(new DefaultResidualNonlinearMagnetostatics(0, std::string("1"), mu_iron));
    add_vector_form(new DefaultResidualNonlinearMagnetostatics(0, std::string("2"), mu_iron));              //FIXME
    add_vector_form(new DefaultVectorFormConst(0, material_copper, -current_density));
  };
};



