#include "weakform/weakform.h"
#include "weakform_library/maxwell.h"
#include "weakform_library/h1.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

using namespace WeakFormsMaxwell::VolumetricMatrixForms;
using namespace WeakFormsMaxwell::VolumetricVectorForms;
using namespace WeakFormsH1::VolumetricVectorForms;

/* Weak forms */

class CustomWeakFormMagnetostatics : public WeakForm
{
public:
  CustomWeakFormMagnetostatics(std::string material_iron_1, std::string material_iron_2, 
                               CubicSpline* mu_inv_iron, std::string material_air,
                               std::string material_copper, double mu_vacuum, 
                               double current_density, int order_inc = 3) 
    : WeakForm(1) {

    // Jacobian.
    add_matrix_form(new DefaultLinearMagnetostatics(0, 0, material_air, 1.0, HERMES_SYM, HERMES_AXISYM_Y, order_inc));
    add_matrix_form(new DefaultLinearMagnetostatics(0, 0, material_copper, 1.0, HERMES_SYM, HERMES_AXISYM_Y, order_inc));
    add_matrix_form(new DefaultJacobianNonlinearMagnetostatics(0, 0, material_iron_1, mu_inv_iron, 1.0, HERMES_SYM, HERMES_AXISYM_Y, order_inc));
    add_matrix_form(new DefaultJacobianNonlinearMagnetostatics(0, 0, material_iron_2, mu_inv_iron, 1.0, HERMES_SYM, HERMES_AXISYM_Y, order_inc));

    // Residual.
    add_vector_form(new DefaultResidualLinearMagnetostatics(0, material_air, 1.0, HERMES_AXISYM_Y, order_inc));
    add_vector_form(new DefaultResidualLinearMagnetostatics(0, material_copper, 1.0, HERMES_AXISYM_Y, order_inc));
    add_vector_form(new DefaultResidualNonlinearMagnetostatics(0, material_iron_1, mu_inv_iron, 1.0, HERMES_AXISYM_Y, order_inc));
    add_vector_form(new DefaultResidualNonlinearMagnetostatics(0, material_iron_2, mu_inv_iron, 1.0, HERMES_AXISYM_Y, order_inc));
    add_vector_form(new DefaultVectorFormConst(0, material_copper, -current_density * mu_vacuum, HERMES_AXISYM_Y));
  };
};



