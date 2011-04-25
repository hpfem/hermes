#include "weakform/weakform.h"
#include "weakform_library/h1.h"
#include "weakform_library/elasticity.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"

using namespace std;

class CustomWeakFormLinearElasticity : public WeakForm
{
public:
  CustomWeakFormLinearElasticity(double E, double nu, double rho_g,
                                 std::string surface_force_bdy, double f0, double f1)
            : WeakForm(2)
  {
    double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
    double mu = E / (2*(1 + nu));

    // Jacobian.
    add_matrix_form(new WeakFormsElasticity::DefaultJacobianElasticity_0_0(0, 0, lambda, mu));
    add_matrix_form(new WeakFormsElasticity::DefaultJacobianElasticity_0_1(0, 1, lambda, mu));
    add_matrix_form(new WeakFormsElasticity::DefaultJacobianElasticity_1_1(1, 1, lambda, mu));

    // Residual - first equation.
    add_vector_form(new WeakFormsElasticity::DefaultResidualElasticity_0_0(0, HERMES_ANY, lambda, mu));
    add_vector_form(new WeakFormsElasticity::DefaultResidualElasticity_0_1(0, HERMES_ANY, lambda, mu));
    // Surface force (first component).
    add_vector_form_surf(new WeakFormsH1::SurfaceVectorForms::DefaultVectorFormSurf(0, surface_force_bdy, -f0)); 

    // Residual - second equation.
    add_vector_form(new WeakFormsElasticity::DefaultResidualElasticity_1_0(1, HERMES_ANY, lambda, mu));
    add_vector_form(new WeakFormsElasticity::DefaultResidualElasticity_1_1(1, HERMES_ANY, lambda, mu));
    // Gravity loading in the second vector component.
    add_vector_form(new WeakFormsH1::VolumetricVectorForms::DefaultVectorFormVol(1, HERMES_ANY, -rho_g));
    // Surface force (second component).
    add_vector_form_surf(new WeakFormsH1::SurfaceVectorForms::DefaultVectorFormSurf(1, surface_force_bdy, -f1)); 


    /*

    DOES NOT WORK YET

    // Jacobian matrix.
    // There is one multi-component and one single-component form since we want to exploit symmetry of the forms.
    add_multicomponent_matrix_form(new WeakFormsElasticity::DefaultMultiComponentVolumetricMatrixFormLinearSym(
         Hermes::vector<std::pair<unsigned int, unsigned int> >(make_pair(0, 0), make_pair(1, 1)), HERMES_ANY, lambda, mu));

    add_matrix_form(new WeakFormsElasticity::DefaultVolumetricMatrixFormLinear_x_y(0, 1, lambda, mu));

    // Residual.
    add_multicomponent_vector_form(new WeakFormsElasticity::DefaultMultiComponentVolumetricResidualFormLinearSym(
         Hermes::vector<unsigned int>(0, 1), lambda, mu));

    add_vector_form(new WeakFormsElasticity::DefaultVolumetricResidualFormLinear_x_y(0, lambda, mu));
    add_vector_form(new WeakFormsElasticity::DefaultVolumetricResidualFormLinear_y_x(1, lambda, mu));

    // Gravity loading.
    add_vector_form(new WeakFormsH1::VolumetricVectorForms::DefaultVectorFormVol(1, HERMES_ANY, -rho_g));

    // External forces.
    add_multicomponent_vector_form_surf(new WeakFormsH1::SurfaceVectorForms::DefaultMultiComponentVectorFormSurf(
         Hermes::vector<unsigned int>(0, 1), non_zero_neumann_bnd, Hermes::vector<double>(-f0, -f1)));
    */
  };
};
