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
                                 std::string non_zero_neumann_bnd, double f0, double f1) 
            : WeakForm(2)
  {
    double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
    double mu = E / (2*(1 + nu));

    // There is one multi-component and one single-component form since we want to exploit symmetry of the forms.
    add_multicomponent_matrix_form(new WeakFormsElasticity::MultiComponentDefaultVolumetricMatrixFormLinearSym(Hermes::vector<std::pair<unsigned int, unsigned int> >(make_pair(0, 0), make_pair(1, 1)), lambda, mu));

    add_matrix_form(new WeakFormsElasticity::DefaultVolumetricMatrixFormLinear_x_y(0, 1, lambda, mu));
    
    add_multicomponent_vector_form(new WeakFormsElasticity::MultiComponentDefaultVolumetricResidualFormLinearSym(Hermes::vector<unsigned int>(0, 1), lambda, mu));

    add_vector_form(new WeakFormsElasticity::DefaultVolumetricResidualFormLinear_x_y(0, lambda, mu));

    // gravity loading
    add_vector_form(new WeakFormsH1::VolumetricVectorForms::DefaultVectorFormConst(1, -rho_g));

    // external forces
    add_multicomponent_vector_form_surf(new WeakFormsH1::SurfaceVectorForms::MultiComponentDefaultVectorFormSurf(Hermes::vector<unsigned int>(0, 1), non_zero_neumann_bnd, Hermes::vector<double>(-f0, -f1)));
  };
};
