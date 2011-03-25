#include "weakform/weakform.h"
#include "weakform_library/laplace.h"
#include "weakform_library/linear_elasticity.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

class CustomWeakFormLinearElasticity : public WeakForm
{
public:
  CustomWeakFormLinearElasticity(double E, double nu, double rho_g, 
                                 std::string non_zero_neumann_bnd, double f0, double f1) 
            : WeakForm(2)
  {
    double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
    double mu = E / (2*(1 + nu));

    add_matrix_form(new MatrixFormVolLinearElasticity_x_x(0, 0, lambda, mu));
    add_matrix_form(new MatrixFormVolLinearElasticity_x_y(0, 1, lambda, mu));
    add_matrix_form(new MatrixFormVolLinearElasticity_y_y(1, 1, lambda, mu));

    // gravity loading
    add_vector_form(new DefaultVectorFormVolConst(1, rho_g));

    // external forces
    add_vector_form_surf(new DefaultVectorFormSurf(0, non_zero_neumann_bnd, f0));
    add_vector_form_surf(new DefaultVectorFormSurf(1, non_zero_neumann_bnd, f1));
  };
};
