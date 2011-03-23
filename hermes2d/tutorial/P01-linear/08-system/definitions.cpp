#include "weakform/weakform.h"
#include "weakform_library/linear_elasticity.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

class CustomWeakFormLinearElasticity : public DefaultWeakFormLinearElasticity
{
public:
  CustomWeakFormLinearElasticity(double E, double nu, double rho_g, 
                                 std::string non_zero_neumann_bnd, double f0, double f1) 
            : DefaultWeakFormLinearElasticity(E, nu, rho_g) 
  {
    add_vector_form_surf(new VectorFormSurfForce_0(non_zero_neumann_bnd, f0));
    add_vector_form_surf(new VectorFormSurfForce_1(non_zero_neumann_bnd, f1));
  };
};
