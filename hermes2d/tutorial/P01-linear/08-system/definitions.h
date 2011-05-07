#include "hermes2d.h"
#include "weakform/weakform.h"
#include "weakform_library/h1.h"
#include "weakform_library/elasticity.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"

//#define USE_MULTICOMPONENT_FORMS

class CustomWeakFormLinearElasticity : public WeakForm
{
public:
  CustomWeakFormLinearElasticity(double E, double nu, double rho_g,
                                 std::string surface_force_bdy, double f0, double f1);
};
