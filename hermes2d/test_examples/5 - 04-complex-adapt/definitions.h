#include "hermes2d.h"
#include "boundaryconditions/essential_bcs.h"

/* Weak forms */

class CustomWeakForm : public WeakForm
{ 
public:
  CustomWeakForm(std::string mat_air,  double mu_air,
                 std::string mat_iron, double mu_iron, double gamma_iron,
                 std::string mat_wire, double mu_wire, scalar j_ext, double omega);
};
