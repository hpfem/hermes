#include "hermes2d.h"

/* Weak forms */

class CustomWeakFormPoissonNewton : public WeakForm
{
public:
  CustomWeakFormPoissonNewton(double lambda, double alpha, double T0, std::string bdy_heat_flux);
};