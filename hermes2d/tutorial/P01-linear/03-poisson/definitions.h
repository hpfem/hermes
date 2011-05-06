#include "hermes2d.h"
#include "weakform/weakform.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/h1.h"

/* Weak forms */

class CustomWeakFormPoisson : public WeakForm
{
public:
  CustomWeakFormPoisson(std::string mat_al, double lambda_al,
                        std::string mat_cu, double lambda_cu,
                        double vol_heat_src);
};
