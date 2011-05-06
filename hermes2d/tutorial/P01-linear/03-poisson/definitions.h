#include "hermes2d.h"

/* Weak forms */

class CustomWeakFormPoisson : public WeakForm
{
public:
  CustomWeakFormPoisson(std::string mat_al, double lambda_al,
                        std::string mat_cu, double lambda_cu,
                        double vol_heat_src);
};
