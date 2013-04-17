#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::WeakFormsH1;
using namespace Hermes::Hermes2D::Views;

class CustomWeakFormPoisson : public WeakForm<double>
{
public:
  CustomWeakFormPoisson(const std::string& mat_motor, double eps_motor, 
                        const std::string& mat_air, double eps_air);
};
