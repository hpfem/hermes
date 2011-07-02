#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

/* Weak forms */

class CustomWeakForm : public WeakForm<std::complex<double> >
{ 
public:
  CustomWeakForm(std::string mat_air,  double mu_air,
                 std::string mat_iron, double mu_iron, double gamma_iron,
                 std::string mat_wire, double mu_wire, std::complex<double> j_ext, double omega);
};
