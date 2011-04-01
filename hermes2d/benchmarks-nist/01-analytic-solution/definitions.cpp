#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/h1.h"

using namespace WeakFormsH1;
using namespace WeakFormsH1::RightHandSides;
using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::VolumetricVectorForms;

/* Right-hand side */

class CustomRightHandSide: public DefaultNonConstRightHandSide
{
public:
  CustomRightHandSide(double pol_deg) 
    : DefaultNonConstRightHandSide(), pol_deg(pol_deg) {};

  virtual double value(double x, double y) const {
    double a = pow(2.0, 4.0*pol_deg);
    double b = pow(x-1.0, 8.0);
    double c = (38.0*pow(x, 2.0) - 38.0*x + 9.0);
    double d = pow(y-1.0, pol_deg);
    double e = pow(y-1.0, 8.0);
    double f = (38.0*pow(y, 2.0) - 38.0*y + 9.0);
    double g = pow(x-1.0, pol_deg);

    return -(pol_deg*a*pow(x, 8.0)*b*c*pow(y, pol_deg)*d 
	     + pol_deg*a*pow(y, 8.0)*e*f*pow(x, pol_deg)*g);
  }

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(std::max(8, (int)(pol_deg+0.51)));
  }

  double pol_deg;
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double poly_deg) 
            : ExactSolutionScalar(mesh), poly_deg(poly_deg) {};

  double value(double x, double y) const {
    return pow(2, 4 * poly_deg) * pow(x, poly_deg) * pow(1 - x, poly_deg) 
           * pow(y, poly_deg) * pow(1 - y, poly_deg);
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    double A = pow((1.0-y), poly_deg);
    double B = pow((1.0-x), poly_deg);
    double C = pow(y, poly_deg);
    double D = pow(x, poly_deg);

    dx = ((poly_deg * pow(16.0, poly_deg)*A*C) / (x-1.0) 
         + (poly_deg*pow(16.0, poly_deg)*A*C)/x)*B*D;
    dy = ((poly_deg*pow(16.0, poly_deg)*B*D)/(y-1.0)
         + (poly_deg*pow(16.0, poly_deg)*B*D)/y)*A*C;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(poly_deg);
  }

  double poly_deg;
};

/* Weak forms */

class CustomWeakFormPoisson : public WeakForm
{
public:
  CustomWeakFormPoisson(DefaultNonConstRightHandSide* rhs) : WeakForm(1) {
    add_matrix_form(new DefaultLinearDiffusion(0, 0));
    add_vector_form(new DefaultVectorFormNonConst(0, rhs));
  };
};

