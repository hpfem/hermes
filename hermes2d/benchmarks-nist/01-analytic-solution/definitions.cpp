#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/laplace.h"

/* Right-hand side */

class CustomRightHandSide: public DefaultNonConstRightHandSide
{
public:
  CustomRightHandSide(double coeff1) : DefaultNonConstRightHandSide(coeff1) {};

  virtual double value(double x, double y) {
    double a = pow(2.0, 4.0*coeff1);
    double b = pow(x-1.0, 8.0);
    double c = (38.0*pow(x, 2.0) - 38.0*x + 9.0);
    double d = pow(y-1.0, coeff1);
    double e = pow(y-1.0, 8.0);
    double f = (38.0*pow(y, 2.0) - 38.0*y + 9.0);
    double g = pow(x-1.0, coeff1);

    return -(coeff1*a*pow(x, 8.0)*b*c*pow(y, coeff1)*d 
	     + coeff1*a*pow(y, 8.0)*e*f*pow(x, coeff1)*g);
  }

  virtual Ord ord(Ord x, Ord y) {
    return Ord(std::max(8, (int)(coeff1+0.5)));
  }
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double poly_deg) 
            : ExactSolutionScalar(mesh), poly_deg(poly_deg) {};

  double value(double x, double y) {
    return pow(2, 4 * poly_deg) * pow(x, poly_deg) * pow(1 - x, poly_deg) 
           * pow(y, poly_deg) * pow(1 - y, poly_deg);
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) {
    double A = pow((1.0-y), poly_deg);
    double B = pow((1.0-x), poly_deg);
    double C = pow(y, poly_deg);
    double D = pow(x, poly_deg);

    dx = ((poly_deg * pow(16.0, poly_deg)*A*C) / (x-1.0) 
         + (poly_deg*pow(16.0, poly_deg)*A*C)/x)*B*D;
    dy = ((poly_deg*pow(16.0, poly_deg)*B*D)/(y-1.0)
         + (poly_deg*pow(16.0, poly_deg)*B*D)/y)*A*C;
  };

  virtual Ord ord(Ord x, Ord y) {
    return Ord(poly_deg);
  }

  double poly_deg;
};

/* Weak forms */

class CustomWeakFormPoisson : public WeakForm
{
public:
  CustomWeakFormPoisson(DefaultNonConstRightHandSide* rhs) : WeakForm(1) {
    add_matrix_form(new DefaultMatrixFormVolConst(0, 0));
    add_vector_form(new DefaultVectorFormVolNonConst(0, rhs));
  };
};

