#include "weakform/weakform.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/h1.h"

using namespace WeakFormsH1;

/* Right-hand side */

class CustomRightHandSide: public DefaultFunction
{
public:
  CustomRightHandSide(double alpha, double x_loc, double y_loc, double r_zero)
    : DefaultFunction(), alpha(alpha), x_loc(x_loc), y_loc(y_loc), r_zero(r_zero) { };

  virtual double value(double x, double y) const {
    double a = pow(x - x_loc, 2);
    double b = pow(y - y_loc, 2);
    double c = sqrt(a + b);
    double d = ((alpha*x - (alpha * x_loc)) * (2*x - (2 * x_loc)));
    double e = ((alpha*y - (alpha * y_loc)) * (2*y - (2 * y_loc)));
    double f = (pow(alpha*c - (alpha * r_zero), 2) + 1.0);
    double g = (alpha * c - (alpha * r_zero));

    return -(((alpha/(c * f)) - (d/(2 * pow(a + b, 1.5) * f))
           - ((alpha * d * g)/((a + b) * pow(f, 2))) +
        (alpha/(c * f)) - (e/(2 * pow(a + b, 1.5) * f))
           - ((alpha * e * g)/((a + b) * pow(f, 2)))));
  }

  virtual Ord ord (Ord x, Ord y) const {
    return Ord(8);
  }
  double alpha, x_loc, y_loc, r_zero;
};


/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double alpha, double x_loc, double
                      y_loc, double r_zero)
             : ExactSolutionScalar(mesh), alpha(alpha), x_loc(x_loc),
                                   y_loc(y_loc), r_zero(r_zero) { }

  virtual scalar value(double x, double y) const {
    return atan(alpha * (sqrt(pow(x - x_loc, 2) + pow(y - y_loc, 2)) - r_zero));
  };

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    double a = pow(x - x_loc, 2);
    double b = pow(y - y_loc, 2);
    double c = sqrt(a + b);
    double d = (alpha*x - (alpha * x_loc));
    double e = (alpha*y - (alpha * y_loc));
    double f = (pow(alpha*c - (alpha * r_zero), 2) + 1.0);

    dx = (d/(c * f));
    dy = (e/(c * f));
  };

  virtual Ord ord (Ord x, Ord y) const {
    return Ord(8);
  }

  double alpha, x_loc, y_loc, r_zero;
};

/* Weak forms */

class CustomWeakFormPoisson : public WeakForm
{
public:
  CustomWeakFormPoisson(CustomRightHandSide* rhs) : WeakForm(1) {
    add_matrix_form(new DefaultJacobianDiffusion(0, 0));
    add_vector_form(new DefaultVectorFormVol(0, HERMES_ANY, 1.0, rhs));
  };
};

