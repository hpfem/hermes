#include "weakform/weakform.h"
#include "weakform_library/h1.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::VolumetricVectorForms;
using namespace WeakFormsH1::RightHandSides;

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double slope)
        : ExactSolutionScalar(mesh), slope(slope) { }

  virtual scalar value (double x, double y) const {
    return atan(slope * (sqrt(sqr(x-1.25) + sqr(y+0.25)) - M_PI/3));
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    double t = sqrt(sqr(x - 1.25) + sqr(y + 0.25));
    double u = t * (sqr(slope) * sqr(t - M_PI/3) + 1);
    dx = slope * (x-1.25) / u;
    dy = slope * (y+0.25) / u;
  }

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(20);
  }

  // Member.
  double slope;
};

/* Right-hand side */

class CustomRightHandSide: public DefaultNonConstRightHandSide
{
public:
  CustomRightHandSide(double coeff1) : DefaultNonConstRightHandSide(), coeff1(coeff1) { }

  virtual scalar value(double x, double y) const {
    double t2 = sqr(y + 0.25) + sqr(x - 1.25);
    double t = sqrt(t2);
    double u = (sqr(M_PI - 3.0*t)*sqr(coeff1) + 9.0);
    return 27.0/2.0 * sqr(2.0*y + 0.5) * (M_PI - 3.0*t) * pow(coeff1, 3.0) / (sqr(u) * t2) +
           27.0/2.0 * sqr(2.0*x - 2.5) * (M_PI - 3.0*t) * pow(coeff1, 3.0) / (sqr(u) * t2) -
           9.0/4.0 * sqr(2.0*y + 0.5) * coeff1 / (u * pow(t,3.0)) -
           9.0/4.0 * sqr(2.0*x - 2.5) * coeff1 / (u * pow(t,3.0)) +
           18.0 * coeff1 / (u * t);
  }

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(20);
  }

  double coeff1;
};

/* Weak forms */

class CustomWeakFormPoisson : public WeakForm
{
public:
  CustomWeakFormPoisson(double slope) : WeakForm(1) {
      add_matrix_form(new DefaultLinearDiffusion(0, 0));
      add_vector_form(new DefaultVectorFormNonConst(0, new CustomRightHandSide(-slope)));
  }
};

