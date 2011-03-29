#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/h1.h"

using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::VolumetricVectorForms;
using namespace WeakFormsH1::RightHandSides;

/*  Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh) : ExactSolutionScalar(mesh) { };

  virtual scalar value (double x, double y) const {
    double r = sqrt(x*x + y*y);
    double a = atan2(x, y);
    return pow(r, 2.0/3.0) * sin(2.0*a/3.0 + M_PI/3);
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    double t1 = 2.0/3.0*atan2(x, y) + M_PI/3;
    double t2 = pow(x*x + y*y, 1.0/3.0);
    double t3 = x*x * ((y*y)/(x*x) + 1);
    dx = 2.0/3.0*x*sin(t1)/(t2*t2) + 2.0/3.0*y*t2*cos(t1)/t3;
    dy = 2.0/3.0*y*sin(t1)/(t2*t2) - 2.0/3.0*x*t2*cos(t1)/t3;
  };

  virtual Ord ord(Ord x, Ord y) const {
    Ord r = x;
    Ord a = Ord(20);
    return pow(r, 2.0/3.0) * sin(2.0*a/3.0 + M_PI/3);
  }
};

/* Weak forms */

class CustomWeakFormPoisson : public WeakForm
{
public:
    CustomWeakFormPoisson() : WeakForm(1) {
        add_matrix_form(new DefaultLinearDiffusion(0, 0));
        add_vector_form(new DefaultVectorFormConst(0, 0));
    }
};
