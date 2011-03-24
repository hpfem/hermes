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
    return - coeff1 * (coeff1 - 1.) * pow(x, coeff1 - 2.); 
  }

  virtual Ord ord(Ord x, Ord y) {
    return Ord((int)(coeff1 + 3));
  }
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double alpha) 
        : ExactSolutionScalar(mesh), alpha(alpha) {};

  virtual scalar value(double x, double y) {
    return pow(x, alpha);
  };

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) {
    dx = alpha * pow(x, alpha - 1.);
    dy = 0;
  };

  virtual Ord ord(Ord x, Ord y) {
    return Ord((int)(alpha + 1));
  }

  double alpha;
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

