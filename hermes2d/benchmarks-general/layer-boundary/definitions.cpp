#include "hermes2d.h"

using namespace WeakFormsH1;

/* Exact solution */

class CustomExactFunction
{
public:
  CustomExactFunction(double K) : K(K) {};

  double uhat(double x) {
    return 1. - (exp(K*x) + exp(-K*x)) / (exp(K) + exp(-K));
  }
  double duhat_dx(double x) {
    return -K * (exp(K*x) - exp(-K*x)) / (exp(K) + exp(-K));
  }
  double dduhat_dxx(double x) {
    return -K*K * (exp(K*x) + exp(-K*x)) / (exp(K) + exp(-K));
  }

  // Member.
  double K;
};

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double K) : ExactSolutionScalar(mesh)
  {
    cef = new CustomExactFunction(K);
  };

  virtual scalar value (double x, double y) const {
    return cef->uhat(x) * cef->uhat(y);
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = cef->duhat_dx(x) * cef->uhat(y);
    dy = cef->uhat(x) * cef->duhat_dx(y);
  }

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(20);
  };

  ~CustomExactSolution() { delete cef;}

  CustomExactFunction* cef;
};

/* Right-hand side */

class CustomRightHandSide: public DefaultFunction
{
public:
  CustomRightHandSide(double coeff1) : DefaultFunction(), coeff1(coeff1) {
    cef = new CustomExactFunction(coeff1);
  };

  virtual scalar value(double x, double y) const {
    return -(cef->dduhat_dxx(x) * cef->uhat(y) + cef->uhat(x) * cef->dduhat_dxx(y))
           + coeff1 * coeff1 * cef->uhat(x) * cef->uhat(y);
  }

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(5);
  }

  ~CustomRightHandSide() { delete cef;}

  CustomExactFunction* cef;
  double coeff1;
};

/* Weak forms */

class CustomWeakFormPerturbedPoisson : public WeakForm
{
public:
  CustomWeakFormPerturbedPoisson(CustomRightHandSide* rhs) : WeakForm(1) {
    add_matrix_form(new DefaultJacobianDiffusion(0, 0));
    add_matrix_form(new DefaultMatrixFormVol(0, 0, HERMES_ANY, rhs->coeff1*rhs->coeff1));
    add_vector_form(new DefaultVectorFormVol(0, HERMES_ANY, 1.0, rhs));
  };
};
