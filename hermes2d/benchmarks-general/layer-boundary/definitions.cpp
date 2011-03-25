#include "weakform/weakform.h"
#include "weakform_library/laplace.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

/* Custom function that is used in the exact solution and in right-hand side */

// Exact solution to the 1D problem -u'' + K*K*u = K*K in (-1,1) with zero Dirichlet BC.
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

/* Exact solution */

// Exact solution to the 2D equation.
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

// Right-hand side for the 2D equation -Laplace u + K*K*u = 0 with zero Dirichlet BC.

class CustomRightHandSide: public DefaultNonConstRightHandSide
{
public:
  CustomRightHandSide(double coeff1) : DefaultNonConstRightHandSide(), coeff1(coeff1) {
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

// Weak forms for the 2D equation with Dirichlet boundary conditions.
class CustomWeakFormPerturbedPoisson : public WeakForm
{
public:
  CustomWeakFormPerturbedPoisson(CustomRightHandSide* rhs) : WeakForm(1) {
    add_matrix_form(new DefaultMatrixFormStiffness(0, 0));
    add_matrix_form(new DefaultMatrixFormMass(0, 0, rhs->coeff1*rhs->coeff1));
    add_vector_form(new DefaultVectorFormVolNonConst(0, rhs));
  };
};
