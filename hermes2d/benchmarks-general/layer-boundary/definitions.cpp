#include "weakform/weakform.h"
#include "weakform_library/h1.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::VolumetricVectorForms;
using namespace WeakFormsH1::RightHandSides;

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

class CustomWeakFormPerturbedPoisson : public WeakForm
{
public:
  CustomWeakFormPerturbedPoisson(CustomRightHandSide* rhs) : WeakForm(1) {
    add_matrix_form(new DefaultLinearDiffusion(0, 0));
    add_matrix_form(new DefaultLinearMass(0, 0, rhs->coeff1*rhs->coeff1));
    add_vector_form(new DefaultVectorFormNonConst(0, rhs));
  };
};
