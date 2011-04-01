#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/h1.h"

using namespace WeakFormsH1;
using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::VolumetricVectorForms;
using namespace WeakFormsH1::RightHandSides;

/* Right-hand side */

class CustomExactFunction {
public:
  CustomExactFunction(double k, double alpha): k(k), alpha(alpha) { };

  double fn(double x, double y) {
    if (x <= 0) return cos(k * y);
    else return cos(k * y) + pow(x, alpha);
  }

  double k, alpha;
};

class CustomRightHandSide : public DefaultNonConstRightHandSide
{
public:
  CustomRightHandSide(double k, double alpha) 
    : DefaultNonConstRightHandSide(), k(k), alpha(alpha) {
    cef = new CustomExactFunction(k, alpha);
  };

  virtual double value(double x, double y) const {
    if (x < 0) return cef->fn(x, y) * k * k;
    else return cef->fn(x, y) * k * k 
                - alpha *(alpha - 1) * pow(x, alpha - 2.) 
                - k * k * pow(x, alpha);
  }

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(20);
  }

  ~CustomRightHandSide() { delete cef; }

  CustomExactFunction* cef; 
  double k, alpha;
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double k, double alpha) 
    : ExactSolutionScalar(mesh), k(k), alpha(alpha) {
    cef = new CustomExactFunction(k, alpha);
  };

  virtual double value(double x, double y) const {
    return cef->fn(x, y);
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    if (x <= 0) dx = 0;
    else dx = alpha * pow(x, alpha - 1);
    dy = -sin(k * y) * k;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(20);
  }

  ~CustomExactSolution() { delete cef; }

  CustomExactFunction* cef; 
  double k, alpha;
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

