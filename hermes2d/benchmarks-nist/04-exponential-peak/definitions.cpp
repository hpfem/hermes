#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/laplace.h"

/* Right-hand side */

class CustomRightHandSide: public DefaultNonConstRightHandSide
{
public:
  CustomRightHandSide(double coeff1, double coeff2, double coeff3) 
    : DefaultNonConstRightHandSide(coeff1, coeff2, coeff3) {};

  virtual double value(double x, double y) {
    double a_P = (-coeff1 * pow((x - coeff2), 2) - coeff1 * pow((y - coeff2), 2));
  
    return -(4 * exp(a_P) * coeff1 * (coeff1 * (x - coeff2) * (x - coeff2) 
				      + ALPHA_P * (y - coeff3) * (y - coeff3) - 1));
  }

  virtual Ord ord(Ord x, Ord y) {
    return Ord(8);
  }
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double ALPHA_P, double X_LOC, double Y_LOC) 
         : ExactSolutionScalar(mesh), ALPHA_P(ALPHA_P), X_LOC(X_LOC), Y_LOC(Y_LOC) {};

  virtual scalar value (double x, double y) {
    return exp(-ALPHA_P * (pow((x - X_LOC), 2) + pow((y - Y_LOC), 2)));
  };

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) {
    double a = -ALPHA_P * ( (x - X_LOC) * (x - X_LOC) + (y - Y_LOC) * (y - Y_LOC));
    dx = -exp(a) * (2 * ALPHA_P * (x - X_LOC));
    dy = -exp(a) * (2 * ALPHA_P * (y - Y_LOC));
  };

  virtual Ord ord (Ord x, Ord y) {
    return Ord(8);
  };

  double ALPHA_P;
  double X_LOC;
  double Y_LOC;
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
