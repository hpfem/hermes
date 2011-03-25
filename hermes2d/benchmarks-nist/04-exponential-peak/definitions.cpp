#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/laplace.h"

/* Right-hand side */

class CustomRightHandSide: public DefaultNonConstRightHandSide
{
public:
  CustomRightHandSide(double alpha, double x_loc, double y_loc) 
    : DefaultNonConstRightHandSide(), alpha(alpha), x_loc(x_loc), y_loc(y_loc) {};

  virtual double value(double x, double y) const {
    double a_P = (-alpha * pow((x - x_loc), 2) - alpha * pow((y - y_loc), 2));
  
    return -(4 * exp(a_P) * alpha * (alpha * (x - x_loc) * (x - x_loc) 
				      + alpha * (y - y_loc) * (y - y_loc) - 1));
  }

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(8);
  }

  double alpha, x_loc, y_loc;
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double ALPHA_P, double X_LOC, double Y_LOC) 
         : ExactSolutionScalar(mesh), ALPHA_P(ALPHA_P), X_LOC(X_LOC), Y_LOC(Y_LOC) {};

  virtual scalar value (double x, double y) const {
    return exp(-ALPHA_P * (pow((x - X_LOC), 2) + pow((y - Y_LOC), 2)));
  };

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    double a = -ALPHA_P * ( (x - X_LOC) * (x - X_LOC) + (y - Y_LOC) * (y - Y_LOC));
    dx = -exp(a) * (2 * ALPHA_P * (x - X_LOC));
    dy = -exp(a) * (2 * ALPHA_P * (y - Y_LOC));
  };

  virtual Ord ord (Ord x, Ord y) const {
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
    add_matrix_form(new DefaultMatrixFormStiffness(0, 0));
    add_vector_form(new DefaultVectorFormVolNonConst(0, rhs));
  };
};
