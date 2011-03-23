#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/laplace.h"

/* Right-hand side */

// Right-hand side for the 2D equation -Laplace u = f with Dirichlet BC.
class CustomRightHandSide: public DefaultNonConstRightHandSide
{
public:
  CustomRightHandSide(double coeff1, double coeff2, double coeff3) 
    : DefaultNonConstRightHandSide(coeff1, coeff2, coeff3) {};

  virtual double value(double x, double y) {
    double a_P = (-coeff1 * pow((x - coeff2), 2) - coeff1 * pow((y - coeff2), 2));
  
    return 4 * exp(a_P) * coeff1 * (coeff1 * (x - coeff2) * (x - coeff2) 
           + ALPHA_P * (y - coeff3) * (y - coeff3) - 1);
  }

  virtual Ord ord(Ord x, Ord y) {
    return Ord(5);
  }
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double ALPHA_P, double X_LOC, double Y_LOC) 
         : ExactSolutionScalar(mesh), ALPHA_P(ALPHA_P), X_LOC(X_LOC), Y_LOC(Y_LOC) {};

  double fn(double x, double y) {
    return exp(-ALPHA_P * (pow((x - X_LOC), 2) + pow((y - Y_LOC), 2)));
  };

  // Function representing an exact scalar-valued solution.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    double a = -ALPHA_P * ( (x - X_LOC) * (x - X_LOC) + (y - Y_LOC) * (y - Y_LOC));

    dx = -exp(a) * (2 * ALPHA_P * (x - X_LOC));
    dy = -exp(a) * (2 * ALPHA_P * (y - Y_LOC));

    return fn(x, y);
  };

  // Members.
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

/* Essential boundary conditions (use values provided by exact solution) */

class EssentialBCNonConstExact : public EssentialBC
{
public:
  EssentialBCNonConstExact(std::string marker, CustomExactSolution* exact_solution) : 
        EssentialBC(Hermes::vector<std::string>()), exact_solution(exact_solution) 
  {
    markers.push_back(marker);
  };
  
  ~EssentialBCNonConstExact() {};

  virtual EssentialBCValueType get_value_type() const { 
    return BC_FUNCTION; 
  };

  virtual scalar function(double x, double y) const {
    return exact_solution->fn(x, y);
  };

  CustomExactSolution* exact_solution;
};
