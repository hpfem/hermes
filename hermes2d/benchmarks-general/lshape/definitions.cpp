#include "weakform_library/laplace.h"

/* Exact solution */

// Exact solution to the L-Shape problem.
class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh) : ExactSolutionScalar(mesh) {};

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

/* Essential boundary conditions */

// Dirichlet boundary conditions (uses the exact solution).
#include "boundaryconditions/essential_bcs.h"
class EssentialBCNonConst : public EssentialBC {
public:
  EssentialBCNonConst(std::string marker, CustomExactSolution* exact_solution)
        : EssentialBC(Hermes::vector<std::string>()), exact_solution(exact_solution) {
    markers.push_back(marker);
  }

  ~EssentialBCNonConst() { }

  inline EssentialBCValueType get_value_type() const { return EssentialBC::BC_FUNCTION; }

  scalar function(double x, double y) const {
    return exact_solution->value(x, y);
  }

  // Member.
  CustomExactSolution* exact_solution;
};

