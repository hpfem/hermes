// Exact solution to the L-Shape problem.
class MyExactSolution : public ExactSolutionScalar
{
public:
  MyExactSolution(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  static double value(double x, double y) {
    double r = sqrt(x*x + y*y);
    double a = atan2(x, y);
    return pow(r, 2.0/3.0) * sin(2.0*a/3.0 + M_PI/3);
  }

  // Function representing an exact scalar-valued solution.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    double t1 = 2.0/3.0*atan2(x, y) + M_PI/3;
    double t2 = pow(x*x + y*y, 1.0/3.0);
    double t3 = x*x * ((y*y)/(x*x) + 1);
    dx = 2.0/3.0*x*sin(t1)/(t2*t2) + 2.0/3.0*y*t2*cos(t1)/t3;
    dy = 2.0/3.0*y*sin(t1)/(t2*t2) - 2.0/3.0*x*t2*cos(t1)/t3;
    return value(x, y);
  };
};

// Dirichlet boundary conditions (uses the exact solution).
#include "boundaryconditions/essential_bcs.h"
class EssentialBCNonConstant : public EssentialBC {
public:
  EssentialBCNonConstant(std::string marker, MyExactSolution* exact_solution)
        : EssentialBC(Hermes::vector<std::string>()), exact_solution(exact_solution) {
    markers.push_back(marker);
  }

  ~EssentialBCNonConstant() { }

  inline EssentialBCValueType get_value_type() const { return EssentialBC::BC_FUNCTION; }

  scalar function(double x, double y) const {
    return exact_solution->value(x, y);
  }

  // Member.
  MyExactSolution* exact_solution;
};

