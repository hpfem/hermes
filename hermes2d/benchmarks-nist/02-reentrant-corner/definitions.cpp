// Exact solution.
class MyExactSolution : public ExactSolutionScalar
{
public:
  MyExactSolution(Mesh* mesh, int param) : ExactSolutionScalar(mesh), param(param) {
    if (param == 0) {
      OMEGA = ((5.0 * M_PI)/ 4.0);
      ALPHA = (M_PI/ OMEGA);
     }
    else if (param == 1) {
      OMEGA = ((3.0 * M_PI)/ 2.0);
      ALPHA = (M_PI/ OMEGA);
    }
    else if (param == 2) {
      OMEGA = ((7.0 * M_PI)/ 4.0);
      ALPHA = (M_PI/ OMEGA);
    }
    else{
      OMEGA = (2.0 * M_PI);
      ALPHA = (M_PI/ OMEGA);
    }
  };

  double get_angle(double y, double x) {
    double theta = atan2(y, x);
    if (theta < 0)
      theta += 2 * M_PI;
    return theta;
  };

  double value(double x, double y)
  {
    return (pow(sqrt(x*x + y*y), ALPHA) * sin(ALPHA * get_angle(y, x)));
  };

  // Function representing an exact one-dimension valued solution.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    double a = sqrt(x*x + y*y);
    double b = pow(a, (ALPHA - 1.0));
    double c = pow(a, ALPHA);
    double d = ((y*y)/(x*x) + 1.0 );

    dx = (((ALPHA* x* sin(ALPHA * get_angle(y,x)) *b)/a) 
         - ((ALPHA *y *cos(ALPHA * get_angle(y, x)) * c)/(pow(x, 2.0) *d)));
    dy = (((ALPHA* cos(ALPHA* get_angle(y, x)) *c)/(x * d)) 
         + ((ALPHA* y* sin(ALPHA* get_angle(y, x)) *b)/a));

    return value(x, y);
  };

  // Members.
  int param;
  double OMEGA;
  double ALPHA;
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

