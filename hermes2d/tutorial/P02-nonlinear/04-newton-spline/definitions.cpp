#include "hermes2d.h"

using namespace WeakFormsH1;

/* Initial condition for the Newton's method */

class CustomInitialCondition : public ExactSolutionScalar
{
public:
  CustomInitialCondition(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = (y+10)/100.;
    dy = (x+10)/100.;
  };

  virtual scalar value (double x, double y) const {
    return (x+10)*(y+10)/100. + 2;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return x*y;
  }
};

/* Essential boundary conditions */

class CustomEssentialBCNonConst : public EssentialBoundaryCondition {
public:
  CustomEssentialBCNonConst(std::string marker)
           : EssentialBoundaryCondition(Hermes::vector<std::string>(marker)) { }

  ~CustomEssentialBCNonConst() {};

  inline EssentialBCValueType get_value_type() const {
    return EssentialBoundaryCondition::BC_FUNCTION;
  }

  virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const
  {
    return (x+10)*(y+10)/100.;
  }
};
