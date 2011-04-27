#include "weakform/weakform.h"
#include "weakform_library/h1.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"

using namespace WeakFormsH1;

/* Initial condition */

class InitialSolutionHeatTransfer : public ExactSolutionScalar
{
public:
  InitialSolutionHeatTransfer(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual scalar value (double x, double y) const {
    return (x+10)*(y+10)/100. + 2.;
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = (y+10)/10.;
    dy = (x+10)/10.;
  }

  virtual Ord ord(Ord x, Ord y) const {
    return (x+10)*(y+10)/100. + 2.;
  }
};

/* Essential BC */

class CustomEssentialBCNonConst : public EssentialBoundaryCondition {
public:
  CustomEssentialBCNonConst(std::string marker) : EssentialBoundaryCondition(Hermes::vector<std::string>(marker)) { }

  ~CustomEssentialBCNonConst() { };

  inline EssentialBCValueType get_value_type() const {
    return EssentialBoundaryCondition::BC_FUNCTION;
  }

  virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const {
    return (x+10)*(y+10)/100.;
  }
};
