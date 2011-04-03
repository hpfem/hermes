#include "weakform/weakform.h"
#include "weakform_library/h1.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

using namespace WeakFormsH1;
using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::VolumetricVectorForms;

/* Weak forms */

class CustomWeakFormHeatTransferNewton : public WeakForm
{
public:
  CustomWeakFormHeatTransferNewton(CubicSpline* cspline, double heat_src) 
        : WeakForm(1) {
    // Jacobian.
    add_matrix_form(new DefaultJacobianNonlinearDiffusion(0, 0, cspline));
    // Residual.
    add_vector_form(new DefaultResidualNonlinearDiffusion(0, cspline));
    add_vector_form(new DefaultVectorFormConst(0, -heat_src));
  };
};

/* Initial condition for the Newton's method */

class CustomInitialSolutionHeatTransfer : public ExactSolutionScalar
{
public:
  CustomInitialSolutionHeatTransfer(Mesh* mesh) : ExactSolutionScalar(mesh) {};

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
           : EssentialBoundaryCondition(Hermes::vector<std::string>())
  {
    markers.push_back(marker);
  }

  ~CustomEssentialBCNonConst() {};

  inline EssentialBCValueType get_value_type() const { 
    return EssentialBoundaryCondition::BC_FUNCTION; 
  }

  virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const
  {
    return (x+10)*(y+10)/100.;
  }
};


