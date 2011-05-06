#include "hermes2d.h"
#include "weakform/weakform.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/h1.h"

/* Weak forms */

class CustomWeakFormPoissonNeumann : public WeakForm
{
public:
  CustomWeakFormPoissonNeumann(std::string mat_al, double lambda_al,
                               std::string mat_cu, double lambda_cu,
                               double vol_heat_src, std::string bdy_heat_flux,
                               double heat_flux);
};

/* Custom non-constant Dirichlet condition */

class CustomDirichletCondition : public EssentialBoundaryCondition {
public:
  CustomDirichletCondition(Hermes::vector<std::string> markers, double A, double B, double C);

  virtual EssentialBoundaryCondition::EssentialBCValueType get_value_type() const;

  virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;

  protected:
    double A, B, C;
};

