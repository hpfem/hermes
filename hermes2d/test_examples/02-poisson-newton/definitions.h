#include "hermes2d.h"

/* Weak forms */

class CustomWeakFormPoissonNewton : public Hermes::Hermes2D::WeakForm<double>
{
public:
  CustomWeakFormPoissonNewton(std::string mat_al, Hermes::Hermes1DFunction<double>* lambda_al,
                              std::string mat_cu, Hermes::Hermes1DFunction<double>* lambda_cu,
                              Hermes::Hermes2DFunction<double>* vol_src_term, std::string bdy_heat_flux,
                              double alpha, double t_exterior);
};

/* Custom non-constant Dirichlet condition */

class CustomDirichletCondition : public Hermes::Hermes2D::EssentialBoundaryCondition<double>
{
public:
  CustomDirichletCondition(Hermes::vector<std::string> markers, double A, double B, double C);

  virtual EssentialBCValueType get_value_type() const;

  virtual double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;

  protected:
    double A, B, C;
};

