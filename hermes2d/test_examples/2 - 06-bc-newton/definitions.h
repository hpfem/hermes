#include "hermes2d.h"

/* Weak forms */

class CustomWeakFormPoissonNewton : public WeakForm<double>
{
public:
  CustomWeakFormPoissonNewton(std::string mat_al, HermesFunction<double>* lambda_al,
                              std::string mat_cu, HermesFunction<double>* lambda_cu,
                              HermesFunction<double>* vol_src_term, std::string bdy_heat_flux,
                              double alpha, double t_exterior);
};

/* Custom non-constant Dirichlet condition */

class CustomDirichletCondition : public EssentialBoundaryCondition<double> 
{
public:
  CustomDirichletCondition(Hermes::vector<std::string> markers, double A, double B, double C);

  virtual EssentialBCValueType get_value_type() const;

  virtual double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;

  protected:
    double A, B, C;
};

