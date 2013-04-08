#include "hermes2d.h"

/* Weak forms */

using namespace Hermes;
using namespace Hermes::Hermes2D;

class CustomWeakFormHeatRK : public WeakForm<double>
{
public:
  CustomWeakFormHeatRK(std::string bdy_air, double alpha, double lambda, double heatcap, double rho,
                       double* current_time_ptr, double temp_init, double t_final);

private:
  // This form is custom since it contains time-dependent exterior temperature.
  class CustomFormResidualSurf : public VectorFormSurf<double>
  {
  private:
      double h;
  public:
    CustomFormResidualSurf(int i, std::string area, double alpha, double rho,
                           double heatcap, double* current_time_ptr, double temp_init, double t_final)
          : VectorFormSurf<double>(i), alpha(alpha), rho(rho),
                                     heatcap(heatcap), current_time_ptr(current_time_ptr),
                                     temp_init(temp_init), t_final(t_final) 
    {
      this->set_area(area);
    };

    template<typename Real, typename Scalar>
    Scalar vector_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e,
                         Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const;

    virtual VectorFormSurf<double>* clone() const;

    // Time-dependent exterior temperature.
    template<typename Real>
    Real temp_ext(Real t) const;

    // Members.
    double alpha, rho, heatcap, *current_time_ptr, temp_init, t_final;
  };
};