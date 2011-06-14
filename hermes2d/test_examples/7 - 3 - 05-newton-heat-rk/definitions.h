#include "hermes2d.h"
#include "runge_kutta.h"

/* Weak forms */

class CustomWeakFormHeatRK : public WeakForm
{
public:
  CustomWeakFormHeatRK(std::string bdy_air, double alpha, double lambda, double heatcap, double rho,
                       double* current_time_ptr, double temp_init, double t_final);

private:
  // This form is custom since it contains time-dependent exterior temperature.
  class CustomFormResidualSurf : public WeakForm::VectorFormSurf
  {
  private:
      double h;
  public:
    CustomFormResidualSurf(int i, std::string area, double alpha, double rho,
                           double heatcap, double* current_time_ptr, double temp_init, double t_final)
          : WeakForm::VectorFormSurf(i, area), alpha(alpha), rho(rho),
                                     heatcap(heatcap), current_time_ptr(current_time_ptr),
                                     temp_init(temp_init), t_final(t_final) {};

    template<typename Real, typename Scalar>
    Scalar vector_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                            Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e,
                         ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::VectorFormSurf* clone();

    // Time-dependent exterior temperature.
    template<typename Real>
    Real temp_ext(Real t) const;

    // Members.
    double alpha, rho, heatcap, *current_time_ptr, temp_init, t_final;
  };
};

