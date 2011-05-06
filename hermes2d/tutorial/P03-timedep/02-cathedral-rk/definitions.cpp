#include "hermes2d.h"

using namespace WeakFormsH1;

/* Weak forms */

class CustomWeakFormHeatRK : public WeakForm
{
public:
  CustomWeakFormHeatRK(std::string bdy_air, double alpha, double lambda, double heatcap, double rho,
                       double* current_time_ptr, double temp_init, double t_final) : WeakForm(1)
  {
    // Jacobian volumetric part.
    add_matrix_form(new DefaultJacobianDiffusion(0, 0, HERMES_ANY, -lambda / (heatcap * rho)));

    // Jacobian surface part.
    add_matrix_form_surf(new DefaultMatrixFormSurf(0, 0, bdy_air, -alpha / (heatcap * rho)));

    // Residual - volumetric.
    add_vector_form(new DefaultResidualDiffusion(0, HERMES_ANY, -lambda / (heatcap * rho)));

    // Residual - surface.
    add_vector_form_surf(new CustomFormResidualSurf(0, bdy_air, alpha, rho, heatcap,
                current_time_ptr, temp_init, t_final));
  };

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
                                 temp_init(temp_init), t_final(t_final) { }

    template<typename Real, typename Scalar>
    Scalar vector_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                            Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar T_ext = temp_ext(get_current_stage_time());
      Scalar result = 0;

      for (int i = 0; i < n; i++) {
        result += wt[i] * (T_ext - u_ext[0]->val[i]) * v->val[i];
      }

      return alpha / (rho * heatcap) * result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e,
                         ExtData<scalar> *ext) const {
        return vector_form_surf<double, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
        return vector_form_surf<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    virtual WeakForm::VectorFormSurf* clone() {
      return new CustomFormResidualSurf(*this);
    }

    // Time-dependent exterior temperature.
    template<typename Real>
    Real temp_ext(Real t) const {
      return temp_init + 10. * sin(2*M_PI*t/t_final);
    }

    // Members.
    double alpha, rho, heatcap, *current_time_ptr, temp_init, t_final;
  };
};

