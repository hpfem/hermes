#include "hermes2d.h"

using namespace WeakFormsH1;

class CustomWeakFormHeatRK1 : public WeakForm
{
public:
  CustomWeakFormHeatRK1(std::string bdy_air, double alpha, double lambda, double heatcap, double rho,
                        double time_step, double* current_time_ptr, double temp_init, double t_final,
                        Solution* prev_time_sln) : WeakForm(1)
  {
    /* Jacobian */
    // Contribution of the time derivative term.
    add_matrix_form(new DefaultMatrixFormVol(0, 0, HERMES_ANY, 1.0 / time_step));
    // Contribution of the diffusion term.
    add_matrix_form(new DefaultJacobianDiffusion(0, 0, HERMES_ANY, lambda / (rho * heatcap)));
    // Contribution of the Newton boundary condition.
    add_matrix_form_surf(new DefaultMatrixFormSurf(0, 0, bdy_air, alpha / (rho * heatcap)));

    // Residual.
    // Contribution of the time derivative term.
    add_vector_form(new DefaultResidualVol(0, HERMES_ANY, 1.0 / time_step));
    // Contribution of the diffusion term.
    add_vector_form(new DefaultResidualDiffusion(0, HERMES_ANY, lambda / (rho * heatcap)));
    CustomVectorFormVol* vec_form_vol = new CustomVectorFormVol(0, time_step);
    vec_form_vol->ext.push_back(prev_time_sln);
    add_vector_form(vec_form_vol);
    // Contribution of the Newton boundary condition.
    add_vector_form_surf(new DefaultResidualSurf(0, bdy_air, alpha / (rho * heatcap)));
    // Contribution of the Newton boundary condition.
    add_vector_form_surf(new CustomVectorFormSurf(0, bdy_air, alpha, rho, heatcap,
                         current_time_ptr, temp_init, t_final));
  };

private:
  // This form is custom since it contains previous time-level solution.
  class CustomVectorFormVol : public WeakForm::VectorFormVol
  {
  public:
    CustomVectorFormVol(int i, double time_step)
      : WeakForm::VectorFormVol(i), time_step(time_step) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Func<Real>* temp_prev_time = ext->fn[0];
      return -int_u_v<Real, Scalar>(n, wt, temp_prev_time, v) / time_step;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    double time_step;
  };

  // This form is custom since it contains time-dependent exterior temperature.
  class CustomVectorFormSurf : public WeakForm::VectorFormSurf
  {
  public:
    CustomVectorFormSurf(int i, std::string area, double alpha, double rho, double heatcap,
                                double* current_time_ptr, double temp_init, double t_final)
      : WeakForm::VectorFormSurf(i, area), alpha(alpha), rho(rho), heatcap(heatcap), current_time_ptr(current_time_ptr),
                                 temp_init(temp_init), t_final(t_final) { }

    template<typename Real, typename Scalar>
    Scalar vector_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      return -alpha / (rho * heatcap) * temp_ext(*current_time_ptr + time_step) * int_v<Real>(n, wt, v);
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
        return vector_form_surf<double, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
        return vector_form_surf<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    // Time-dependent exterior temperature.
    template<typename Real>
    Real temp_ext(Real t) const {
      return temp_init + 10. * sin(2*M_PI*t/t_final);
    }

    double alpha, rho, heatcap, *current_time_ptr, temp_init, t_final;
  };
};

