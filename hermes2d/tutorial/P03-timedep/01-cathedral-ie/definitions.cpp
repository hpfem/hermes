#include "weakform/weakform.h"
#include "weakform_library/h1.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

using namespace WeakFormsH1;
using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::SurfaceMatrixForms;

class CustomWeakFormHeatRK1 : public WeakForm
{
public:
  CustomWeakFormHeatRK1(std::string bdy_air, double alpha, double lambda, double heatcap, double rho, 
                        double time_step, double* current_time_ptr, double temp_init, double t_final, 
                        Solution* prev_time_sln) : WeakForm(1)
  {
    add_matrix_form(new DefaultLinearDiffusion(0, 0, lambda));
    add_matrix_form(new DefaultLinearMass(0, 0, heatcap * rho / time_step));
    CustomVectorFormVolHeatRK1* vec_form_vol = new CustomVectorFormVolHeatRK1(0, heatcap, rho, time_step);
    vec_form_vol->ext.push_back(prev_time_sln);
    add_vector_form(vec_form_vol);

    add_matrix_form_surf(new DefaultMatrixFormSurf(0, 0, bdy_air, alpha * lambda));
    add_vector_form_surf(new CustomVectorFormSurfHeatRK1(0, bdy_air, alpha, lambda, current_time_ptr, temp_init, t_final));
  };

private:
  // This form is custom since it contains previous time-level solution.
  class CustomVectorFormVolHeatRK1 : public WeakForm::VectorFormVol
  {
  public:
    CustomVectorFormVolHeatRK1(int i, double heatcap, double rho, double time_step) 
      : WeakForm::VectorFormVol(i), heatcap(heatcap), rho(rho), time_step(time_step) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Func<Real>* temp_prev_time = ext->fn[0];
      return heatcap * rho * int_u_v<Real, Scalar>(n, wt, temp_prev_time, v) / time_step;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    double alpha, heatcap, rho, time_step;
  };

  // This form is custom since it contains time-dependent exterior temperature.
  class CustomVectorFormSurfHeatRK1 : public WeakForm::VectorFormSurf
  {
  private:
      double h;
  public:
    CustomVectorFormSurfHeatRK1(int i, std::string area, double alpha, double lambda, 
                                double* current_time_ptr, double temp_init, double t_final) 
      : WeakForm::VectorFormSurf(i, area), alpha(alpha), lambda(lambda), current_time_ptr(current_time_ptr), 
                                 temp_init(temp_init), t_final(t_final) { }

    template<typename Real, typename Scalar>
    Scalar vector_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      return lambda * alpha * temp_ext(*current_time_ptr + time_step) * int_v<Real, Scalar>(n, wt, v);
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
        return vector_form_surf<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
        return vector_form_surf<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    // Time-dependent exterior temperature.
    template<typename Real>
    Real temp_ext(Real t) const {
      return temp_init + 10. * sin(2*M_PI*t/t_final);
    }

    double alpha, lambda, *current_time_ptr, temp_init, t_final;
  };
};

