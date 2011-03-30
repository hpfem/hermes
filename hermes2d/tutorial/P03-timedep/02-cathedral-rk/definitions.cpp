#include "weakform/weakform.h"
#include "weakform_library/h1.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

using namespace WeakFormsH1;
using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::SurfaceMatrixForms;

/* Weak forms */

class CustomWeakFormHeatRK : public WeakForm
{
public:
  CustomWeakFormHeatRK(std::string bdy_air, double alpha, double lambda, double heatcap, double rho, 
                       double* current_time_ptr, double temp_init, double t_final, 
                       Solution* prev_time_sln) : WeakForm(1)
  {
    // Jacobian volumetric part.
    add_matrix_form(new CustomFormJacobianVol(0, 0, lambda / (heatcap * rho)));

    // Jacobian surface part.
    add_matrix_form_surf(new CustomFormJacobianSurf(0, 0, bdy_air, alpha * lambda / (heatcap * rho)));

    // Residual - volumetric.
    CustomFormResidualVol* vec_form_vol = new CustomFormResidualVol(0, heatcap, rho, lambda);
    vec_form_vol->ext.push_back(prev_time_sln);
    add_vector_form(vec_form_vol);

    // Residual - surface.
    CustomFormResidualSurf* vec_form_surf = new CustomFormResidualSurf(0, bdy_air, alpha, lambda, rho, heatcap, 
                         current_time_ptr, temp_init, t_final);
    vec_form_surf->ext.push_back(prev_time_sln);
    add_vector_form_surf(vec_form_surf);
  };

private:
  // Overriding method clone() in class DefaultLinearDiffusion.
  class CustomFormJacobianVol : public DefaultLinearDiffusion
  {
  public:
    CustomFormJacobianVol(int i, int j, double coeff) 
      : DefaultLinearDiffusion(i, j, coeff) { }

      virtual WeakForm::MatrixFormVol* clone() {
        return new CustomFormJacobianVol(*this);
    }
  };

  // Overriding method clone() in class DefaultMatrixFormSurf.
  class CustomFormJacobianSurf : public DefaultMatrixFormSurf
  {
  public:
    CustomFormJacobianSurf(int i, int j, std::string area, double coeff) 
      : DefaultMatrixFormSurf(i, j, area, coeff) { }

      virtual WeakForm::MatrixFormSurf* clone() {
        return new CustomFormJacobianSurf(*this);
    }
  };

  // This form is custom since it contains previous time-level solution.
  class CustomFormResidualVol : public WeakForm::VectorFormVol
  {
  public:
    CustomFormResidualVol(int i, double heatcap, double rho, double lambda) 
      : WeakForm::VectorFormVol(i), heatcap(heatcap), rho(rho), lambda(lambda) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const {
      Func<Real>* K_sln = u_ext[0];
      Func<Real>* sln_prev_time = ext->fn[0];
      Scalar result = 0;
      for (int i = 0; i < n; i++) {
        Scalar sln_dx_i = sln_prev_time->dx[i] + K_sln->dx[i];
        Scalar sln_dy_i = sln_prev_time->dy[i] + K_sln->dy[i];
        result += -wt[i] * (sln_dx_i * v->dx[i] +  sln_dy_i * v->dy[i]);	       
      }

      return lambda / (heatcap * rho) * result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    virtual WeakForm::VectorFormVol* clone() {
      return new CustomFormResidualVol(*this);
    }

    double alpha, heatcap, rho, lambda;
  };

  // This form is custom since it contains time-dependent exterior temperature.
  class CustomFormResidualSurf : public WeakForm::VectorFormSurf
  {
  private:
      double h;
  public:
    CustomFormResidualSurf(int i, std::string area, double alpha, double lambda, double rho, 
                           double heatcap, double* current_time_ptr, double temp_init, double t_final) 
      : WeakForm::VectorFormSurf(i, area), alpha(alpha), lambda(lambda), rho(rho), 
	                         heatcap(heatcap), current_time_ptr(current_time_ptr), 
                                 temp_init(temp_init), t_final(t_final) { }

    template<typename Real, typename Scalar>
    Scalar vector_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                            Geom<Real> *e, ExtData<Scalar> *ext) const {
      Func<Scalar>* K_sln = u_ext[0];
      Func<Scalar>* sln_prev_time = ext->fn[0];

      Scalar result1 = get_current_stage_time() * int_v<Real, Scalar>(n, wt, v);
      Scalar result2 = 0;
      for (int i = 0; i < n; i++) {
        Scalar sln_val_i = sln_prev_time->val[i] + K_sln->val[i];
        result2 += wt[i] * sln_val_i * v->val[i];		       
      }
      
      return lambda * alpha / (rho * heatcap) * (result1 - result2);
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, 
                 ExtData<scalar> *ext) const {
        return vector_form_surf<scalar, scalar>(n, wt, u_ext, v, e, ext);
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
    double alpha, lambda, rho, heatcap, *current_time_ptr, temp_init, t_final;
  };
};

