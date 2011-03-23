#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

class CustomWeakFormHeatRK1 : public WeakForm
{
public:
  // Problem parameters.
  CustomWeakFormHeatRK1(std::string bdy_air, double alpha, double lambda, double heatcap, double rho, 
                    double time_step, double* current_time_ptr, double temp_init, double t_final, Solution* prev_time_sln) : WeakForm(1)
  {
    add_matrix_form(new MyMatrixFormVolHeatRK1(0, 0, alpha, lambda, heatcap, rho, time_step));
    MyVectorFormVolHeatRK1* vec_form_vol = new MyVectorFormVolHeatRK1(0, alpha, heatcap, rho, time_step);
    vec_form_vol->ext.push_back(prev_time_sln);
    add_vector_form(vec_form_vol);

    add_matrix_form_surf(new MyMatrixFormSurfHeatRK1(0, 0, bdy_air, alpha, lambda));
    add_vector_form_surf(new MyVectorFormSurfHeatRK1(0, bdy_air, alpha, lambda, current_time_ptr, temp_init, t_final));
  };

private:
  class MyMatrixFormVolHeatRK1 : public WeakForm::MatrixFormVol
  {
  public:
    MyMatrixFormVolHeatRK1(int i, int j, double alpha, double lambda, double heatcap, double rho, double time_step) 
      : WeakForm::MatrixFormVol(i, j), alpha(alpha), lambda(lambda), heatcap(heatcap), rho(rho), time_step(time_step) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return heatcap * rho * int_u_v<Real, Scalar>(n, wt, u, v) / time_step +
         lambda * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Members.
    double alpha, lambda, heatcap, rho, time_step;
  };

  class MyVectorFormVolHeatRK1 : public WeakForm::VectorFormVol
  {
  public:
    MyVectorFormVolHeatRK1(int i, double alpha, double heatcap, double rho, double time_step) 
      : WeakForm::VectorFormVol(i), alpha(alpha), heatcap(heatcap), rho(rho), time_step(time_step) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Func<Real>* temp_prev_time = ext->fn[0];
      return heatcap * rho * int_u_v<Real, Scalar>(n, wt, temp_prev_time, v) / time_step;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    // Members.
    double alpha, heatcap, rho, time_step;
  };

  class MyMatrixFormSurfHeatRK1 : public WeakForm::MatrixFormSurf
  {
  public:
    MyMatrixFormSurfHeatRK1(int i, int j, std::string bdy_air, double alpha, double lambda) 
      : WeakForm::MatrixFormSurf(i, j, bdy_air), alpha(alpha), lambda(lambda) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return lambda * alpha * int_u_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
        return matrix_form_surf<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
        return matrix_form_surf<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Members.
    double alpha, lambda;
  };

  class MyVectorFormSurfHeatRK1 : public WeakForm::VectorFormSurf
  {
  private:
      double h;
  public:
    MyVectorFormSurfHeatRK1(int i, std::string area, double alpha, double lambda, 
                            double* current_time_ptr, double temp_init, double t_final) 
      : WeakForm::VectorFormSurf(i, area), alpha(alpha), lambda(lambda), current_time_ptr(current_time_ptr), 
                                 temp_init(temp_init), t_final(t_final) { }

    template<typename Real, typename Scalar>
    Scalar vector_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return lambda * alpha * temp_ext(*current_time_ptr + time_step) * int_v<Real, Scalar>(n, wt, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
        return vector_form_surf<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
        return vector_form_surf<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    // Time-dependent exterior temperature.
    template<typename Real>
    Real temp_ext(Real t) {
      return temp_init + 10. * sin(2*M_PI*t/t_final);
    }

    // Members.
    double alpha, lambda, *current_time_ptr, temp_init, t_final;
  };
};

