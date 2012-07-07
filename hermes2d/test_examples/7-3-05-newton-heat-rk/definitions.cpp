#include "definitions.h"

CustomWeakFormHeatRK::CustomWeakFormHeatRK(std::string bdy_air, double alpha, double lambda, double heatcap, double rho,
                                           double* current_time_ptr, double temp_init, double t_final) : Hermes::Hermes2D::WeakForm<double>(1)
{
  // Jacobian volumetric part.
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(0, 0, HERMES_ANY, new Hermes1DFunction<double>(-lambda / (heatcap * rho))));

  // Jacobian surface part.
  add_matrix_form_surf(new WeakFormsH1::DefaultMatrixFormSurf<double>(0, 0, bdy_air, new Hermes2DFunction<double>(-alpha / (heatcap * rho))));

  // Residual - volumetric.
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(0, HERMES_ANY, new Hermes1DFunction<double>(-lambda / (heatcap * rho))));

  // Residual - surface.
  add_vector_form_surf(new CustomFormResidualSurf(0, bdy_air, alpha, rho, heatcap,
                       current_time_ptr, temp_init, t_final));
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormHeatRK::CustomFormResidualSurf::vector_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                                                                      Geom<Real> *e, ExtData<Scalar> *ext) const
{
  Scalar T_ext = Scalar(temp_ext(get_current_stage_time()));
  Scalar result = Scalar(0);

  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (T_ext - u_ext[0]->val[i]) * v->val[i];
  }

  return alpha / (rho * heatcap) * result;
}

double CustomWeakFormHeatRK::CustomFormResidualSurf::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e,
                                                           ExtData<double> *ext) const
{
  return vector_form_surf<double, double>(n, wt, u_ext, v, e, ext);
}

Ord CustomWeakFormHeatRK::CustomFormResidualSurf::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
{
  return vector_form_surf<Ord, Ord>(n, wt, u_ext, v, e, ext);
}

VectorFormSurf<double>* CustomWeakFormHeatRK::CustomFormResidualSurf::clone()
{
  return new CustomFormResidualSurf(*this);
}

template<typename Real>
Real CustomWeakFormHeatRK::CustomFormResidualSurf::temp_ext(Real t) const
{
  return temp_init + 10. * Hermes::sin(2*M_PI*t/t_final);
}