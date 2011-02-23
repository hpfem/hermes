template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                     Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return HEATCAP * RHO * int_u_v<Real, Scalar>(n, wt, u, v) / time_step +
         LAMBDA * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                   Geom<Real> *e, ExtData<Scalar> *ext)
{
  Func<Real> *temp_prev = ext->fn[0];
  return HEATCAP * RHO * int_u_v<Real, Scalar>(n, wt, temp_prev, v) / time_step;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                          Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return LAMBDA * ALPHA * int_u_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                        Geom<Real> *e, ExtData<Scalar> *ext)
{
  return LAMBDA * ALPHA * temp_ext(current_time + time_step) * int_v<Real, Scalar>(n, wt, v);
}
