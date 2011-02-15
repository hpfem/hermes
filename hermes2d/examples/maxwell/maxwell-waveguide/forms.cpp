// bilinear form
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  cplx ikappa = cplx(0.0, kappa);
  return 1.0/mu_r * int_curl_e_curl_f<Real, Scalar>(n, wt, u, v) -
         ikappa * sqrt(mu_0 / e_0) * int_F_e_f<Real, Scalar>(n, wt, gam, u, v, e) -
         sqr(kappa) * int_F_e_f<Real, Scalar>(n, wt, er, u, v, e);
}

// surface linear form
template<typename Real, typename Scalar>
Scalar linear_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  cplx ii = cplx(0.0, 1.0);
  return ii * omega * J * int_v1<Real, Scalar>(n, wt, v); // just second component of v, since J = (0, J)
}

// error calculation
template<typename Real, typename Scalar>
Scalar hcurl_form_kappa(int n, double *wt, Func<Scalar> *u_ext[], Func<Scalar> *u, Func<Scalar> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_curl_e_curl_f<Scalar, Scalar>(n, wt, u, v) + sqr(kappa) * int_e_f<Scalar, Scalar>(n, wt, u, v);
}
