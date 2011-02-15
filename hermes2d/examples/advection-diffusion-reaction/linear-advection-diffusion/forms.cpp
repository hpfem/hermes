// Bilinear form.
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i=0; i < n; i++)
  {
    result += wt[i] * (EPSILON * (u->dx[i]*v->dx[i] + u->dy[i]*v->dy[i])
                               - (B1 * u->val[i] * v->dx[i] + B2 * u->val[i] * v->dy[i])
                      );
  }
  return result;
}

// Variational multiscale stabilization.
template<typename Real, typename Scalar>
Scalar bilinear_form_stabilization(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
  Real h_e = e->diam;
  Scalar result = 0;
  for (int i=0; i < n; i++) {
    double b_norm = sqrt(B1*B1 + B2*B2);
    Real tau = 1. / sqrt(9*pow(4*EPSILON/pow(h_e, 2), 2) + pow(2*b_norm/h_e, 2));
    result += wt[i] * tau * (-B1 * v->dx[i] - B2 * v->dy[i] + EPSILON * v->laplace[i])
                          * (-B1 * u->dx[i] - B2 * u->dy[i] + EPSILON * u->laplace[i]);
  }
  return result;
#else
  error("Define H2D_SECOND_DERIVATIVES_ENABLED in h2d_common.h if you want to use second derivatives of shape functions in weak forms.");
#endif
}

// Shock capturing.
template<typename Real, typename Scalar>
Scalar bilinear_form_shock_capturing(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v,
        Geom<Real> *e, ExtData<Scalar> *ext)
{
  Real h_e = e->diam;
  Real s_c = 0.9;
  Real result = 0;
  for (int i=0; i < n; i++) {
    // This R makes it nonlinear! So we need to use the Newton method:
    Real R_squared = pow(B1 * u->dx[i] + B2 * u->dy[i], 2.);
    Real R = sqrt(R_squared); //This just does fabs(B1 * u->dx[i] + B2 * u->dy[i]); but it can be parsed
    result += wt[i] * s_c * 0.5 * h_e * R *
              (u->dx[i]*v->dx[i] + u->dy[i]*v->dy[i]) /
              (sqrt(pow(u->dx[i], 2) + pow(u->dy[i], 2)) + 1.e-8);
  }
  return result;
}
