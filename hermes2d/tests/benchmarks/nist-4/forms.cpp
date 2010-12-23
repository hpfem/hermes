template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real>
Real rhs(Real x, Real y)
{
  Real a_P = (-ALPHA_P * pow((x - X_LOC), 2) - ALPHA_P * pow((y - Y_LOC), 2));
  
  return 4 * exp(a_P) * ALPHA_P * (ALPHA_P * (x - X_LOC) * (x - X_LOC) + ALPHA_P * (y - Y_LOC) * (y - Y_LOC) - 1);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return -int_F_v<Real, Scalar>(n, wt, rhs, v, e);
}
