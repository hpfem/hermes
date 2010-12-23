template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real>
Real rhs(Real x, Real y)
{
  Real a = pow(2.0, 4.0*EXACT_SOL_P);
  Real b = pow(x-1.0, 8.0);
  Real c = (38.0*pow(x, 2.0) - 38.0*x + 9.0);
  Real d = pow(y-1.0, EXACT_SOL_P);
  Real e = pow(y-1.0, 8.0);
  Real f = (38.0*pow(y, 2.0) - 38.0*y + 9.0);
  Real g = pow(x-1.0, EXACT_SOL_P);

  return EXACT_SOL_P*a*pow(x, 8.0)*b*c*pow(y, EXACT_SOL_P)*d + EXACT_SOL_P*a*pow(y, 8.0)*e*f*pow(x,EXACT_SOL_P)*g;
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return -int_F_v<Real, Scalar>(n, wt, rhs, v, e);
}
