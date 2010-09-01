template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real>
Real rhs(Real x, Real y)
{
  Real a = (-ALPHA * pow((x - X_LOC), 2) - ALPHA*pow((y - Y_LOC), 2));
  Real b = (2 * ALPHA * x - ALPHA);
  Real c = (2 * ALPHA * y - ALPHA);
  
  return exp(a) * pow(b,2) - 2 * ALPHA * exp(a) + exp(a) * pow(c,2) - 2 * ALPHA * exp(a);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return -int_F_v<Real, Scalar>(n, wt, rhs, v, e);
}
