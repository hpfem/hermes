
template<typename Real>
Real rhs(Real x, Real y)
{
  Real u = atan(K*x);
  Real dudx = 1./(1 + (K*x)*(K*x)) * K;
  Real ddudxx = - K / (1 + (K*x)*(K*x)) / (1 + (K*x)*(K*x)) * 2. * K * K * x;
  return - ddudxx + u;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_vol(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) + 
         int_u_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf_right(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar dfdx_at_1 = 1./(1 + (K*1.)*(K*1.)) * K;
  return - dfdx_at_1 * int_v<Real, Scalar>(n, wt, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf_left(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar dfdx_at_minus_1 = -1./(1 + (-K*1.)*(-K*1.)) * K;
  return - dfdx_at_minus_1 * int_v<Real, Scalar>(n, wt, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_vol(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_F_v<Real, Scalar>(n, wt, rhs, v, e);
}


