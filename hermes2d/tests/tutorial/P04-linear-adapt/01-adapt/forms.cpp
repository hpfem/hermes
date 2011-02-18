template<typename Real, typename Scalar>
Scalar biform1(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
               Geom<Real> *e, ExtData<Scalar> *ext)
{
  return EPS_1 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar biform2(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
               Geom<Real> *e, ExtData<Scalar> *ext)
{
  return EPS_2 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}
