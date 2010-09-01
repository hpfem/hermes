// return the value \int \nabla u . \nabla v dx
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
                     Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

// return the value \int v dx
template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, 
                   Geom<Real> *e, ExtData<Scalar> *ext)
{
  return CONST_F*int_v<Real, Scalar>(n, wt, v);
}
