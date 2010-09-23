scalar rhs(scalar x, scalar y)
{
  if (x < 0) return fn(x, y)*K*K;
  else return fn(x, y)*K*K-ALPHA*(ALPHA-1)*pow(x, ALPHA - 2.) - K*K*pow(x, ALPHA);
}

template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt,Func<Scalar> *u_ext[],  Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

scalar linear_form(int n, double *wt, Func<scalar> *u_ext[], Func<scalar> *v, Geom<scalar> *e, ExtData<scalar> *ext)
{
  return int_F_v<scalar, scalar>(n, wt, rhs, v, e);
}

// integration order for the linear_form
Ord linear_form_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}
