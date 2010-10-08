template<typename real, typename scalar>
scalar bilinear_form(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, Func<real> *v, 
              Geom<real> *e, ExtData<scalar> *ext)
{
  return int_grad_u_grad_v<real, scalar>(n, wt, u, v, e) + K*K * int_u_v<real, scalar>(n, wt, u, v, e);
}

template<typename real, typename scalar>
scalar linear_form(int n, double *wt, Func<scalar> *u_ext[], Func<real> *v, 
              Geom<real> *e, ExtData<scalar> *ext)
{
  return int_F_v<real, scalar>(n, wt, rhs, v, e);
}

// integration order for the linear_form
Ord linear_form_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(24);
}
