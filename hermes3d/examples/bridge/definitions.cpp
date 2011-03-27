template<typename real, typename scalar>
scalar bilinear_form(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, 
                     Func<real> *v, Geom<real> *e, ExtData<scalar> *ext) 
{
  return int_grad_u_grad_v<real, scalar>(n, wt, u, v, e);
}

template<typename real, typename scalar>
scalar linear_form(int n, double *wt, Func<scalar> *u_ext[], Func<real> *v, 
                   Geom<real> *e, ExtData<scalar> *ext) 
{
  return - CONST_F * int_v<real, scalar>(n, wt, v, e);
}
