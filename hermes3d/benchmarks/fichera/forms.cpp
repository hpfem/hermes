template<typename real, typename scalar>
scalar bilinear_form(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, Func<real> *v, Geom<real> *e, ExtData<scalar> *data) 
{
  return int_grad_u_grad_v<real, scalar>(n, wt, u, v, e);
}

template<typename real>
real rhs(real x, real y, real z) 
{
  return -0.75 * pow(x*x + y*y + z*z, -0.75);
}

template<typename real, typename scalar>
scalar linear_form(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, Geom<real> *e, ExtData<scalar> *data) 
{
  return int_F_v<real, scalar>(n, wt, rhs, u, e);
}
