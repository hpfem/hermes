template<typename real, typename scalar>
scalar bilinear_form(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *u, fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data) 
{
  return int_grad_u_grad_v<real, scalar>(n, wt, u, v, e) + int_u_v<real, scalar>(n, wt, u, v, e) / TAU;
}

template<typename real> real rhs(real x, real y, real z) 
{
  real ddxx = -2 * sin(TIME) * (1 - y*y) * (1 - z*z);
  real ddyy = -2 * sin(TIME) * (1 - x*x) * (1 - z*z);
  real ddzz = -2 * sin(TIME) * (1 - x*x) * (1 - y*y);
  real dt = cos(TIME) * (1 - x*x) * (1 - y*y) * (1 - z*z);

  return dt - (ddxx + ddyy + ddzz);
}

template<typename real, typename scalar>
scalar linear_form(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data) 
{
  return int_F_v<real, scalar>(n, wt, rhs, v, e) + int_u_v<real, scalar>(n, wt, data->ext + 0, v, e) / TAU;
}
