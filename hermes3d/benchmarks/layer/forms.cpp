template<typename real, typename scalar>
scalar biform(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *u, fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data)
{
  return int_grad_u_grad_v<real, scalar>(n, wt, u, v, e);
}

template<typename real>
real rhs(real x, real y, real z)
{
  real t2 = sqr(z + 0.25) + sqr(y + 0.25) + sqr(x + 0.25);
  real t = sqrt(t2);
  real u = sqr(SLOPE) * sqr(t - M_PI/3) + 1;
  real v = 2 * pow(SLOPE, 3) * (t - M_PI/3) / (t2 * sqr(u));
  real w = SLOPE / (pow(t2, 1.5) * u);

  return (3 * SLOPE) / (t * u) - t2 * (v + w);
}

template<typename real, typename scalar>
scalar liform(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data)
{
  return -int_F_v<real, scalar>(n, wt, rhs, v, e);
}


