/* Exact solution */

double fn(double x, double y, double z)
{
  return atan(SLOPE * (sqrt(sqr(x + 0.25) + sqr(y + 0.25) + sqr(z + 0.25)) - M_PI/3));
}

double fndd(double x, double y, double z, double &dx, double &dy, double &dz)
{
  double t = sqrt(sqr(z + 0.25) + sqr(y + 0.25) + sqr(x + 0.25));
  double u = t * (sqr(SLOPE) * sqr(t - M_PI/3) + 1);

  dx = (SLOPE * (x + 0.25)) / u;
  dy = (SLOPE * (y + 0.25)) / u;
  dz = (SLOPE * (z + 0.25)) / u;

  return fn(x, y, z);
}

/* Weak forms */

template<typename real, typename scalar>
scalar biform(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, Func<real> *v, Geom<real> *e, ExtData<scalar> *data)
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
scalar liform(int n, double *wt, Func<scalar> *u_ext[], Func<real> *v, Geom<real> *e, ExtData<scalar> *data)
{
  return -int_F_v<real, scalar>(n, wt, rhs, v, e);
}


