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

