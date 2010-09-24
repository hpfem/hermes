static double fn(double x, double y)
{
  double r = sqrt(x*x + y*y);
  return sin(1/(ALPHA + r));
}

static double fndd(double x, double y, double& dx, double& dy)
{
  double r = sqrt(x*x + y*y);
  double h = 1/(ALPHA + r);
  dx = -cos(h) * h * h * x / r;
  dy = -cos(h) * h * h * y / r;
  return fn(x, y);
}
