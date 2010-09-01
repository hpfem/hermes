static double fn(double x, double y)
{
  return pow(x, ALPHA);
}

static double fndd(double x, double y, double& dx, double& dy)
{
  dx = (ALPHA/(pow(x, 0.4)));
  dy = 0;

  return fn(x, y);
}
