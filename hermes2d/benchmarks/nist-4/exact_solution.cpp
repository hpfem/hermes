static double fn(double x, double y)
{
  return exp(-ALPHA_P * (pow((x - X_LOC), 2) + pow((y - Y_LOC), 2)));
}

static double fndd(double x, double y, double& dx, double& dy)
{
  double a = -ALPHA_P * ( (x - X_LOC) * (x - X_LOC) + (y - Y_LOC) * (y - Y_LOC));

  dx = -exp(a) * (2 * ALPHA_P * (x - X_LOC));
  dy = -exp(a) * (2 * ALPHA_P * (y - Y_LOC));

  return fn(x, y);
}
