static double fn(double x, double y)
{
  return exp(-ALPHA*(pow((x - X_LOC), 2) + pow((y - Y_LOC), 2)));
}

static double fndd(double x, double y, double& dx, double& dy)
{
  double a = (-ALPHA * pow((x - X_LOC), 2) - ALPHA * pow((y - Y_LOC), 2));
  double b = (2 * ALPHA * x - ALPHA);
  double c = (2 * ALPHA * y - ALPHA);

  dx = -exp(a)*b;
  dy = -exp(a)*c;

  return fn(x, y);
}
