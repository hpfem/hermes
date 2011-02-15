static double fn(double x, double y)
{
  return atan(K * x);
}

static double exact(double x, double y, double& dx, double& dy)
{
  dx = 1./(1 + (K*x)*(K*x)) * K;
  dy = 0;
  return fn(x, y);
}
