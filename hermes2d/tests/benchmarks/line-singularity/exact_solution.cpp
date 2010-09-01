static double fn(double x, double y)
{
  if (x <= 0) return cos(K*y);
  else return cos(K*y) + pow(x, ALPHA);
}

static double fndd(double x, double y, double& dx, double& dy)
{
  if (x <= 0) dx = 0;
  else dx = ALPHA*pow(x, ALPHA - 1);
  dy = -sin(K*y)*K;
  return fn(x, y);
}
