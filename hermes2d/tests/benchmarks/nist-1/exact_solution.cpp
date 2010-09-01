static double fn(double x, double y)
{
  return pow(2, 4 * EXACT_SOL_P) * pow(x, EXACT_SOL_P) * pow(1 - x, EXACT_SOL_P) * pow(y, EXACT_SOL_P) * pow(1 - y, EXACT_SOL_P);
}

static double fndd(double x, double y, double& dx, double& dy)
{
  double A = pow((1.0-y), EXACT_SOL_P);
  double B = pow((1.0-x), EXACT_SOL_P);
  double C = pow(y, EXACT_SOL_P);
  double D = pow(x, EXACT_SOL_P);

  dx = ((EXACT_SOL_P * pow(16.0, EXACT_SOL_P)*A*C) / (x-1.0) + (EXACT_SOL_P*pow(16.0, EXACT_SOL_P)*A*C)/x)*B*D;
  dy = ((EXACT_SOL_P*pow(16.0, EXACT_SOL_P)*B*D)/(y-1.0)+(EXACT_SOL_P*pow(16.0, EXACT_SOL_P)*B*D)/y)*A*C;

  return fn(x, y);
}
