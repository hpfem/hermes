double fn(double x, double y, double z)
{
  return sin(TIME) * (1 - x*x) * (1 - y*y) * (1 - z*z);
}

double fndd(double x, double y, double z, double &dx, double &dy, double &dz) 
{
  dx = -2 * sin(TIME) * x * (1 - y*y) * (1 - z*z);
  dy = -2 * sin(TIME) * (1 - x*x) * y * (1 - z*z);
  dz = -2 * sin(TIME) * (1 - x*x) * (1 - y*y) * z;

  return fn(x, y, z);
}
