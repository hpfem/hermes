double fn(double x, double y, double z) 
{
  return pow(x*x + y*y + z*z, .25);
}

double fndd(double x, double y, double z, double &dx, double &dy, double &dz) {
  dx = 0.5 * x * pow(x*x + y*y + z*z, -.75);
  dy = 0.5 * y * pow(x*x + y*y + z*z, -.75);
  dz = 0.5 * z * pow(x*x + y*y + z*z, -.75);

  return fn(x, y, z);
}
