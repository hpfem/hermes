// Constant function.
static double fn_const(double x, double y, double& dx, double& dy)
{
  dx = dy = 0;
  return 3;
}

// Function of a plane that rotates into the screen plane.
static double fn_plane(double x, double y, double& dx, double& dy)
{
  const double k = tan(50*M_PI/180.);
  dx = dy = 0; // Don't bother with derivatives, we just want to display the function.
  return k*y;
}

// Function that looks almost like a cuboid.
// (A Cartesian product of exact solutions to the 1D problem -u'' + K*K*u = K*K
// in (-1,1) with zero Dirichlet BC.)
static double fn_cuboid(double x, double y, double& dx, double& dy)
{
  static double K = 1e2;
  dx = dy = 0;
  return  (1. - (exp(K*x) + exp(-K*x))/(exp(K) + exp(-K))) *
          (1. - (exp(K*y) + exp(-K*y))/(exp(K) + exp(-K)));
}


// Function that for (x,y) in (-1,-1)x(1,1) spans a large range of values. It
// requires manual setting of limits on the vertical axis to be displayed reasonably.
static double fn_bigrange(double x, double y, double& dx, double& dy)
{
  dx = dy = 0;
  return cos(4*sqr(M_PI)*(sqr(x)+sqr(y)))/exp(-4*(sqr(x)+sqr(y)));
}

// Paraboloid.
static double fn_paraboloid(double x, double y, double& dx, double& dy)
{
  dx = dy = 0;
  return sqr(x) + sqr(y);
}

