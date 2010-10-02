static double fn(double x, double y)
{
  return (1 - exp(-(1-x)/EPSILON)) * (1 - exp(-(1-y)/EPSILON))
    * cos(M_PI * (x + y)); 
}

// SYMPY
static double fndd(double x, double y, double& dx, double& dy)
{
  dx = -M_PI*(1 - exp(-(1 - x)/EPSILON))*(1 - exp(-(1 - y)/EPSILON))*sin(M_PI*(x + y)) - (1 - exp(-(1 - y)/EPSILON))*cos(M_PI*(x + y))*exp(-(1 - x)/EPSILON)/EPSILON;
  dy = -M_PI*(1 - exp(-(1 - x)/EPSILON))*(1 - exp(-(1 - y)/EPSILON))*sin(M_PI*(x + y)) - (1 - exp(-(1 - x)/EPSILON))*cos(M_PI*(x + y))*exp(-(1 - y)/EPSILON)/EPSILON;
  return fn(x, y);
}

// Right-hand side.
template<typename Real>
Real rhs(Real x, Real y)
{
  return -EPSILON*(-2*pow(M_PI,2)*(1 - exp(-(1 - x)/EPSILON))*(1 - exp(-(1 - y)/EPSILON))*cos(M_PI*(x + y)) + 2*M_PI*(1 - exp(-(1 - x)/EPSILON))*exp(-(1 - y)/EPSILON)*sin(M_PI*(x + y))/EPSILON + 2*M_PI*(1 - exp(-(1 - y)/EPSILON))*exp(-(1 - x)/EPSILON)*sin(M_PI*(x + y))/EPSILON - (1 - exp(-(1 - y)/EPSILON))*cos(M_PI*(x + y))*exp(-(1 - x)/EPSILON)/pow(EPSILON,2) - (1 - exp(-(1 - x)/EPSILON))*cos(M_PI*(x + y))*exp(-(1 - y)/EPSILON)/pow(EPSILON,2)) - 3*M_PI*(1 - exp(-(1 - x)/EPSILON))*(1 - exp(-(1 - y)/EPSILON))*sin(M_PI*(x + y)) - 2*(1 - exp(-(1 - y)/EPSILON))*cos(M_PI*(x + y))*exp(-(1 - x)/EPSILON)/EPSILON - (1 - exp(-(1 - x)/EPSILON))*cos(M_PI*(x + y))*exp(-(1 - y)/EPSILON)/EPSILON;
}

/*
static double fndd(double x, double y, double& dx, double& dy)
{
  double a = exp(-(1-x)/EPSILON);
  double b = exp(-(1-y)/EPSILON);
  double c = cos(M_PI*(x+y));
  double d = sin(M_PI*(x+y));
  
  dx = -a*(1-b)*c/EPSILON - (1-a) * (1-b) * d;
  dy = -(1-a)*b*c/EPSILON - (1-a) * (1-b) * d;
  
  return fn(x, y);
}
*/
/*
// Right-hand side.
template<typename Real>
Real rhs(Real x, Real y) {
  Real a = exp(-(1-x)/EPSILON);
  Real b = exp(-(1-y)/EPSILON);
  Real c = cos(M_PI*(x+y));
  Real d = sin(M_PI*(x+y));
  
  return a*(1-b)*c/(EPSILON) - 2*a*(1-a)*(1-b)*d + 2*EPSILON*(1-a)*(1-b)*c*M_PI*M_PI 
         + (1-a)*b*c/EPSILON - 2*(1-a)*b*d - 2*a*(1-b)*c/EPSILON - 3*(1-a)*(1-b)*d*M_PI - (1-a)*b*c/EPSILON;
}
*/
