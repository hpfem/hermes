#include "definitions.h"

double CustomExactSolution::value (double x, double y) const 
{
  return Hermes::atan(slope * (Hermes::sqrt(Hermes::sqr(x-1.25) + Hermes::sqr(y+0.25)) - M_PI/3));
}

void CustomExactSolution::derivatives (double x, double y, double& dx, double& dy) const 
{
  double t = Hermes::sqrt(Hermes::sqr(x - 1.25) + Hermes::sqr(y + 0.25));
  double u = t * (Hermes::sqr(slope) * Hermes::sqr(t - M_PI/3) + 1);
  dx = slope * (x - 1.25) / u;
  dy = slope * (y + 0.25) / u;
}

Ord CustomExactSolution::ord(double x, double y) const 
{
  return Ord(20);
}


double CustomFunction::value(double x, double y) const 
{
  double t2 = Hermes::sqr(y + 0.25) + Hermes::sqr(x - 1.25);
  double t = Hermes::sqrt(t2);
  double u = (Hermes::sqr(M_PI - 3.0*t) * Hermes::sqr(slope) + 9.0);

  return (   27.0/2.0 * Hermes::sqr(2.0*y + 0.5) * (M_PI - 3.0*t) * Hermes::pow(slope, 3.0) / (Hermes::sqr(u) * t2) 
           + 27.0/2.0 * Hermes::sqr(2.0*x - 2.5) * (M_PI - 3.0*t) * Hermes::pow(slope, 3.0) / (Hermes::sqr(u) * t2) 
           - 9.0/4.0 * Hermes::sqr(2.0*y + 0.5) * slope / (u * Hermes::pow(t,3.0)) 
           - 9.0/4.0 * Hermes::sqr(2.0*x - 2.5) * slope / (u * Hermes::pow(t,3.0)) + 18.0 * slope / (u * t)
		);
}

Ord CustomFunction::value(Ord x, Ord y) const 
{
  return Ord(20);
}
