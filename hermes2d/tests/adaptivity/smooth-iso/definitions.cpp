#include "definitions.h"

void CustomExactSolution::derivatives(double x, double y, double& dx, double& dy) const
{
  dx = Hermes::cos(x)*Hermes::sin(y);
  dy = Hermes::sin(x)*Hermes::cos(y);
}

double CustomExactSolution::value(double x, double y) const
{
  return Hermes::sin(x)*Hermes::sin(y);
}

Ord CustomExactSolution::ord(Ord x, Ord y) const
{
  return Ord(7);
}

double CustomFunction::value(double x, double y) const
{
  return -2*Hermes::sin(x)*Hermes::sin(y);
}

Ord CustomFunction::value(Ord x, Ord y) const
{
  return Ord(7);
}