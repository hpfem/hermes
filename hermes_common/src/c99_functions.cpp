#include <cmath>
#include "common.h"
#include "c99_functions.h"

#ifdef IMPLEMENT_C99

double Hermes::C_99::exp2(double x)
{
  return pow(2.0, x);
}

double Hermes::C_99::log2(double x)
{
  return log(x) / M_LN2;
}

double Hermes::C_99::cbrt(double x)
{
  if (!_isnan(x))
    return pow(x, 1.0 / 3.0);
  else
    return x;
}

#endif
