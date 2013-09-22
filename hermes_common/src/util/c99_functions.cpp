#include <cmath>
#include "common.h"
#include "c99_functions.h"
/*! \file c99_functions.cpp
\brief File containing definitions from the C99 standard that are missing in MSVC.
*/
#ifdef IMPLEMENT_C99

double exp2(double x)
{
  return pow(2.0, x);
}

double log2(double x)
{
  return log(x) / M_LN2;
}

double cbrt(double x)
{
  if(!_isnan(x))
    return pow(x, 1.0 / 3.0);
  else
    return x;
}

#endif