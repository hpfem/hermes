#include <cmath>
#include "../common.h"
#include "c99_functions.h"

#ifdef IMPLEMENT_C99

/* constants */
const unsigned long long _NAN = 0x7fffffffffffffffL; //NAN according to IEEE specification
H3D_API const double NAN = *(double*)&_NAN;

/* functions */
H3D_API double exp2(double x)
{
	return pow(2.0, x);
}

H3D_API double log2(double x)
{
	return log(x) / M_LN2;
}

H3D_API double cbrt(double x)
{
	if (!_isnan(x))
		return pow(x, 1.0 / 3.0);
	else
		return x;
}

#endif
