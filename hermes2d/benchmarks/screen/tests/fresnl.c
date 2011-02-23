/*							fresnl.c
 *
 *	Fresnel integral
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, S, C;
 * void fresnl();
 *
 * fresnl( x, _&S, _&C );
 *
 *
 * DESCRIPTION:
 *
 * Evaluates the Fresnel integrals
 *
 *           x
 *           -
 *          | |
 * C(x) =   |   cos(pi/2 t**2) dt,
 *        | |
 *         -
 *          0
 *
 *           x
 *           -
 *          | |
 * S(x) =   |   sin(pi/2 t**2) dt.
 *        | |
 *         -
 *          0
 *
 *
 * The integrals are evaluated by power series for small x.
 * For x >= 1 auxiliary functions f(x) and g(x) are employed
 * such that
 *
 * C(x) = 0.5 + f(x) sin( pi/2 x**2 ) - g(x) cos( pi/2 x**2 )
 * S(x) = 0.5 - f(x) cos( pi/2 x**2 ) - g(x) sin( pi/2 x**2 )
 *
 *
 *
 * ACCURACY:
 *
 *  Relative error.
 *
 * Arithmetic  function   domain     # trials      peak         rms
 *   IEEE       S(x)      0, 10       30000       1.1e-6      1.9e-7
 *   IEEE       C(x)      0, 10       30000       1.1e-6      2.0e-7
 */

/*
Cephes Math Library Release 2.1:  January, 1989
Copyright 1984, 1987, 1989 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <math.h>


/* S(x) for small x */
static double sn[7] = {
 1.647629463788700E-009,
-1.522754752581096E-007,
 8.424748808502400E-006,
-3.120693124703272E-004,
 7.244727626597022E-003,
-9.228055941124598E-002,
 5.235987735681432E-001
};

/* C(x) for small x */
static double cn[7] = {
 1.416802502367354E-008,
-1.157231412229871E-006,
 5.387223446683264E-005,
-1.604381798862293E-003,
 2.818489036795073E-002,
-2.467398198317899E-001,
 9.999999760004487E-001
};


/* Auxiliary function f(x) */
static double fn[8] = {
-1.903009855649792E+012,
 1.355942388050252E+011,
-4.158143148511033E+009,
 7.343848463587323E+007,
-8.732356681548485E+005,
 8.560515466275470E+003,
-1.032877601091159E+002,
 2.999401847870011E+000
};

/* Auxiliary function g(x) */
static double gn[8] = {
-1.860843997624650E+011,
 1.278350673393208E+010,
-3.779387713202229E+008,
 6.492611570598858E+006,
-7.787789623358162E+004,
 8.602931494734327E+002,
-1.493439396592284E+001,
 9.999841934744914E-001
};


#define PIF M_PI
#define PIO2F PIF/2.0



double polevl( double xx, double *coef, int N )
{
  double ans, x;
  double *p;
  int i;

  x = xx;
  p = coef;
  ans = *p++;

  i = N;
  do
        ans = ans * x  +  *p++;
  while( --i );

  return ans;
}




void fresnl( double xxa, double *ssa, double *cca )
{
  double f, g, cc, ss, c, s, t, u, x, x2;
  /*debug double t1;*/

  x = xxa;
  x = fabs(x);
  x2 = x * x;
  if( x2 < 2.5625 )
	{
	t = x2 * x2;
	ss = x * x2 * polevl( t, sn, 6);
	cc = x * polevl( t, cn, 6);
	goto done;
	}

  if( x > 36974.0 )
	{
	cc = 0.5;
	ss = 0.5;
	goto done;
	}

/*		Asymptotic power series auxiliary functions
 *		for large argument
 */
	x2 = x * x;
	t = PIF * x2;
	u = 1.0/(t * t);
	t = 1.0/t;
	f = 1.0 - u * polevl( u, fn, 7);
	g = t * polevl( u, gn, 7);

	t = PIO2F * x2;
	c = cos(t);
	s = sin(t);
	t = PIF * x;
	cc = 0.5  +  (f * s  -  g * c)/t;
	ss = 0.5  -  (f * c  +  g * s)/t;

done:
  if( xxa < 0.0 )
	{
	cc = -cc;
	ss = -ss;
	}

  *cca = cc;
  *ssa = ss;
}



/*#include <stdio.h>

int main(int argc, char* argv[])
{
  double c, s;
  fresnl(1, &s, &c);
  printf("c=%18e s=%18e\n", c, s);
  return 0;
}
*/


