#include <math.h>
#include <complex>

/****************************************************************
***  For the Bessel and Gamma functions,                      ***
***  copy the following into your cpp file.                   ***
***  (Taken from cmath packege from NetLib and modified.)     ***
*****************************************************************/

/*              mconf.h
*
*  Common include file for math routines
*
*
*
* SYNOPSIS:
*
* #include "mconf.h"
*
*
*
* DESCRIPTION:
*
* This file contains definitions for error codes that are
* passed to the common error handling routine mtherr()
* (which see).
*
* The file also includes a conditional assembly definition
* for the type of computer arithmetic (IEEE, DEC, Motorola
* IEEE, or UNKnown).
*
* For Digital Equipment PDP-11 and VAX computers, certain
* IBM systems, and others that use numbers with a 56-bit
* significand, the symbol DEC should be defined.  In this
* mode, most floating point constants are given as arrays
* of octal integers to eliminate decimal to binary conversion
* errors that might be introduced by the compiler.
*
* For little-endian computers, such as IBM PC, that follow the
* IEEE Standard for Binary Floating Point Arithmetic (ANSI/IEEE
* Std 754-1985), the symbol IBMPC should be defined.  These
* numbers have 53-bit significands.  In this mode, constants
* are provided as arrays of hexadecimal 16 bit integers.
*
* Big-endian IEEE format is denoted MIEEE.  On some RISC
* systems such as Sun SPARC, double precision constants
* must be stored on 8-byte address boundaries.  std::since integer
* arrays may be aligned differently, the MIEEE configuration
* may fail on such machines.
*
* To accommodate other types of computer arithmetic, all
* constants are also provided in a normal decimal radix
* which one can hope are correctly converted to a suitable
* format by the available C language compiler.  To invoke
* this mode, define the symbol UNK.
*
* An important difference among these modes is a predefined
* set of machine arithmetic constants for each.  The numbers
* MACHEP (the machine roundoff error), MAXNUM (largest number
* represented), and several other parameters are preset by
* the configuration symbol.  Check the file const.c to
* ensure that these values are correct for your computer.
*
* Configurations NANS, INFINITIES, MINUSZERO, and DENORMAL
* may fail on many systems.  Verify that they are supposed
* to work on your computer.
*/
/*
Cephes Math Library Release 2.3:  June, 1995
Copyright 1984, 1987, 1989, 1995 by Stephen L. Moshier
*/


/* Define if the `long double' type works.  */
#define HAVE_LONG_DOUBLE 1

/* Define as the return type of signal handlers (int or void).  */
#define RETSIGTYPE void

/* Define if you have the ANSI C header files.  */
#define STDC_HEADERS 1

/* Define if your processor stores words with the most significant
byte first (like Motorola and SPARC, unlike Intel and VAX).  */
/* #undef WORDS_BIGENDIAN */

/* Define if floating point words are bigendian.  */
/* #undef FLOAT_WORDS_BIGENDIAN */

/* The number of bytes in a int.  */
#define SIZEOF_INT 4

/* Define if you have the <string.h> header file.  */
#define HAVE_STRING_H 1

/* Name of package */
#define PACKAGE "cephes"

/* Version number of package */
#define VERSION "2.7"

/* Constant definitions for math error conditions
*/

#define DOMAIN    1 /* argument domain error */
#define SING    2 /* argument std::singularity */
#define OVERFLOW  3 /* overflow range error */
#define UNDERFLOW 4 /* underflow range error */
#define TLOSS   5 /* total loss of precision */
#define PLOSS   6 /* partial loss of precision */

#define EDOM    33
#define ERANGE    34
/* Complex numeral.  */
typedef struct
{
  double r;
  double i;
} cmplx;

#ifdef HAVE_LONG_DOUBLE
/* Long double complex numeral.  */
typedef struct
{
  long double r;
  long double i;
} cmplxl;
#endif


/* Type of computer arithmetic */

/* PDP-11, Pro350, VAX:
*/
/* #define DEC 1 */

/* Intel IEEE, low order words come first:
*/
/* #define IBMPC 1 */

/* Motorola IEEE, high order words come first
* (Sun 680x0 workstation):
*/
/* #define MIEEE 1 */

/* UNKnown arithmetic, invokes coefficients given in
* normal decimal format.  Beware of range boundary
* problems (MACHEP, MAXLOG, etc. in const.c) and
* roundoff problems in std::pow.c:
* (Sun SPARCstation)
*/
#define UNK 1

/* If you define UNK, then be sure to set BIGENDIAN properly. */
#ifdef FLOAT_WORDS_BIGENDIAN
#define BIGENDIAN 1
#else
#define BIGENDIAN 0
#endif
/* Define this `volatile' if your compiler thinks
* that floating point arithmetic obeys the associative
* and distributive laws.  It will defeat some optimizations
* (but probably not enough of them).
*
* #define VOLATILE volatile
*/
#define VOLATILE

/* For 12-byte long doubles on an i386, pad a 16-bit short 0
* to the end of real constants initialized by integer arrays.
*
* #define XPD 0,
*
* Otherwise, the type is 10 bytes long and XPD should be
* defined blank (e.g., Microsoft C).
*
* #define XPD
*/
#define XPD 0,

/* Define to support tiny denormal numbers, else undefine. */
#define DENORMAL 1

/* Define to ask for infinity support, else undefine. */
// #define INFINITIES 1

/* Define to ask for support of numbers that are Not-a-Number,
else undefine.  This may automatically define INFINITIES in some files. */
// #define NANS 1

/* Define to distinguish between -0.0 and +0.0.  */
#define MINUSZERO 1

/* Define 1 for ANSI C atan2() function
See atan.c and cstd::log.c. */
#define ANSIC 1

/* Get ANSI function prototypes, if you want them. */
#if 1
/* #ifdef __STDC__ */
#define ANSIPROT 1
int mtherr ( const char *, int );
#else
int mtherr();
#endif

/* Variable for error reporting.  See mtherr.c.  */
extern int merror;


///////////////////////// END of   mconf.h ///////////////////////










/*              const.c
*
*  Globally declared constants
*
*
*
* SYNOPSIS:
*
* extern double nameofconstant;
*
*
*
*
* DESCRIPTION:
*
* This file contains a number of mathematical constants and
* also some needed size parameters of the computer arithmetic.
* The values are supplied as arrays of hexadecimal integers
* for IEEE arithmetic; arrays of octal constants for DEC
* arithmetic; and in a normal decimal scientific notation for
* other machines.  The particular notation used is determined
* by a symbol (DEC, IBMPC, or UNK) defined in the include file
* mconf.h.
*
* The default size parameters are as follows.
*
* For DEC and UNK modes:
* MACHEP =  1.38777878078144567553E-17       2**-56
* MAXLOG =  8.8029691931113054295988E1       std::log(2**127)
* MINLOG = -8.872283911167299960540E1        std::log(2**-128)
* MAXNUM =  1.701411834604692317316873e38    2**127
*
* For IEEE arithmetic (IBMPC):
* MACHEP =  1.11022302462515654042E-16       2**-53
* MAXLOG =  7.09782712893383996843E2         std::log(2**1024)
* MINLOG = -7.08396418532264106224E2         std::log(2**-1022)
* MAXNUM =  1.7976931348623158E308           2**1024
*
* The global symbols for mathematical constants are
* PI     =  3.14159265358979323846           pi
* PIO2   =  1.57079632679489661923           pi/2
* PIO4   =  7.85398163397448309616E-1        pi/4
* SQRT2  =  1.41421356237309504880           std::sqrt(2)
* SQRTH  =  7.07106781186547524401E-1        std::sqrt(2)/2
* LOG2E  =  1.4426950408889634073599         1/std::log(2)
* SQ2OPI =  7.9788456080286535587989E-1      std::sqrt( 2/pi )
* LOGE2  =  6.93147180559945309417E-1        std::log(2)
* LOGSQ2 =  3.46573590279972654709E-1        std::log(2)/2
* THPIO4 =  2.35619449019234492885           3*pi/4
* TWOOPI =  6.36619772367581343075535E-1     2/pi
*
* These lists are subject to change.
*/

/*              const.c */

/*
Cephes Math Library Release 2.3:  March, 1995
Copyright 1984, 1995 by Stephen L. Moshier
*/

#ifdef UNK
#if 1
double MACHEP =  1.11022302462515654042E-16;   /* 2**-53 */
#else
double MACHEP =  1.38777878078144567553E-17;   /* 2**-56 */
#endif
double UFLOWTHRESH =  2.22507385850720138309E-308; /* 2**-1022 */
#ifdef DENORMAL
double MAXLOG =  7.09782712893383996732E2;     /* std::log(MAXNUM) */
/* double MINLOG = -7.44440071921381262314E2; */     /* std::log(2**-1074) */
double MINLOG = -7.451332191019412076235E2;     /* std::log(2**-1075) */
#else
double MAXLOG =  7.08396418532264106224E2;     /* std::log 2**1022 */
double MINLOG = -7.08396418532264106224E2;     /* std::log 2**-1022 */
#endif
double MAXNUM =  1.79769313486231570815E308;    /* 2**1024*(1-MACHEP) */
double PI     =  3.14159265358979323846;       /* pi */
double PIO2   =  1.57079632679489661923;       /* pi/2 */
double PIO4   =  7.85398163397448309616E-1;    /* pi/4 */
double SQRT2  =  1.41421356237309504880;       /* std::sqrt(2) */
double SQRTH  =  7.07106781186547524401E-1;    /* std::sqrt(2)/2 */
double LOG2E  =  1.4426950408889634073599;     /* 1/std::log(2) */
double SQ2OPI =  7.9788456080286535587989E-1;  /* std::sqrt( 2/pi ) */
double LOGE2  =  6.93147180559945309417E-1;    /* std::log(2) */
double LOGSQ2 =  3.46573590279972654709E-1;    /* std::log(2)/2 */
double THPIO4 =  2.35619449019234492885;       /* 3*pi/4 */
double TWOOPI =  6.36619772367581343075535E-1; /* 2/pi */
#ifdef INFINITIES
double INFINITY_BESSEL = 1.0/0.0;  /* 99e999; */
#else
double INFINITY_BESSEL =  1.79769313486231570815E308;    /* 2**1024*(1-MACHEP) */
#endif
#ifdef NANS
double NAN_BESSEL = 1.0/0.0 - 1.0/0.0;
#else
double NAN_BESSEL = 0.0;
#endif
#ifdef MINUSZERO
double NEGZERO = -0.0;
#else
double NEGZERO = 0.0;
#endif
#endif

#ifdef IBMPC
/* 2**-53 =  1.11022302462515654042E-16 */
unsigned short MACHEP[4] = {0x0000,0x0000,0x0000,0x3ca0};
unsigned short UFLOWTHRESH[4] = {0x0000,0x0000,0x0000,0x0010};
#ifdef DENORMAL
/* std::log(MAXNUM) =  7.09782712893383996732224E2 */
unsigned short MAXLOG[4] = {0x39ef,0xfefa,0x2e42,0x4086};
/* std::log(2**-1074) = - -7.44440071921381262314E2 */
/*unsigned short MINLOG[4] = {0x71c3,0x446d,0x4385,0xc087};*/
unsigned short MINLOG[4] = {0x3052,0xd52d,0x4910,0xc087};
#else
/* std::log(2**1022) =   7.08396418532264106224E2 */
unsigned short MAXLOG[4] = {0xbcd2,0xdd7a,0x232b,0x4086};
/* std::log(2**-1022) = - 7.08396418532264106224E2 */
unsigned short MINLOG[4] = {0xbcd2,0xdd7a,0x232b,0xc086};
#endif
/* 2**1024*(1-MACHEP) =  1.7976931348623158E308 */
unsigned short MAXNUM[4] = {0xffff,0xffff,0xffff,0x7fef};
unsigned short PI[4]     = {0x2d18,0x5444,0x21fb,0x4009};
unsigned short PIO2[4]   = {0x2d18,0x5444,0x21fb,0x3ff9};
unsigned short PIO4[4]   = {0x2d18,0x5444,0x21fb,0x3fe9};
unsigned short SQRT2[4]  = {0x3bcd,0x667f,0xa09e,0x3ff6};
unsigned short SQRTH[4]  = {0x3bcd,0x667f,0xa09e,0x3fe6};
unsigned short LOG2E[4]  = {0x82fe,0x652b,0x1547,0x3ff7};
unsigned short SQ2OPI[4] = {0x3651,0x33d4,0x8845,0x3fe9};
unsigned short LOGE2[4]  = {0x39ef,0xfefa,0x2e42,0x3fe6};
unsigned short LOGSQ2[4] = {0x39ef,0xfefa,0x2e42,0x3fd6};
unsigned short THPIO4[4] = {0x21d2,0x7f33,0xd97c,0x4002};
unsigned short TWOOPI[4] = {0xc883,0x6dc9,0x5f30,0x3fe4};
#ifdef INFINITIES
unsigned short INFINITY_BESSEL[4] = {0x0000,0x0000,0x0000,0x7ff0};
#else
unsigned short INFINITY_BESSEL[4] = {0xffff,0xffff,0xffff,0x7fef};
#endif
#ifdef NANS
unsigned short NAN_BESSEL[4] = {0x0000,0x0000,0x0000,0x7ffc};
#else
unsigned short NAN_BESSEL[4] = {0x0000,0x0000,0x0000,0x0000};
#endif
#ifdef MINUSZERO
unsigned short NEGZERO[4] = {0x0000,0x0000,0x0000,0x8000};
#else
unsigned short NEGZERO[4] = {0x0000,0x0000,0x0000,0x0000};
#endif
#endif

#ifdef MIEEE
/* 2**-53 =  1.11022302462515654042E-16 */
unsigned short MACHEP[4] = {0x3ca0,0x0000,0x0000,0x0000};
unsigned short UFLOWTHRESH[4] = {0x0010,0x0000,0x0000,0x0000};
#ifdef DENORMAL
/* std::log(2**1024) =   7.09782712893383996843E2 */
unsigned short MAXLOG[4] = {0x4086,0x2e42,0xfefa,0x39ef};
/* std::log(2**-1074) = - -7.44440071921381262314E2 */
/* unsigned short MINLOG[4] = {0xc087,0x4385,0x446d,0x71c3}; */
unsigned short MINLOG[4] = {0xc087,0x4910,0xd52d,0x3052};
#else
/* std::log(2**1022) =  7.08396418532264106224E2 */
unsigned short MAXLOG[4] = {0x4086,0x232b,0xdd7a,0xbcd2};
/* std::log(2**-1022) = - 7.08396418532264106224E2 */
unsigned short MINLOG[4] = {0xc086,0x232b,0xdd7a,0xbcd2};
#endif
/* 2**1024*(1-MACHEP) =  1.7976931348623158E308 */
unsigned short MAXNUM[4] = {0x7fef,0xffff,0xffff,0xffff};
unsigned short PI[4]     = {0x4009,0x21fb,0x5444,0x2d18};
unsigned short PIO2[4]   = {0x3ff9,0x21fb,0x5444,0x2d18};
unsigned short PIO4[4]   = {0x3fe9,0x21fb,0x5444,0x2d18};
unsigned short SQRT2[4]  = {0x3ff6,0xa09e,0x667f,0x3bcd};
unsigned short SQRTH[4]  = {0x3fe6,0xa09e,0x667f,0x3bcd};
unsigned short LOG2E[4]  = {0x3ff7,0x1547,0x652b,0x82fe};
unsigned short SQ2OPI[4] = {0x3fe9,0x8845,0x33d4,0x3651};
unsigned short LOGE2[4]  = {0x3fe6,0x2e42,0xfefa,0x39ef};
unsigned short LOGSQ2[4] = {0x3fd6,0x2e42,0xfefa,0x39ef};
unsigned short THPIO4[4] = {0x4002,0xd97c,0x7f33,0x21d2};
unsigned short TWOOPI[4] = {0x3fe4,0x5f30,0x6dc9,0xc883};
#ifdef INFINITIES
unsigned short INFINITY_BESSEL[4] = {0x7ff0,0x0000,0x0000,0x0000};
#else
unsigned short INFINITY_BESSEL[4] = {0x7fef,0xffff,0xffff,0xffff};
#endif
#ifdef NANS
unsigned short NAN_BESSEL[4] = {0x7ff8,0x0000,0x0000,0x0000};
#else
unsigned short NAN_BESSEL[4] = {0x0000,0x0000,0x0000,0x0000};
#endif
#ifdef MINUSZERO
unsigned short NEGZERO[4] = {0x8000,0x0000,0x0000,0x0000};
#else
unsigned short NEGZERO[4] = {0x0000,0x0000,0x0000,0x0000};
#endif
#endif

#ifdef DEC
/* 2**-56 =  1.38777878078144567553E-17 */
unsigned short MACHEP[4] = {0022200,0000000,0000000,0000000};
unsigned short UFLOWTHRESH[4] = {0x0080,0x0000,0x0000,0x0000};
/* std::log 2**127 = 88.029691931113054295988 */
unsigned short MAXLOG[4] = {041660,007463,0143742,025733,};
/* std::log 2**-128 = -88.72283911167299960540 */
unsigned short MINLOG[4] = {0141661,071027,0173721,0147572,};
/* 2**127 = 1.701411834604692317316873e38 */
unsigned short MAXNUM[4] = {077777,0177777,0177777,0177777,};
unsigned short PI[4]     = {040511,007732,0121041,064302,};
unsigned short PIO2[4]   = {040311,007732,0121041,064302,};
unsigned short PIO4[4]   = {040111,007732,0121041,064302,};
unsigned short SQRT2[4]  = {040265,002363,031771,0157145,};
unsigned short SQRTH[4]  = {040065,002363,031771,0157144,};
unsigned short LOG2E[4]  = {040270,0125073,024534,013761,};
unsigned short SQ2OPI[4] = {040114,041051,0117241,0131204,};
unsigned short LOGE2[4]  = {040061,071027,0173721,0147572,};
unsigned short LOGSQ2[4] = {037661,071027,0173721,0147572,};
unsigned short THPIO4[4] = {040426,0145743,0174631,007222,};
unsigned short TWOOPI[4] = {040042,0174603,067116,042025,};
/* Approximate infinity by MAXNUM.  */
unsigned short INFINITY_BESSEL[4] = {077777,0177777,0177777,0177777,};
unsigned short NAN_BESSEL[4] = {0000000,0000000,0000000,0000000};
#ifdef MINUSZERO
unsigned short NEGZERO[4] = {0000000,0000000,0000000,0100000};
#else
unsigned short NEGZERO[4] = {0000000,0000000,0000000,0000000};
#endif
#endif

#ifndef UNK
extern unsigned short MACHEP[];
extern unsigned short UFLOWTHRESH[];
extern unsigned short MAXLOG[];
extern unsigned short UNDstd::log[];
extern unsigned short MINLOG[];
extern unsigned short MAXNUM[];
extern unsigned short PI[];
extern unsigned short PIO2[];
extern unsigned short PIO4[];
extern unsigned short SQRT2[];
extern unsigned short SQRTH[];
extern unsigned short LOG2E[];
extern unsigned short SQ2OPI[];
extern unsigned short LOGE2[];
extern unsigned short LOGSQ2[];
extern unsigned short THPIO4[];
extern unsigned short TWOOPI[];
extern unsigned short INFINITY_BESSEL[];
extern unsigned short NAN_BESSEL[];
extern unsigned short NEGZERO[];
#endif


//////////////////////////// end of    const.c ///////////////////////



/*              mtherr.c
*
*  Library common error handling routine
*
*
*
* SYNOPSIS:
*
* char *fctnam;
* int code;
* int mtherr();
*
* mtherr( fctnam, code );
*
*
*
* DESCRIPTION:
*
* This routine may be called to report one of the following
* error conditions (in the include file mconf.h).
*
*   Mnemonic        Value          Significance
*
*    DOMAIN            1       argument domain error
*    SING              2       function std::singularity
*    OVERFLOW          3       overflow range error
*    UNDERFLOW         4       underflow range error
*    TLOSS             5       total loss of precision
*    PLOSS             6       partial loss of precision
*    EDOM             33       Unix domain error code
*    ERANGE           34       Unix range error code
*
* The default version of the file prints the function name,
* passed to it by the pointer fctnam, followed by the
* error condition.  The display is directed to the standard
* output device.  The routine then returns to the calling
* program.  Users may wish to modify the program to abort by
* calling exit() under severe error conditions such as domain
* errors.
*
* std::since all error conditions pass control to this function,
* the display may be easily changed, eliminated, or directed
* to an error std::logging device.
*
* SEE ALSO:
*
* mconf.h
*
*/

/*
Cephes Math Library Release 2.0:  April, 1987
Copyright 1984, 1987 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <stdio.h>

int merror = 0;

/* Notice: the order of appearance of the following
* messages is bound to the error codes defined
* in mconf.h.
*/
static const char *ermsg[7] = {
  "unknown",      /* error code 0 */
  "domain",       /* error code 1 */
  "std::singularity",  /* et seq.      */
  "overflow",
  "underflow",
  "total loss of precision",
  "partial loss of precision"
};


int mtherr(const char *name, int code)
{

  /* Display string passed by calling program,
  * which is supposed to be the name of the
  * function in which the error occurred:
  */
  printf( "\n%s ", name );

  /* Set global error message word */
  merror = code;

  /* Display error message defined
  * by the code argument.
  */
  if( (code <= 0) || (code >= 7) )
    code = 0;
  printf( "%s error\n", ermsg[code] );

  /* Return to calling
  * program
  */
  return( 0 );
}


//////////////////////////// end of    mtherr.c ///////////////////////


/*              polevl.c
*              p1evl.c
*
*  Evaluate polynomial
*
*
*
* SYNOPSIS:
*
* int N;
* double x, y, coef[N+1], polevl[];
*
* y = polevl( x, coef, N );
*
*
*
* DESCRIPTION:
*
* Evaluates polynomial of degree N:
*
*                     2          N
* y  =  C  + C x + C x  +...+ C x
*        0    1     2          N
*
* Coefficients are stored in reverse order:
*
* coef[0] = C  , ..., coef[N] = C  .
*            N                   0
*
*  The function p1evl() assumes that coef[N] = 1.0 and is
* omitted from the array.  Its calling arguments are
* otherwise the same as polevl().
*
*
* SPEED:
*
* In the interest of speed, there are no checks for out
* of bounds arithmetic.  This routine is used by most of
* the functions in the library.  Depending on available
* equipment features, the user may wish to rewrite the
* program in microcode or assembly language.
*
*/


/*
Cephes Math Library Release 2.1:  December, 1988
Copyright 1984, 1987, 1988 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/


double polevl( double x, double coef[], int N)
{
  double ans;
  int i;
  double *p;

  p = coef;
  ans = *p++;
  i = N;

  do
  ans = ans * x  +  *p++;
  while( --i );

  return( ans );
}

/*              p1evl() */
/*                                          N
* Evaluate polynomial when coefficient of x  is 1.0.
* Otherwise same as polevl.
*/

double p1evl( double x, double coef[], int N )
{
  double ans;
  double *p;
  int i;

  p = coef;
  ans = x + *p++;
  i = N-1;

  do
  ans = ans * x  + *p++;
  while( --i );

  return( ans );
}


//////////////////////////// end of    polevl.c ///////////////////////






////////////////////////////   airy.c      ///////////////////////

static double c1 = 0.35502805388781723926;
static double c2 = 0.258819403792806798405;
static double sqrt3 = 1.732050807568877293527;
static double sqpii = 5.64189583547756286948E-1;
extern double PI;

extern double MAXNUM, MACHEP;
#ifdef UNK
#define MAXAIRY 25.77
#endif
#ifdef DEC
#define MAXAIRY 25.77
#endif
#ifdef IBMPC
#define MAXAIRY 103.892
#endif
#ifdef MIEEE
#define MAXAIRY 103.892
#endif


#ifdef UNK
static double AN[8] = {
  3.46538101525629032477E-1,
  1.20075952739645805542E1,
  7.62796053615234516538E1,
  1.68089224934630576269E2,
  1.59756391350164413639E2,
  7.05360906840444183113E1,
  1.40264691163389668864E1,
  9.99999999999999995305E-1,
};
static double AD[8] = {
  5.67594532638770212846E-1,
  1.47562562584847203173E1,
  8.45138970141474626562E1,
  1.77318088145400459522E2,
  1.64234692871529701831E2,
  7.14778400825575695274E1,
  1.40959135607834029598E1,
  1.00000000000000000470E0,
};
#endif
#ifdef DEC
static unsigned short AN[32] = {
  0037661,0066561,0024675,0131301,
  0041100,0017434,0034324,0101466,
  0041630,0107450,0067427,0007430,
  0042050,0013327,0071000,0034737,
  0042037,0140642,0156417,0167366,
  0041615,0011172,0075147,0051165,
  0041140,0066152,0160520,0075146,
  0040200,0000000,0000000,0000000,
};
static unsigned short AD[32] = {
  0040021,0046740,0011422,0064606,
  0041154,0014640,0024631,0062450,
  0041651,0003435,0101152,0106401,
  0042061,0050556,0034605,0136602,
  0042044,0036024,0152377,0151414,
  0041616,0172247,0072216,0115374,
  0041141,0104334,0124154,0166007,
  0040200,0000000,0000000,0000000,
};
#endif
#ifdef IBMPC
static unsigned short AN[32] = {
  0xb658,0x2537,0x2dae,0x3fd6,
  0x9067,0x871a,0x03e3,0x4028,
  0xe1e3,0x0de2,0x11e5,0x4053,
  0x073c,0xee40,0x02da,0x4065,
  0xfddf,0x5ba1,0xf834,0x4063,
  0xea4f,0x4f4c,0xa24f,0x4051,
  0x0f4d,0x5c2a,0x0d8d,0x402c,
  0x0000,0x0000,0x0000,0x3ff0,
};
static unsigned short AD[32] = {
  0x4d31,0x0262,0x29bc,0x3fe2,
  0x2ca5,0x0533,0x8334,0x402d,
  0x51a0,0xb04d,0x20e3,0x4055,
  0xb7b0,0xc730,0x2a2d,0x4066,
  0xfa61,0x9a9f,0x8782,0x4064,
  0xd35f,0xee91,0xde94,0x4051,
  0x9d81,0x950d,0x311b,0x402c,
  0x0000,0x0000,0x0000,0x3ff0,
};
#endif
#ifdef MIEEE
static unsigned short AN[32] = {
  0x3fd6,0x2dae,0x2537,0xb658,
  0x4028,0x03e3,0x871a,0x9067,
  0x4053,0x11e5,0x0de2,0xe1e3,
  0x4065,0x02da,0xee40,0x073c,
  0x4063,0xf834,0x5ba1,0xfddf,
  0x4051,0xa24f,0x4f4c,0xea4f,
  0x402c,0x0d8d,0x5c2a,0x0f4d,
  0x3ff0,0x0000,0x0000,0x0000,
};
static unsigned short AD[32] = {
  0x3fe2,0x29bc,0x0262,0x4d31,
  0x402d,0x8334,0x0533,0x2ca5,
  0x4055,0x20e3,0xb04d,0x51a0,
  0x4066,0x2a2d,0xc730,0xb7b0,
  0x4064,0x8782,0x9a9f,0xfa61,
  0x4051,0xde94,0xee91,0xd35f,
  0x402c,0x311b,0x950d,0x9d81,
  0x3ff0,0x0000,0x0000,0x0000,
};
#endif

#ifdef UNK
static double APN[8] = {
  6.13759184814035759225E-1,
  1.47454670787755323881E1,
  8.20584123476060982430E1,
  1.71184781360976385540E2,
  1.59317847137141783523E2,
  6.99778599330103016170E1,
  1.39470856980481566958E1,
  1.00000000000000000550E0,
};
static double APD[8] = {
  3.34203677749736953049E-1,
  1.11810297306158156705E1,
  7.11727352147859965283E1,
  1.58778084372838313640E2,
  1.53206427475809220834E2,
  6.86752304592780337944E1,
  1.38498634758259442477E1,
  9.99999999999999994502E-1,
};
#endif
#ifdef DEC
static unsigned short APN[32] = {
  0040035,0017522,0065145,0054755,
  0041153,0166556,0161471,0057174,
  0041644,0016750,0034445,0046462,
  0042053,0027515,0152316,0046717,
  0042037,0050536,0067023,0023264,
  0041613,0172252,0007240,0131055,
  0041137,0023503,0052472,0002305,
  0040200,0000000,0000000,0000000,
};
static unsigned short APD[32] = {
  0037653,0016276,0112106,0126625,
  0041062,0162577,0067111,0111761,
  0041616,0054160,0140004,0137455,
  0042036,0143460,0104626,0157206,
  0042031,0032330,0067131,0114260,
  0041611,0054667,0147207,0134564,
  0041135,0114412,0070653,0146015,
  0040200,0000000,0000000,0000000,
};
#endif
#ifdef IBMPC
static unsigned short APN[32] = {
  0xab3e,0x4d4c,0xa3ea,0x3fe3,
  0x2bcf,0xdc67,0x7dad,0x402d,
  0xa9a6,0x0724,0x83bd,0x4054,
  0xc9ba,0xba99,0x65e9,0x4065,
  0x64d7,0xcdc2,0xea2b,0x4063,
  0x1646,0x41d4,0x7e95,0x4051,
  0x4099,0x6aa7,0xe4e8,0x402b,
  0x0000,0x0000,0x0000,0x3ff0,
};
static unsigned short APD[32] = {
  0xd5b3,0xd288,0x6397,0x3fd5,
  0x327e,0xedc9,0x5caf,0x4026,
  0x97e6,0x1800,0xcb0e,0x4051,
  0xdbd1,0x1132,0xd8e6,0x4063,
  0x3316,0x0dcb,0x269b,0x4063,
  0xf72f,0xf9d0,0x2b36,0x4051,
  0x7982,0x4e35,0xb321,0x402b,
  0x0000,0x0000,0x0000,0x3ff0,
};
#endif
#ifdef MIEEE
static unsigned short APN[32] = {
  0x3fe3,0xa3ea,0x4d4c,0xab3e,
  0x402d,0x7dad,0xdc67,0x2bcf,
  0x4054,0x83bd,0x0724,0xa9a6,
  0x4065,0x65e9,0xba99,0xc9ba,
  0x4063,0xea2b,0xcdc2,0x64d7,
  0x4051,0x7e95,0x41d4,0x1646,
  0x402b,0xe4e8,0x6aa7,0x4099,
  0x3ff0,0x0000,0x0000,0x0000,
};
static unsigned short APD[32] = {
  0x3fd5,0x6397,0xd288,0xd5b3,
  0x4026,0x5caf,0xedc9,0x327e,
  0x4051,0xcb0e,0x1800,0x97e6,
  0x4063,0xd8e6,0x1132,0xdbd1,
  0x4063,0x269b,0x0dcb,0x3316,
  0x4051,0x2b36,0xf9d0,0xf72f,
  0x402b,0xb321,0x4e35,0x7982,
  0x3ff0,0x0000,0x0000,0x0000,
};
#endif

#ifdef UNK
static double BN16[5] = {
  -2.53240795869364152689E-1,
  5.75285167332467384228E-1,
  -3.29907036873225371650E-1,
  6.44404068948199951727E-2,
  -3.82519546641336734394E-3,
};
static double BD16[5] = {
  /* 1.00000000000000000000E0,*/
  -7.15685095054035237902E0,
  1.06039580715664694291E1,
  -5.23246636471251500874E0,
  9.57395864378383833152E-1,
  -5.50828147163549611107E-2,
};
#endif
#ifdef DEC
static unsigned short BN16[20] = {
  0137601,0124307,0010213,0035210,
  0040023,0042743,0101621,0016031,
  0137650,0164623,0036056,0074511,
  0037203,0174525,0000473,0142474,
  0136172,0130041,0066726,0064324,
};
static unsigned short BD16[20] = {
  /*0040200,0000000,0000000,0000000,*/
  0140745,0002354,0044335,0055276,
  0041051,0124717,0170130,0104013,
  0140647,0070135,0046473,0103501,
  0040165,0013745,0033324,0127766,
  0137141,0117204,0076164,0033107,
};
#endif
#ifdef IBMPC
static unsigned short BN16[20] = {
  0x6751,0xe211,0x3518,0xbfd0,
  0x2383,0x7072,0x68bc,0x3fe2,
  0xcf29,0x6785,0x1d32,0xbfd5,
  0x78a8,0xa027,0x7f2a,0x3fb0,
  0xcd1b,0x2dba,0x5604,0xbf6f,
};
static unsigned short BD16[20] = {
  /*0x0000,0x0000,0x0000,0x3ff0,*/
  0xab58,0x891b,0xa09d,0xc01c,
  0x1101,0xfe0b,0x3539,0x4025,
  0x70e8,0xa9a7,0xee0b,0xc014,
  0x95ff,0xa6da,0xa2fc,0x3fee,
  0x86c9,0x8f8e,0x33d0,0xbfac,
};
#endif
#ifdef MIEEE
static unsigned short BN16[20] = {
  0xbfd0,0x3518,0xe211,0x6751,
  0x3fe2,0x68bc,0x7072,0x2383,
  0xbfd5,0x1d32,0x6785,0xcf29,
  0x3fb0,0x7f2a,0xa027,0x78a8,
  0xbf6f,0x5604,0x2dba,0xcd1b,
};
static unsigned short BD16[20] = {
  /*0x3ff0,0x0000,0x0000,0x0000,*/
  0xc01c,0xa09d,0x891b,0xab58,
  0x4025,0x3539,0xfe0b,0x1101,
  0xc014,0xee0b,0xa9a7,0x70e8,
  0x3fee,0xa2fc,0xa6da,0x95ff,
  0xbfac,0x33d0,0x8f8e,0x86c9,
};
#endif

#ifdef UNK
static double BPPN[5] = {
  4.65461162774651610328E-1,
  -1.08992173800493920734E0,
  6.38800117371827987759E-1,
  -1.26844349553102907034E-1,
  7.62487844342109852105E-3,
};
static double BPPD[5] = {
  /* 1.00000000000000000000E0,*/
  -8.70622787633159124240E0,
  1.38993162704553213172E1,
  -7.14116144616431159572E0,
  1.34008595960680518666E0,
  -7.84273211323341930448E-2,
};
#endif
#ifdef DEC
static unsigned short BPPN[20] = {
  0037756,0050354,0167531,0135731,
  0140213,0101216,0032767,0020375,
  0040043,0104147,0106312,0177632,
  0137401,0161574,0032015,0043714,
  0036371,0155035,0143165,0142262,
};
static unsigned short BPPD[20] = {
  /*0040200,0000000,0000000,0000000,*/
  0141013,0046265,0115005,0161053,
  0041136,0061631,0072445,0156131,
  0140744,0102145,0001127,0065304,
  0040253,0103757,0146453,0102513,
  0137240,0117200,0155402,0113500,
};
#endif
#ifdef IBMPC
static unsigned short BPPN[20] = {
  0x377b,0x9deb,0xca1d,0x3fdd,
  0xe420,0xc6be,0x7051,0xbff1,
  0x5ff3,0xf199,0x710c,0x3fe4,
  0xa8fa,0x8681,0x3c6f,0xbfc0,
  0xb896,0xb8ce,0x3b43,0x3f7f,
};
static unsigned short BPPD[20] = {
  /*0x0000,0x0000,0x0000,0x3ff0,*/
  0xbc45,0xb340,0x6996,0xc021,
  0xbb8b,0x2ea4,0xcc73,0x402b,
  0xed59,0xa04a,0x908c,0xc01c,
  0x70a9,0xf9a5,0x70fd,0x3ff5,
  0x52e8,0x1b60,0x13d0,0xbfb4,
};
#endif
#ifdef MIEEE
static unsigned short BPPN[20] = {
  0x3fdd,0xca1d,0x9deb,0x377b,
  0xbff1,0x7051,0xc6be,0xe420,
  0x3fe4,0x710c,0xf199,0x5ff3,
  0xbfc0,0x3c6f,0x8681,0xa8fa,
  0x3f7f,0x3b43,0xb8ce,0xb896,
};
static unsigned short BPPD[20] = {
  /*0x3ff0,0x0000,0x0000,0x0000,*/
  0xc021,0x6996,0xb340,0xbc45,
  0x402b,0xcc73,0x2ea4,0xbb8b,
  0xc01c,0x908c,0xa04a,0xed59,
  0x3ff5,0x70fd,0xf9a5,0x70a9,
  0xbfb4,0x13d0,0x1b60,0x52e8,
};
#endif

#ifdef UNK
static double AFN[9] = {
  -1.31696323418331795333E-1,
  -6.26456544431912369773E-1,
  -6.93158036036933542233E-1,
  -2.79779981545119124951E-1,
  -4.91900132609500318020E-2,
  -4.06265923594885404393E-3,
  -1.59276496239262096340E-4,
  -2.77649108155232920844E-6,
  -1.67787698489114633780E-8,
};
static double AFD[9] = {
  /* 1.00000000000000000000E0,*/
  1.33560420706553243746E1,
  3.26825032795224613948E1,
  2.67367040941499554804E1,
  9.18707402907259625840E0,
  1.47529146771666414581E0,
  1.15687173795188044134E-1,
  4.40291641615211203805E-3,
  7.54720348287414296618E-5,
  4.51850092970580378464E-7,
};
#endif
#ifdef DEC
static unsigned short AFN[36] = {
  0137406,0155546,0124127,0033732,
  0140040,0057564,0141263,0041222,
  0140061,0071316,0013674,0175754,
  0137617,0037522,0056637,0120130,
  0137111,0075567,0121755,0166122,
  0136205,0020016,0043317,0002201,
  0135047,0001565,0075130,0002334,
  0133472,0051700,0165021,0131551,
  0131620,0020347,0132165,0013215,
};
static unsigned short AFD[36] = {
  /*0040200,0000000,0000000,0000000,*/
  0041125,0131131,0025627,0067623,
  0041402,0135342,0021703,0154315,
  0041325,0162305,0016671,0120175,
  0041022,0177101,0053114,0141632,
  0040274,0153131,0147364,0114306,
  0037354,0166545,0120042,0150530,
  0036220,0043127,0000727,0130273,
  0034636,0043275,0075667,0034733,
  0032762,0112715,0146250,0142474,
};
#endif
#ifdef IBMPC
static unsigned short AFN[36] = {
  0xe6fb,0xd50a,0xdb6c,0xbfc0,
  0x6852,0x9856,0x0bee,0xbfe4,
  0x9f7d,0xc2f7,0x2e59,0xbfe6,
  0xf40b,0x4bb3,0xe7ea,0xbfd1,
  0xbd8a,0xf47d,0x2f6e,0xbfa9,
  0xe090,0xc8d9,0xa401,0xbf70,
  0x009c,0xaf4b,0xe06e,0xbf24,
  0x366d,0x1d42,0x4a78,0xbec7,
  0xa2d2,0xf68e,0x041c,0xbe52,
};
static unsigned short AFD[36] = {
  /*0x0000,0x0000,0x0000,0x3ff0,*/
  0xedf2,0x2572,0xb64b,0x402a,
  0x7b1a,0x4478,0x575c,0x4040,
  0x3410,0xa3b7,0xbc98,0x403a,
  0x9873,0x2ac9,0x5fc8,0x4022,
  0x9319,0x39de,0x9acb,0x3ff7,
  0x5a2b,0xb404,0x9dac,0x3fbd,
  0xf617,0xe03a,0x08ca,0x3f72,
  0xe73b,0xaf76,0xc8d7,0x3f13,
  0x18a7,0xb995,0x52b9,0x3e9e,
};
#endif
#ifdef MIEEE
static unsigned short AFN[36] = {
  0xbfc0,0xdb6c,0xd50a,0xe6fb,
  0xbfe4,0x0bee,0x9856,0x6852,
  0xbfe6,0x2e59,0xc2f7,0x9f7d,
  0xbfd1,0xe7ea,0x4bb3,0xf40b,
  0xbfa9,0x2f6e,0xf47d,0xbd8a,
  0xbf70,0xa401,0xc8d9,0xe090,
  0xbf24,0xe06e,0xaf4b,0x009c,
  0xbec7,0x4a78,0x1d42,0x366d,
  0xbe52,0x041c,0xf68e,0xa2d2,
};
static unsigned short AFD[36] = {
  /*0x3ff0,0x0000,0x0000,0x0000,*/
  0x402a,0xb64b,0x2572,0xedf2,
  0x4040,0x575c,0x4478,0x7b1a,
  0x403a,0xbc98,0xa3b7,0x3410,
  0x4022,0x5fc8,0x2ac9,0x9873,
  0x3ff7,0x9acb,0x39de,0x9319,
  0x3fbd,0x9dac,0xb404,0x5a2b,
  0x3f72,0x08ca,0xe03a,0xf617,
  0x3f13,0xc8d7,0xaf76,0xe73b,
  0x3e9e,0x52b9,0xb995,0x18a7,
};
#endif

#ifdef UNK
static double AGN[11] = {
  1.97339932091685679179E-2,
  3.91103029615688277255E-1,
  1.06579897599595591108E0,
  9.39169229816650230044E-1,
  3.51465656105547619242E-1,
  6.33888919628925490927E-2,
  5.85804113048388458567E-3,
  2.82851600836737019778E-4,
  6.98793669997260967291E-6,
  8.11789239554389293311E-8,
  3.41551784765923618484E-10,
};
static double AGD[10] = {
  /*  1.00000000000000000000E0,*/
  9.30892908077441974853E0,
  1.98352928718312140417E1,
  1.55646628932864612953E1,
  5.47686069422975497931E0,
  9.54293611618961883998E-1,
  8.64580826352392193095E-2,
  4.12656523824222607191E-3,
  1.01259085116509135510E-4,
  1.17166733214413521882E-6,
  4.91834570062930015649E-9,
};
#endif
#ifdef DEC
static unsigned short AGN[44] = {
  0036641,0124456,0167175,0157354,
  0037710,0037250,0001441,0136671,
  0040210,0066031,0150401,0123532,
  0040160,0066545,0003570,0153133,
  0037663,0171516,0072507,0170345,
  0037201,0151011,0007510,0045702,
  0036277,0172317,0104572,0101030,
  0035224,0045663,0000160,0136422,
  0033752,0074753,0047702,0135160,
  0032256,0052225,0156550,0107103,
  0030273,0142443,0166277,0071720,
};
static unsigned short AGD[40] = {
  /*0040200,0000000,0000000,0000000,*/
  0041024,0170537,0117253,0055003,
  0041236,0127256,0003570,0143240,
  0041171,0004333,0172476,0160645,
  0040657,0041161,0055716,0157161,
  0040164,0046226,0006257,0063431,
  0037261,0010357,0065445,0047563,
  0036207,0034043,0057434,0116732,
  0034724,0055416,0130035,0026377,
  0033235,0041056,0154071,0023502,
  0031250,0177071,0167254,0047242,
};
#endif
#ifdef IBMPC
static unsigned short AGN[44] = {
  0xbbde,0xddcf,0x3525,0x3f94,
  0x37b7,0x0064,0x07d5,0x3fd9,
  0x34eb,0x3a20,0x0d83,0x3ff1,
  0x1acb,0xa0ef,0x0dac,0x3fee,
  0xfe1d,0xcea8,0x7e69,0x3fd6,
  0x0978,0x21e9,0x3a41,0x3fb0,
  0x5043,0xf12f,0xfe99,0x3f77,
  0x17a2,0x600e,0x8976,0x3f32,
  0x574e,0x69f8,0x4f3d,0x3edd,
  0x11c8,0xbbad,0xca92,0x3e75,
  0xee7a,0x7d97,0x78a4,0x3df7,
};
static unsigned short AGD[40] = {
  /*0x0000,0x0000,0x0000,0x3ff0,*/
  0x6b40,0xf3d5,0x9e2b,0x4022,
  0x18d4,0xc0ef,0xd5d5,0x4033,
  0xdc35,0x7ea7,0x211b,0x402f,
  0xdbce,0x2b79,0xe84e,0x4015,
  0xece3,0xc195,0x8992,0x3fee,
  0xa9ee,0xed64,0x221d,0x3fb6,
  0x93bb,0x6be3,0xe704,0x3f70,
  0xa5a0,0xd603,0x8b61,0x3f1a,
  0x24e8,0xdb07,0xa845,0x3eb3,
  0x89d4,0x3dd5,0x1fc7,0x3e35,
};
#endif
#ifdef MIEEE
static unsigned short AGN[44] = {
  0x3f94,0x3525,0xddcf,0xbbde,
  0x3fd9,0x07d5,0x0064,0x37b7,
  0x3ff1,0x0d83,0x3a20,0x34eb,
  0x3fee,0x0dac,0xa0ef,0x1acb,
  0x3fd6,0x7e69,0xcea8,0xfe1d,
  0x3fb0,0x3a41,0x21e9,0x0978,
  0x3f77,0xfe99,0xf12f,0x5043,
  0x3f32,0x8976,0x600e,0x17a2,
  0x3edd,0x4f3d,0x69f8,0x574e,
  0x3e75,0xca92,0xbbad,0x11c8,
  0x3df7,0x78a4,0x7d97,0xee7a,
};
static unsigned short AGD[40] = {
  /*0x3ff0,0x0000,0x0000,0x0000,*/
  0x4022,0x9e2b,0xf3d5,0x6b40,
  0x4033,0xd5d5,0xc0ef,0x18d4,
  0x402f,0x211b,0x7ea7,0xdc35,
  0x4015,0xe84e,0x2b79,0xdbce,
  0x3fee,0x8992,0xc195,0xece3,
  0x3fb6,0x221d,0xed64,0xa9ee,
  0x3f70,0xe704,0x6be3,0x93bb,
  0x3f1a,0x8b61,0xd603,0xa5a0,
  0x3eb3,0xa845,0xdb07,0x24e8,
  0x3e35,0x1fc7,0x3dd5,0x89d4,
};
#endif

#ifdef UNK
static double APFN[9] = {
  1.85365624022535566142E-1,
  8.86712188052584095637E-1,
  9.87391981747398547272E-1,
  4.01241082318003734092E-1,
  7.10304926289631174579E-2,
  5.90618657995661810071E-3,
  2.33051409401776799569E-4,
  4.08718778289035454598E-6,
  2.48379932900442457853E-8,
};
static double APFD[9] = {
  /*  1.00000000000000000000E0,*/
  1.47345854687502542552E1,
  3.75423933435489594466E1,
  3.14657751203046424330E1,
  1.09969125207298778536E1,
  1.78885054766999417817E0,
  1.41733275753662636873E-1,
  5.44066067017226003627E-3,
  9.39421290654511171663E-5,
  5.65978713036027009243E-7,
};
#endif
#ifdef DEC
static unsigned short APFN[36] = {
  0037475,0150174,0071752,0166651,
  0040142,0177621,0164246,0101757,
  0040174,0142670,0106760,0006573,
  0037715,0067570,0116274,0022404,
  0037221,0074157,0053341,0117207,
  0036301,0104257,0015075,0004777,
  0035164,0057502,0164034,0001313,
  0033611,0022254,0176000,0112565,
  0031725,0055523,0025153,0166057,
};
static unsigned short APFD[36] = {
  /*0040200,0000000,0000000,0000000,*/
  0041153,0140334,0130506,0061402,
  0041426,0025551,0024440,0070611,
  0041373,0134750,0047147,0176702,
  0041057,0171532,0105430,0017674,
  0040344,0174416,0001726,0047754,
  0037421,0021207,0020167,0136264,
  0036262,0043621,0151321,0124324,
  0034705,0001313,0163733,0016407,
  0033027,0166702,0150440,0170561,
};
#endif
#ifdef IBMPC
static unsigned short APFN[36] = {
  0x5db5,0x8e7d,0xba0f,0x3fc7,
  0xd07e,0x3d14,0x5ff2,0x3fec,
  0x01af,0x11be,0x98b7,0x3fef,
  0x84a1,0x1397,0xadef,0x3fd9,
  0x33d1,0xeadc,0x2f0d,0x3fb2,
  0xa140,0xe347,0x3115,0x3f78,
  0x8059,0x5d03,0x8be8,0x3f2e,
  0x12af,0x9f80,0x2495,0x3ed1,
  0x7d86,0x654d,0xab6a,0x3e5a,
};
static unsigned short APFD[36] = {
  /*0x0000,0x0000,0x0000,0x3ff0,*/
  0xcc60,0x9628,0x781b,0x402d,
  0x0e31,0x2524,0xc56d,0x4042,
  0xffb8,0x09cc,0x773d,0x403f,
  0x03f7,0x5163,0xfe6b,0x4025,
  0xc9fd,0xc07a,0x9f21,0x3ffc,
  0xf796,0xe40e,0x2450,0x3fc2,
  0x351a,0x3a5a,0x48f2,0x3f76,
  0x63a1,0x7cfb,0xa059,0x3f18,
  0x1e2e,0x5a24,0xfdb8,0x3ea2,
};
#endif
#ifdef MIEEE
static unsigned short APFN[36] = {
  0x3fc7,0xba0f,0x8e7d,0x5db5,
  0x3fec,0x5ff2,0x3d14,0xd07e,
  0x3fef,0x98b7,0x11be,0x01af,
  0x3fd9,0xadef,0x1397,0x84a1,
  0x3fb2,0x2f0d,0xeadc,0x33d1,
  0x3f78,0x3115,0xe347,0xa140,
  0x3f2e,0x8be8,0x5d03,0x8059,
  0x3ed1,0x2495,0x9f80,0x12af,
  0x3e5a,0xab6a,0x654d,0x7d86,
};
static unsigned short APFD[36] = {
  /*0x3ff0,0x0000,0x0000,0x0000,*/
  0x402d,0x781b,0x9628,0xcc60,
  0x4042,0xc56d,0x2524,0x0e31,
  0x403f,0x773d,0x09cc,0xffb8,
  0x4025,0xfe6b,0x5163,0x03f7,
  0x3ffc,0x9f21,0xc07a,0xc9fd,
  0x3fc2,0x2450,0xe40e,0xf796,
  0x3f76,0x48f2,0x3a5a,0x351a,
  0x3f18,0xa059,0x7cfb,0x63a1,
  0x3ea2,0xfdb8,0x5a24,0x1e2e,
};
#endif

#ifdef UNK
static double APGN[11] = {
  -3.55615429033082288335E-2,
  -6.37311518129435504426E-1,
  -1.70856738884312371053E0,
  -1.50221872117316635393E0,
  -5.63606665822102676611E-1,
  -1.02101031120216891789E-1,
  -9.48396695961445269093E-3,
  -4.60325307486780994357E-4,
  -1.14300836484517375919E-5,
  -1.33415518685547420648E-7,
  -5.63803833958893494476E-10,
};
static double APGD[11] = {
  /*  1.00000000000000000000E0,*/
  9.85865801696130355144E0,
  2.16401867356585941885E1,
  1.73130776389749389525E1,
  6.17872175280828766327E0,
  1.08848694396321495475E0,
  9.95005543440888479402E-2,
  4.78468199683886610842E-3,
  1.18159633322838625562E-4,
  1.37480673554219441465E-6,
  5.79912514929147598821E-9,
};
#endif
#ifdef DEC
static unsigned short APGN[44] = {
  0137021,0124372,0176075,0075331,
  0140043,0023330,0177672,0161655,
  0140332,0131126,0010413,0171112,
  0140300,0044263,0175560,0054070,
  0140020,0044206,0142603,0073324,
  0137321,0015130,0066144,0144033,
  0136433,0061243,0175542,0103373,
  0135361,0053721,0020441,0053203,
  0134077,0141725,0160277,0130612,
  0132417,0040372,0100363,0060200,
  0130432,0175052,0171064,0034147,
};
static unsigned short APGD[40] = {
  /*0040200,0000000,0000000,0000000,*/
  0041035,0136420,0030124,0140220,
  0041255,0017432,0034447,0162256,
  0041212,0100456,0154544,0006321,
  0040705,0134026,0127154,0123414,
  0040213,0051612,0044470,0172607,
  0037313,0143362,0053273,0157051,
  0036234,0144322,0054536,0007264,
  0034767,0146170,0054265,0170342,
  0033270,0102777,0167362,0073631,
  0031307,0040644,0167103,0021763,
};
#endif
#ifdef IBMPC
static unsigned short APGN[44] = {
  0xaf5b,0x5f87,0x351f,0xbfa2,
  0x5c76,0x1ff7,0x64db,0xbfe4,
  0x7e49,0xc221,0x564a,0xbffb,
  0x0b07,0x7f6e,0x0916,0xbff8,
  0x6edb,0xd8b0,0x0910,0xbfe2,
  0x9903,0x0d8c,0x234b,0xbfba,
  0x50df,0x7f6c,0x6c54,0xbf83,
  0x2ad0,0x2424,0x2afa,0xbf3e,
  0xf631,0xbc17,0xf87a,0xbee7,
  0x6c10,0x501e,0xe81f,0xbe81,
  0x870d,0x5e46,0x5f45,0xbe03,
};
static unsigned short APGD[40] = {
  /*0x0000,0x0000,0x0000,0x3ff0,*/
  0x9812,0x060a,0xb7a2,0x4023,
  0xfc96,0x4724,0xa3e3,0x4035,
  0x819a,0xdb2c,0x5025,0x4031,
  0x94e2,0xd5cd,0xb702,0x4018,
  0x1eb1,0x4927,0x6a71,0x3ff1,
  0x7bc5,0x4ad7,0x78de,0x3fb9,
  0xc1d7,0x4b2b,0x991a,0x3f73,
  0xbe1c,0x0b16,0xf98f,0x3f1e,
  0x4ef3,0xfdde,0x10bf,0x3eb7,
  0x647e,0x9dc8,0xe834,0x3e38,
};
#endif
#ifdef MIEEE
static unsigned short APGN[44] = {
  0xbfa2,0x351f,0x5f87,0xaf5b,
  0xbfe4,0x64db,0x1ff7,0x5c76,
  0xbffb,0x564a,0xc221,0x7e49,
  0xbff8,0x0916,0x7f6e,0x0b07,
  0xbfe2,0x0910,0xd8b0,0x6edb,
  0xbfba,0x234b,0x0d8c,0x9903,
  0xbf83,0x6c54,0x7f6c,0x50df,
  0xbf3e,0x2afa,0x2424,0x2ad0,
  0xbee7,0xf87a,0xbc17,0xf631,
  0xbe81,0xe81f,0x501e,0x6c10,
  0xbe03,0x5f45,0x5e46,0x870d,
};
static unsigned short APGD[40] = {
  /*0x3ff0,0x0000,0x0000,0x0000,*/
  0x4023,0xb7a2,0x060a,0x9812,
  0x4035,0xa3e3,0x4724,0xfc96,
  0x4031,0x5025,0xdb2c,0x819a,
  0x4018,0xb702,0xd5cd,0x94e2,
  0x3ff1,0x6a71,0x4927,0x1eb1,
  0x3fb9,0x78de,0x4ad7,0x7bc5,
  0x3f73,0x991a,0x4b2b,0xc1d7,
  0x3f1e,0xf98f,0x0b16,0xbe1c,
  0x3eb7,0x10bf,0xfdde,0x4ef3,
  0x3e38,0xe834,0x9dc8,0x647e,
};
#endif

#ifdef ANSIPROT
extern double fabs ( double );
extern double exp ( double );
extern double sqrt ( double );
extern double sin ( double );
extern double cos ( double );
#else
double fabs(), exp(), sqrt();
double polevl(), p1evl(), sin(), cos();
#endif

int airy(double x, double *ai, double *aip, double *bi, double *bip)
{
  double z, zz, t, f, g, uf, ug, k, zeta, theta;
  int domflg;

  domflg = 0;
  if( x > MAXAIRY )
  {
    *ai = 0;
    *aip = 0;
    *bi = MAXNUM;
    *bip = MAXNUM;
    return(-1);
  }

  if( x < -2.09 )
  {
    domflg = 15;
    t = std::sqrt(-x);
    zeta = -2.0 * x * t / 3.0;
    t = std::sqrt(t);
    k = sqpii / t;
    z = 1.0/zeta;
    zz = z * z;
    uf = 1.0 + zz * polevl( zz, AFN, 8 ) / p1evl( zz, AFD, 9 );
    ug = z * polevl( zz, AGN, 10 ) / p1evl( zz, AGD, 10 );
    theta = zeta + 0.25 * PI;
    f = std::sin( theta );
    g = std::cos( theta );
    *ai = k * (f * uf - g * ug);
    *bi = k * (g * uf + f * ug);
    uf = 1.0 + zz * polevl( zz, APFN, 8 ) / p1evl( zz, APFD, 9 );
    ug = z * polevl( zz, APGN, 10 ) / p1evl( zz, APGD, 10 );
    k = sqpii * t;
    *aip = -k * (g * uf + f * ug);
    *bip = k * (f * uf - g * ug);
    return(0);
  }

  if( x >= 2.09 ) /* cbrt(9) */
  {
    domflg = 5;
    t = std::sqrt(x);
    zeta = 2.0 * x * t / 3.0;
    g = exp( zeta );
    t = std::sqrt(t);
    k = 2.0 * t * g;
    z = 1.0/zeta;
    f = polevl( z, AN, 7 ) / polevl( z, AD, 7 );
    *ai = sqpii * f / k;
    k = -0.5 * sqpii * t / g;
    f = polevl( z, APN, 7 ) / polevl( z, APD, 7 );
    *aip = f * k;

    if( x > 8.3203353 ) /* zeta > 16 */
    {
      f = z * polevl( z, BN16, 4 ) / p1evl( z, BD16, 5 );
      k = sqpii * g;
      *bi = k * (1.0 + f) / t;
      f = z * polevl( z, BPPN, 4 ) / p1evl( z, BPPD, 5 );
      *bip = k * t * (1.0 + f);
      return(0);
    }
  }

  f = 1.0;
  g = x;
  t = 1.0;
  uf = 1.0;
  ug = x;
  k = 1.0;
  z = x * x * x;
  while( t > MACHEP )
  {
    uf *= z;
    k += 1.0;
    uf /=k;
    ug *= z;
    k += 1.0;
    ug /=k;
    uf /=k;
    f += uf;
    k += 1.0;
    ug /=k;
    g += ug;
    t = fabs(uf/f);
  }
  uf = c1 * f;
  ug = c2 * g;
  if( (domflg & 1) == 0 )
    *ai = uf - ug;
  if( (domflg & 2) == 0 )
    *bi = sqrt3 * (uf + ug);

  /* the deriviative of ai */
  k = 4.0;
  uf = x * x/2.0;
  ug = z/3.0;
  f = uf;
  g = 1.0 + ug;
  uf /= 3.0;
  t = 1.0;

  while( t > MACHEP )
  {
    uf *= z;
    ug /=k;
    k += 1.0;
    ug *= z;
    uf /=k;
    f += uf;
    k += 1.0;
    ug /=k;
    uf /=k;
    g += ug;
    k += 1.0;
    t = fabs(ug/g);
  }

  uf = c1 * f;
  ug = c2 * g;
  if( (domflg & 4) == 0 )
    *aip = uf - ug;
  if( (domflg & 8) == 0 )
    *bip = sqrt3 * (uf + ug);
  return(0);
}



//////////////////////////// end of    airy.c ///////////////////////





// by mothes:
#define isfinite(x) 1

/*              gamma.c
*
*  Gamma function
*
*
*
* SYNOPSIS:
*
* double x, y, gamma();
* extern int sgngam;
*
* y = gamma( x );
*
*
*
* DESCRIPTION:
*
* Returns gamma function of the argument.  The result is
* correctly signed, and the sign (+1 or -1) is also
* returned in a global (extern) variable named sgngam.
* This variable is also filled in by the std::logarithmic gamma
* function lgam().
*
* Arguments |x| <= 34 are reduced by recurrence and the function
* approximated by a rational function of degree 6/7 in the
* interval (2,3).  Large arguments are handled by Stirling's
* formula. Large negative arguments are made positive ustd::sing
* a reflection formula.
*
*
* ACCURACY:
*
*                      Relative error:
* arithmetic   domain     # trials      peak         rms
*    DEC      -34, 34      10000       1.3e-16     2.5e-17
*    IEEE    -170,-33      20000       2.3e-15     3.3e-16
*    IEEE     -33,  33     20000       9.4e-16     2.2e-16
*    IEEE      33, 171.6   20000       2.3e-15     3.2e-16
*
* Error for arguments outside the test range will be larger
* owing to error amplification by the exponential function.
*
*/
/*              lgam()
*
*  Natural std::logarithm of gamma function
*
*
*
* SYNOPSIS:
*
* double x, y, lgam();
* extern int sgngam;
*
* y = lgam( x );
*
*
*
* DESCRIPTION:
*
* Returns the base e (2.718...) std::logarithm of the absolute
* value of the gamma function of the argument.
* The sign (+1 or -1) of the gamma function is returned in a
* global (extern) variable named sgngam.
*
* For arguments greater than 13, the std::logarithm of the gamma
* function is approximated by the std::logarithmic version of
* Stirling's formula ustd::sing a polynomial approximation of
* degree 4. Arguments between -33 and +33 are reduced by
* recurrence to the interval [2,3] of a rational approximation.
* The std::cosecant reflection formula is employed for arguments
* less than -33.
*
* Arguments greater than MAXLGM return MAXNUM and an error
* message.  MAXLGM = 2.035093e36 for DEC
* arithmetic or 2.556348e305 for IEEE arithmetic.
*
*
*
* ACCURACY:
*
*
* arithmetic      domain        # trials     peak         rms
*    DEC     0, 3                  7000     5.2e-17     1.3e-17
*    DEC     2.718, 2.035e36       5000     3.9e-17     9.9e-18
*    IEEE    0, 3                 28000     5.4e-16     1.1e-16
*    IEEE    2.718, 2.556e305     40000     3.5e-16     8.3e-17
* The error criterion was relative when the function magnitude
* was greater than one but absolute when it was less than one.
*
* The following test used the relative error criterion, though
* at certain points the relative error could be much higher than
* indicated.
*    IEEE    -200, -4             10000     4.8e-16     1.3e-16
*
*/

/*              gamma.c */
/*  gamma function  */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier
*/



#ifdef UNK
static double P[] = {
  1.60119522476751861407E-4,
  1.19135147006586384913E-3,
  1.04213797561761569935E-2,
  4.76367800457137231464E-2,
  2.07448227648435975150E-1,
  4.94214826801497100753E-1,
  9.99999999999999996796E-1
};
static double Q[] = {
  -2.31581873324120129819E-5,
  5.39605580493303397842E-4,
  -4.45641913851797240494E-3,
  1.18139785222060435552E-2,
  3.58236398605498653373E-2,
  -2.34591795718243348568E-1,
  7.14304917030273074085E-2,
  1.00000000000000000320E0
};
#define MAXGAM 171.624376956302725
static double LOGPI = 1.14472988584940017414;
#endif

#ifdef DEC
static unsigned short P[] = {
  0035047,0162701,0146301,0005234,
  0035634,0023437,0032065,0176530,
  0036452,0137157,0047330,0122574,
  0037103,0017310,0143041,0017232,
  0037524,0066516,0162563,0164605,
  0037775,0004671,0146237,0014222,
  0040200,0000000,0000000,0000000
};
static unsigned short Q[] = {
  0134302,0041724,0020006,0116565,
  0035415,0072121,0044251,0025634,
  0136222,0003447,0035205,0121114,
  0036501,0107552,0154335,0104271,
  0037022,0135717,0014776,0171471,
  0137560,0034324,0165024,0037021,
  0037222,0045046,0047151,0161213,
  0040200,0000000,0000000,0000000
};
#define MAXGAM 34.84425627277176174
static unsigned short LPI[4] = {
  0040222,0103202,0043475,0006750,
};
#define std::logPI *(double *)LPI
#endif

#ifdef IBMPC
static unsigned short P[] = {
  0x2153,0x3998,0xfcb8,0x3f24,
  0xbfab,0xe686,0x84e3,0x3f53,
  0x14b0,0xe9db,0x57cd,0x3f85,
  0x23d3,0x18c4,0x63d9,0x3fa8,
  0x7d31,0xdcae,0x8da9,0x3fca,
  0xe312,0x3993,0xa137,0x3fdf,
  0x0000,0x0000,0x0000,0x3ff0
};
static unsigned short Q[] = {
  0xd3af,0x8400,0x487a,0xbef8,
  0x2573,0x2915,0xae8a,0x3f41,
  0xb44a,0xe750,0x40e4,0xbf72,
  0xb117,0x5b1b,0x31ed,0x3f88,
  0xde67,0xe33f,0x5779,0x3fa2,
  0x87c2,0x9d42,0x071a,0xbfce,
  0x3c51,0xc9cd,0x4944,0x3fb2,
  0x0000,0x0000,0x0000,0x3ff0
};
#define MAXGAM 171.624376956302725
static unsigned short LPI[4] = {
  0xa1bd,0x48e7,0x50d0,0x3ff2,
};
#define std::logPI *(double *)LPI
#endif

#ifdef MIEEE
static unsigned short P[] = {
  0x3f24,0xfcb8,0x3998,0x2153,
  0x3f53,0x84e3,0xe686,0xbfab,
  0x3f85,0x57cd,0xe9db,0x14b0,
  0x3fa8,0x63d9,0x18c4,0x23d3,
  0x3fca,0x8da9,0xdcae,0x7d31,
  0x3fdf,0xa137,0x3993,0xe312,
  0x3ff0,0x0000,0x0000,0x0000
};
static unsigned short Q[] = {
  0xbef8,0x487a,0x8400,0xd3af,
  0x3f41,0xae8a,0x2915,0x2573,
  0xbf72,0x40e4,0xe750,0xb44a,
  0x3f88,0x31ed,0x5b1b,0xb117,
  0x3fa2,0x5779,0xe33f,0xde67,
  0xbfce,0x071a,0x9d42,0x87c2,
  0x3fb2,0x4944,0xc9cd,0x3c51,
  0x3ff0,0x0000,0x0000,0x0000
};
#define MAXGAM 171.624376956302725
static unsigned short LPI[4] = {
  0x3ff2,0x50d0,0x48e7,0xa1bd,
};
#define std::logPI *(double *)LPI
#endif

/* Stirling's formula for the gamma function */
#if UNK
static double STIR[5] = {
  7.87311395793093628397E-4,
  -2.29549961613378126380E-4,
  -2.68132617805781232825E-3,
  3.47222221605458667310E-3,
  8.33333333333482257126E-2,
};
#define MAXSTIR 143.01608
static double SQTPI = 2.50662827463100050242E0;
#endif
#if DEC
static unsigned short STIR[20] = {
  0035516,0061622,0144553,0112224,
  0135160,0131531,0037460,0165740,
  0136057,0134460,0037242,0077270,
  0036143,0107070,0156306,0027751,
  0037252,0125252,0125252,0146064,
};
#define MAXSTIR 26.77
static unsigned short SQT[4] = {
  0040440,0066230,0177661,0034055,
};
#define SQTPI *(double *)SQT
#endif
#if IBMPC
static unsigned short STIR[20] = {
  0x7293,0x592d,0xcc72,0x3f49,
  0x1d7c,0x27e6,0x166b,0xbf2e,
  0x4fd7,0x07d4,0xf726,0xbf65,
  0xc5fd,0x1b98,0x71c7,0x3f6c,
  0x5986,0x5555,0x5555,0x3fb5,
};
#define MAXSTIR 143.01608
static unsigned short SQT[4] = {
  0x2706,0x1ff6,0x0d93,0x4004,
};
#define SQTPI *(double *)SQT
#endif
#if MIEEE
static unsigned short STIR[20] = {
  0x3f49,0xcc72,0x592d,0x7293,
  0xbf2e,0x166b,0x27e6,0x1d7c,
  0xbf65,0xf726,0x07d4,0x4fd7,
  0x3f6c,0x71c7,0x1b98,0xc5fd,
  0x3fb5,0x5555,0x5555,0x5986,
};
#define MAXSTIR 143.01608
static unsigned short SQT[4] = {
  0x4004,0x0d93,0x1ff6,0x2706,
};
#define SQTPI *(double *)SQT
#endif

int sgngam = 0;
extern int sgngam;
extern double MAXLOG, MAXNUM, PI;
#ifdef ANSIPROT
extern double pow ( double, double );
extern double log ( double );
extern double exp ( double );
extern double sin ( double );
// extern double polevl ( double, void *, int );
// extern double p1evl ( double, void *, int );
// extern double floor ( double );
extern double fabs ( double );
extern int isnan ( double );
// extern int isfinite ( double );
static double stirf ( double );
double lgam ( double );
#else
double pow(), log(), exp(), sin(), polevl(), p1evl(), floor(), fabs();
int isnan(), isfinite();
static double stirf();
double lgam();
#endif
#ifdef INFINITIES
extern double INFINITY_BESSEL;
#endif
#ifdef NANS
extern double NAN_BESSEL;
#endif

/* Gamma function computed by Stirling's formula.
* The polynomial STIR is valid for 33 <= x <= 172.
*/
static double stirf(double x)
{
  double y, w, v;

  w = 1.0/x;
  w = 1.0 + w * polevl( w, STIR, 4 );
  y = exp(x);
  if( x > MAXSTIR )
  { /* Avoid overflow in std::pow() */
    v = std::pow( x, 0.5 * x - 0.25 );
    y = v * (v / y);
  }
  else
  {
    y = std::pow( x, x - 0.5 ) / y;
  }
  y = SQTPI * y * w;
  return( y );
}



double gamma(double x)
{
  double p, q, z;
  int i;

  sgngam = 1;
#ifdef NANS
  if( isnan(x) )
    return(x);
#endif
#ifdef INFINITIES
#ifdef NANS
  if( x == INFINITY_BESSEL )
    return(x);
  if( x == -INFINITY_BESSEL )
    return(NAN_BESSEL);
#else
  if( !isfinite(x) )
    return(x);
#endif
#endif
  q = fabs(x);

  if( q > 33.0 )
  {
    if( x < 0.0 )
    {
      p = floor(q);
      if( p == q )
      {
#ifdef NANS
gamnan:
        mtherr( "gamma", DOMAIN );
        return (NAN_BESSEL);
#else
        goto goverf;
#endif
      }
      i = (int) p;
      if( (i & 1) == 0 )
        sgngam = -1;
      z = q - p;
      if( z > 0.5 )
      {
        p += 1.0;
        z = q - p;
      }
      z = q * std::sin( PI * z );
      if( z == 0.0 )
      {
#ifdef INFINITIES
        return( sgngam * INFINITY_BESSEL);
#else
goverf:
        mtherr( "gamma", OVERFLOW );
        return( sgngam * MAXNUM);
#endif
      }
      z = fabs(z);
      z = PI/(z * stirf(q) );
    }
    else
    {
      z = stirf(x);
    }
    return( sgngam * z );
  }

  z = 1.0;
  while( x >= 3.0 )
  {
    x -= 1.0;
    z *= x;
  }

  while( x < 0.0 )
  {
    if( x > -1.E-9 )
      goto small_;
    z /= x;
    x += 1.0;
  }

  while( x < 2.0 )
  {
    if( x < 1.e-9 )
      goto small_;
    z /= x;
    x += 1.0;
  }

  if( x == 2.0 )
    return(z);

  x -= 2.0;
  p = polevl( x, P, 6 );
  q = polevl( x, Q, 7 );
  return( z * p / q );

small_:
  if( x == 0.0 )
  {
#ifdef INFINITIES
#ifdef NANS
    goto gamnan;
#else
    return( INFINITY_BESSEL );
#endif
#else
    mtherr( "gamma", SING );
    return( MAXNUM );
#endif
  }
  else
    return( z/((1.0 + 0.5772156649015329 * x) * x) );
}



/* A[]: Stirling's formula expansion of std::log gamma
* B[], C[]: std::log gamma function between 2 and 3
*/
#ifdef UNK
static double A[] = {
  8.11614167470508450300E-4,
  -5.95061904284301438324E-4,
  7.93650340457716943945E-4,
  -2.77777777730099687205E-3,
  8.33333333333331927722E-2
};
static double B[] = {
  -1.37825152569120859100E3,
  -3.88016315134637840924E4,
  -3.31612992738871184744E5,
  -1.16237097492762307383E6,
  -1.72173700820839662146E6,
  -8.53555664245765465627E5
};
static double C[] = {
  /* 1.00000000000000000000E0, */
  -3.51815701436523470549E2,
  -1.70642106651881159223E4,
  -2.20528590553854454839E5,
  -1.13933444367982507207E6,
  -2.53252307177582951285E6,
  -2.01889141433532773231E6
};
/* std::log( std::sqrt( 2*pi ) ) */
static double LS2PI  =  0.91893853320467274178;
#define MAXLGM 2.556348e305
#endif

#ifdef DEC
static unsigned short A[] = {
  0035524,0141201,0034633,0031405,
  0135433,0176755,0126007,0045030,
  0035520,0006371,0003342,0172730,
  0136066,0005540,0132605,0026407,
  0037252,0125252,0125252,0125132
};
static unsigned short B[] = {
  0142654,0044014,0077633,0035410,
  0144027,0110641,0125335,0144760,
  0144641,0165637,0142204,0047447,
  0145215,0162027,0146246,0155211,
  0145322,0026110,0010317,0110130,
  0145120,0061472,0120300,0025363
};
static unsigned short C[] = {
  /*0040200,0000000,0000000,0000000*/
  0142257,0164150,0163630,0112622,
  0143605,0050153,0156116,0135272,
  0144527,0056045,0145642,0062332,
  0145213,0012063,0106250,0001025,
  0145432,0111254,0044577,0115142,
  0145366,0071133,0050217,0005122
};
/* std::log( std::sqrt( 2*pi ) ) */
static unsigned short LS2P[] = {040153,037616,041445,0172645,};
#define LS2PI *(double *)LS2P
#define MAXLGM 2.035093e36
#endif

#ifdef IBMPC
static unsigned short A[] = {
  0x6661,0x2733,0x9850,0x3f4a,
  0xe943,0xb580,0x7fbd,0xbf43,
  0x5ebb,0x20dc,0x019f,0x3f4a,
  0xa5a1,0x16b0,0xc16c,0xbf66,
  0x554b,0x5555,0x5555,0x3fb5
};
static unsigned short B[] = {
  0x6761,0x8ff3,0x8901,0xc095,
  0xb93e,0x355b,0xf234,0xc0e2,
  0x89e5,0xf890,0x3d73,0xc114,
  0xdb51,0xf994,0xbc82,0xc131,
  0xf20b,0x0219,0x4589,0xc13a,
  0x055e,0x5418,0x0c67,0xc12a
};
static unsigned short C[] = {
  /*0x0000,0x0000,0x0000,0x3ff0,*/
  0x12b2,0x1cf3,0xfd0d,0xc075,
  0xd757,0x7b89,0xaa0d,0xc0d0,
  0x4c9b,0xb974,0xeb84,0xc10a,
  0x0043,0x7195,0x6286,0xc131,
  0xf34c,0x892f,0x5255,0xc143,
  0xe14a,0x6a11,0xce4b,0xc13e
};
/* std::log( std::sqrt( 2*pi ) ) */
static unsigned short LS2P[] = {
  0xbeb5,0xc864,0x67f1,0x3fed
};
#define LS2PI *(double *)LS2P
#define MAXLGM 2.556348e305
#endif

#ifdef MIEEE
static unsigned short A[] = {
  0x3f4a,0x9850,0x2733,0x6661,
  0xbf43,0x7fbd,0xb580,0xe943,
  0x3f4a,0x019f,0x20dc,0x5ebb,
  0xbf66,0xc16c,0x16b0,0xa5a1,
  0x3fb5,0x5555,0x5555,0x554b
};
static unsigned short B[] = {
  0xc095,0x8901,0x8ff3,0x6761,
  0xc0e2,0xf234,0x355b,0xb93e,
  0xc114,0x3d73,0xf890,0x89e5,
  0xc131,0xbc82,0xf994,0xdb51,
  0xc13a,0x4589,0x0219,0xf20b,
  0xc12a,0x0c67,0x5418,0x055e
};
static unsigned short C[] = {
  0xc075,0xfd0d,0x1cf3,0x12b2,
  0xc0d0,0xaa0d,0x7b89,0xd757,
  0xc10a,0xeb84,0xb974,0x4c9b,
  0xc131,0x6286,0x7195,0x0043,
  0xc143,0x5255,0x892f,0xf34c,
  0xc13e,0xce4b,0x6a11,0xe14a
};
/* std::log( std::sqrt( 2*pi ) ) */
static unsigned short LS2P[] = {
  0x3fed,0x67f1,0xc864,0xbeb5
};
#define LS2PI *(double *)LS2P
#define MAXLGM 2.556348e305
#endif


/* std::logarithm of gamma function */


double lgam(double x)
{
  double p, q, u, w, z;
  int i;

  sgngam = 1;
#ifdef NANS
  if( isnan(x) )
    return(x);
#endif

#ifdef INFINITIES
  if( !isfinite(x) )
    return(INFINITY_BESSEL);
#endif

  if( x < -34.0 )
  {
    q = -x;
    w = lgam(q); /* note this modifies sgngam! */
    p = floor(q);
    if( p == q )
    {
lgsing:
#ifdef INFINITIES
      mtherr( "lgam", SING );
      return (INFINITY_BESSEL);
#else
      goto loverf;
#endif
    }
    i = (int) p;
    if( (i & 1) == 0 )
      sgngam = -1;
    else
      sgngam = 1;
    z = q - p;
    if( z > 0.5 )
    {
      p += 1.0;
      z = p - q;
    }
    z = q * std::sin( PI * z );
    if( z == 0.0 )
      goto lgsing;
    /*  z = std::log(PI) - std::log( z ) - w;*/
    z = LOGPI - std::log( z ) - w;
    return( z );
  }

  if( x < 13.0 )
  {
    z = 1.0;
    p = 0.0;
    u = x;
    while( u >= 3.0 )
    {
      p -= 1.0;
      u = x + p;
      z *= u;
    }
    while( u < 2.0 )
    {
      if( u == 0.0 )
        goto lgsing;
      z /= u;
      p += 1.0;
      u = x + p;
    }
    if( z < 0.0 )
    {
      sgngam = -1;
      z = -z;
    }
    else
      sgngam = 1;
    if( u == 2.0 )
      return( std::log(z) );
    p -= 2.0;
    x = x + p;
    p = x * polevl( x, B, 5 ) / p1evl( x, C, 6);
    return( std::log(z) + p );
  }

  if( x > MAXLGM )
  {
#ifdef INFINITIES
    return( sgngam * INFINITY_BESSEL );
#else
loverf:
    mtherr( "lgam", OVERFLOW );
    return( sgngam * MAXNUM );
#endif
  }

  q = ( x - 0.5 ) * std::log(x) - x + LS2PI;
  if( x > 1.0e8 )
    return( q );

  p = 1.0/(x*x);
  if( x >= 1000.0 )
    q += ((   7.9365079365079365079365e-4 * p
    - 2.7777777777777777777778e-3) *p
    + 0.0833333333333333333333) / x;
  else
    q += polevl( p, A, 4 ) / x;
  return( q );
}

//////////////////////////// end of    gamma.c ///////////////////////


////////////////////////// bessel funcion ////////////////////////////

/*							jv.c
*
*	Bessel function of noninteger order
*
*
*
* SYNOPSIS:
*
* double v, x, y, jv();
*
* y = jv( v, x );
*
*
*
* DESCRIPTION:
*
* Returns Bessel function of order v of the argument,
* where v is real.  Negative x is allowed if v is an integer.
*
* Several expansions are included: the ascending std::power
* series, the Hankel expansion, and two transitional
* expansions for large v.  If v is not too large, it
* is reduced by recurrence to a region of best accuracy.
* The transitional expansions give 12D accuracy for v > 500.
*
*
*
* ACCURACY:
* Results for integer v are indicated by *, where x and v
* both vary from -125 to +125.  Otherwise,
* x ranges from 0 to 125, v ranges as indicated by "domain."
* Error criterion is absolute, except relative when |jv()| > 1.
*
* arithmetic  v domain  x domain    # trials      peak       rms
*    IEEE      0,125     0,125      100000      4.6e-15    2.2e-16
*    IEEE   -125,0       0,125       40000      5.4e-11    3.7e-13
*    IEEE      0,500     0,500       20000      4.4e-15    4.0e-16
* Integer v:
*    IEEE   -125,125   -125,125      50000      3.5e-15*   1.9e-16*
*
*/


/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier
*/


#define DEBUG 0

#ifdef DEC
#define MAXGAM 34.84425627277176174
#else
#define MAXGAM 171.624376956302725
#endif

#ifdef ANSIPROT
// extern int airy ( double, double *, double *, double *, double * );
extern double fabs ( double );
// extern double floor ( double );
// extern double frexp ( double, int * );
// extern double polevl ( double, void *, int );
// extern double j0 ( double );
// extern double j1 ( double );
extern double sqrt ( double );
extern double cbrt ( double );
extern double exp ( double );
extern double log ( double );
extern double sin ( double );
extern double cos ( double );
extern double acos ( double );
extern double pow ( double, double );
extern double lgam ( double );
static double recur(double *, double, double *, int);
static double jvs(double, double);
static double hankel(double, double);
static double jnx(double, double);
static double jnt(double, double);
#else
int airy();
double fabs(), floor(), frexp(), polevl(), j0(), j1(), sqrt(), cbrt();
double exp(), log(), sin(), cos(), acos(), pow(), gamma(), lgam();
static double recur(), jvs(), hankel(), jnx(), jnt();
#endif

extern double MAXNUM, MACHEP, MINLOG, MAXLOG;
#define BIG  1.44115188075855872E+17

double jv(double n, double x)
{
  double k, q, t, y, an;
  int i, sign, nint;

  nint = 0;	/* Flag for integer n */
  sign = 1;	/* Flag for sign inversion */
  an = fabs( n );
  y = floor( an );
  if( y == an )
  {
    nint = 1;
    i = (int) (an - 16384.0 * floor( an/16384.0 ));
    if( n < 0.0 )
    {
      if( i & 1 )
        sign = -sign;
      n = an;
    }
    if( x < 0.0 )
    {
      if( i & 1 )
        sign = -sign;
      x = -x;
    }
    if( n == 0.0 )
      return( j0(x) );
    if( n == 1.0 )
      return( sign * j1(x) );
  }

  if( (x < 0.0) && (y != an) )
  {
    mtherr( "Jv", DOMAIN );
    y = 0.0;
    goto done;
  }

  y = fabs(x);

  if( y < MACHEP )
    goto underf;

  k = 3.6 * std::sqrt(y);
  t = 3.6 * std::sqrt(an);
  if( (y < t) && (an > 21.0) )
    return( sign * jvs(n,x) );
  if( (an < k) && (y > 21.0) )
    return( sign * hankel(n,x) );

  if( an < 500.0 )
  {
    /* Note: if x is too large, the continued
    * fraction will fail; but then the
    * Hankel expansion can be used.
    */
    if( nint != 0 )
    {
      k = 0.0;
      q = recur( &n, x, &k, 1 );
      if( k == 0.0 )
      {
        y = j0(x)/q;
        goto done;
      }
      if( k == 1.0 )
      {
        y = j1(x)/q;
        goto done;
      }
    }

    if( an > 2.0 * y )
      goto rlarger;

    if( (n >= 0.0) && (n < 20.0)
      && (y > 6.0) && (y < 20.0) )
    {
      /* Recur backwards from a larger value of n
      */
rlarger:
      k = n;

      y = y + an + 1.0;
      if( y < 30.0 )
        y = 30.0;
      y = n + floor(y-n);
      q = recur( &y, x, &k, 0 );
      y = jvs(y,x) * q;
      goto done;
    }

    if( k <= 30.0 )
    {
      k = 2.0;
    }
    else if( k < 90.0 )
    {
      k = (3*k)/4;
    }
    if( an > (k + 3.0) )
    {
      if( n < 0.0 )
        k = -k;
      q = n - floor(n);
      k = floor(k) + q;
      if( n > 0.0 )
        q = recur( &n, x, &k, 1 );
      else
      {
        t = k;
        k = n;
        q = recur( &t, x, &k, 1 );
        k = t;
      }
      if( q == 0.0 )
      {
underf:
        y = 0.0;
        goto done;
      }
    }
    else
    {
      k = n;
      q = 1.0;
    }

    /* boundary between convergence of
    * std::power series and Hankel expansion
    */
    y = fabs(k);
    if( y < 26.0 )
      t = (0.0083*y + 0.09)*y + 12.9;
    else
      t = 0.9 * y;

    if( x > t )
      y = hankel(k,x);
    else
      y = jvs(k,x);
#if DEBUG
    printf( "y = %.16e, recur q = %.16e\n", y, q );
#endif
    if( n > 0.0 )
      y /= q;
    else
      y *= q;
  }

  else
  {
    /* For large n, use the uniform expansion
    * or the transitional expansion.
    * But if x is of the order of n**2,
    * these may blow up, whereas the
    * Hankel expansion will then work.
    */
    if( n < 0.0 )
    {
      mtherr( "Jv", TLOSS );
      y = 0.0;
      goto done;
    }
    t = x/n;
    t /= n;
    if( t > 0.3 )
      y = hankel(n,x);
    else
      y = jnx(n,x);
  }

done:	return( sign * y);
}

/* Reduce the order by backward recurrence.
* AMS55 #9.1.27 and 9.1.73.
*/

static double recur(double *n, double x, double *newn, int cancel)
{
  double pkm2, pkm1, pk, qkm2, qkm1;
  /* double pkp1; */
  double k, ans, qk, xk, yk, r, t, kf;
  static double big = BIG;
  int nflag, ctr;

  /* continued fraction for Jn(x)/Jn-1(x)  */
  if( *n < 0.0 )
    nflag = 1;
  else
    nflag = 0;

fstart:

#if DEBUG
  printf( "recur: n = %.6e, newn = %.6e, cfrac = ", *n, *newn );
#endif

  pkm2 = 0.0;
  qkm2 = 1.0;
  pkm1 = x;
  qkm1 = *n + *n;
  xk = -x * x;
  yk = qkm1;
  ans = 1.0;
  ctr = 0;
  do
  {
    yk += 2.0;
    pk = pkm1 * yk +  pkm2 * xk;
    qk = qkm1 * yk +  qkm2 * xk;
    pkm2 = pkm1;
    pkm1 = pk;
    qkm2 = qkm1;
    qkm1 = qk;
    if( qk != 0 )
      r = pk/qk;
    else
      r = 0.0;
    if( r != 0 )
    {
      t = fabs( (ans - r)/r );
      ans = r;
    }
    else
      t = 1.0;

    if( ++ctr > 1000 )
    {
      mtherr( "jv", UNDERFLOW );
      goto done;
    }
    if( t < MACHEP )
      goto done;

    if( fabs(pk) > big )
    {
      pkm2 /= big;
      pkm1 /= big;
      qkm2 /= big;
      qkm1 /= big;
    }
  }
  while( t > MACHEP );

done:

#if DEBUG
  printf( "%.6e\n", ans );
#endif

  /* Change n to n-1 if n < 0 and the continued fraction is small
  */
  if( nflag > 0 )
  {
    if( fabs(ans) < 0.125 )
    {
      nflag = -1;
      *n = *n - 1.0;
      goto fstart;
    }
  }


  kf = *newn;

  /* backward recurrence
  *              2k
  *  J   (x)  =  --- J (x)  -  J   (x)
  *   k-1         x   k         k+1
  */

  pk = 1.0;
  pkm1 = 1.0/ans;
  k = *n - 1.0;
  r = 2 * k;
  do
  {
    pkm2 = (pkm1 * r  -  pk * x) / x;
    /*	pkp1 = pk; */
    pk = pkm1;
    pkm1 = pkm2;
    r -= 2.0;
    /*
    t = fabs(pkp1) + fabs(pk);
    if( (k > (kf + 2.5)) && (fabs(pkm1) < 0.25*t) )
    {
    k -= 1.0;
    t = x*x;
    pkm2 = ( (r*(r+2.0)-t)*pk - r*x*pkp1 )/t;
    pkp1 = pk;
    pk = pkm1;
    pkm1 = pkm2;
    r -= 2.0;
    }
    */
    k -= 1.0;
  }
  while( k > (kf + 0.5) );

  /* Take the larger of the last two iterates
  * on the theory that it may have less cancellation error.
  */

  if( cancel )
  {
    if( (kf >= 0.0) && (fabs(pk) > fabs(pkm1)) )
    {
      k += 1.0;
      pkm2 = pk;
    }
  }
  *newn = k;
#if DEBUG
  printf( "newn %.6e rans %.6e\n", k, pkm2 );
#endif
  return( pkm2 );
}



/* Ascending std::power series for Jv(x).
* AMS55 #9.1.10.
*/

extern double PI;
extern int sgngam;

static double jvs(double n, double x)
{
  double t, u, y, z, k;
  int ex;

  z = -x * x / 4.0;
  u = 1.0;
  y = u;
  k = 1.0;
  t = 1.0;

  while( t > MACHEP )
  {
    u *= z / (k * (n+k));
    y += u;
    k += 1.0;
    if( y != 0 )
      t = fabs( u/y );
  }
#if DEBUG
  printf( "std::power series=%.5e ", y );
#endif
  t = frexp( 0.5*x, &ex );
  ex = (int) (ex * n);
  if(  (ex > -1023)
    && (ex < 1023)
    && (n > 0.0)
    && (n < (MAXGAM-1.0)) )
  {
    t = std::pow( 0.5*x, n ) / gamma( n + 1.0 );
#if DEBUG
    printf( "std::pow(.5*x, %.4e)/gamma(n+1)=%.5e\n", n, t );
#endif
    y *= t;
  }
  else
  {
#if DEBUG
    z = n * std::log(0.5*x);
    k = lgam( n+1.0 );
    t = z - k;
    printf( "std::log std::pow=%.5e, lgam(%.4e)=%.5e\n", z, n+1.0, k );
#else
    t = n * std::log(0.5*x) - lgam(n + 1.0);
#endif
    if( y < 0 )
    {
      sgngam = -sgngam;
      y = -y;
    }
    t += std::log(y);
#if DEBUG
    printf( "std::log y=%.5e\n", std::log(y) );
#endif
    if( t < -MAXLOG )
    {
      return( 0.0 );
    }
    if( t > MAXLOG )
    {
      mtherr( "Jv", OVERFLOW );
      return( MAXNUM );
    }
    y = sgngam * exp( t );
  }
  return(y);
}

/* Hankel's asymptotic expansion
* for large x.
* AMS55 #9.2.5.
*/

static double hankel(double n, double x)
{
  double t, u, z, k, sign, conv;
  double p, q, j, m, pp, qq;
  int flag;

  m = 4.0*n*n;
  j = 1.0;
  z = 8.0 * x;
  k = 1.0;
  p = 1.0;
  u = (m - 1.0)/z;
  q = u;
  sign = 1.0;
  conv = 1.0;
  flag = 0;
  t = 1.0;
  pp = 1.0e38;
  qq = 1.0e38;

  while( t > MACHEP )
  {
    k += 2.0;
    j += 1.0;
    sign = -sign;
    u *= (m - k * k)/(j * z);
    p += sign * u;
    k += 2.0;
    j += 1.0;
    u *= (m - k * k)/(j * z);
    q += sign * u;
    t = fabs(u/p);
    if( t < conv )
    {
      conv = t;
      qq = q;
      pp = p;
      flag = 1;
    }
    /* stop if the terms start getting larger */
    if( (flag != 0) && (t > conv) )
    {
#if DEBUG
      printf( "Hankel: convergence to %.4E\n", conv );
#endif
      goto hank1;
    }
  }

hank1:
  u = x - (0.5*n + 0.25) * PI;
  t = std::sqrt( 2.0/(PI*x) ) * ( pp * std::cos(u) - qq * std::sin(u) );
#if DEBUG
  printf( "hank: %.6e\n", t );
#endif
  return( t );
}


/* Asymptotic expansion for large n.
* AMS55 #9.3.35.
*/

static double lambda[] = {
  1.0,
  1.041666666666666666666667E-1,
  8.355034722222222222222222E-2,
  1.282265745563271604938272E-1,
  2.918490264641404642489712E-1,
  8.816272674437576524187671E-1,
  3.321408281862767544702647E+0,
  1.499576298686255465867237E+1,
  7.892301301158651813848139E+1,
  4.744515388682643231611949E+2,
  3.207490090890661934704328E+3
};
static double mu[] = {
  1.0,
  -1.458333333333333333333333E-1,
  -9.874131944444444444444444E-2,
  -1.433120539158950617283951E-1,
  -3.172272026784135480967078E-1,
  -9.424291479571202491373028E-1,
  -3.511203040826354261542798E+0,
  -1.572726362036804512982712E+1,
  -8.228143909718594444224656E+1,
  -4.923553705236705240352022E+2,
  -3.316218568547972508762102E+3
};
static double P1[] = {
  -2.083333333333333333333333E-1,
  1.250000000000000000000000E-1
};
static double P2[] = {
  3.342013888888888888888889E-1,
  -4.010416666666666666666667E-1,
  7.031250000000000000000000E-2
};
static double P3[] = {
  -1.025812596450617283950617E+0,
  1.846462673611111111111111E+0,
  -8.912109375000000000000000E-1,
  7.324218750000000000000000E-2
};
static double P4[] = {
  4.669584423426247427983539E+0,
  -1.120700261622299382716049E+1,
  8.789123535156250000000000E+0,
  -2.364086914062500000000000E+0,
  1.121520996093750000000000E-1
};
static double P5[] = {
  -2.8212072558200244877E1,
  8.4636217674600734632E1,
  -9.1818241543240017361E1,
  4.2534998745388454861E1,
  -7.3687943594796316964E0,
  2.27108001708984375E-1
};
static double P6[] = {
  2.1257013003921712286E2,
  -7.6525246814118164230E2,
  1.0599904525279998779E3,
  -6.9957962737613254123E2,
  2.1819051174421159048E2,
  -2.6491430486951555525E1,
  5.7250142097473144531E-1
};
static double P7[] = {
  -1.9194576623184069963E3,
  8.0617221817373093845E3,
  -1.3586550006434137439E4,
  1.1655393336864533248E4,
  -5.3056469786134031084E3,
  1.2009029132163524628E3,
  -1.0809091978839465550E2,
  1.7277275025844573975E0
};


static double jnx(double n, double x)
{
  double zeta, sqz, zz, zp, np;
  double cbn, n23, t, z, sz;
  double pp, qq, z32i, zzi;
  double ak, bk, akl, bkl;
  int sign, doa, dob, nflg, k, s, tk, tkp1, m;
  static double u[8];
  static double ai, aip, bi, bip;

  /* Test for x very close to n.
  * Use expansion for transition region if so.
  */
  cbn = cbrt(n);
  z = (x - n)/cbn;
  if( fabs(z) <= 0.7 )
    return( jnt(n,x) );

  z = x/n;
  zz = 1.0 - z*z;
  if( zz == 0.0 )
    return(0.0);

  if( zz > 0.0 )
  {
    sz = std::sqrt( zz );
    t = 1.5 * (std::log( (1.0+sz)/z ) - sz );	/* zeta ** 3/2		*/
    zeta = cbrt( t * t );
    nflg = 1;
  }
  else
  {
    sz = std::sqrt(-zz);
    t = 1.5 * (sz - acos(1.0/z));
    zeta = -cbrt( t * t );
    nflg = -1;
  }
  z32i = fabs(1.0/t);
  sqz = cbrt(t);

  /* Airy function */
  n23 = cbrt( n * n );
  t = n23 * zeta;

#if DEBUG
  printf("zeta %.5E, Airy(%.5E)\n", zeta, t );
#endif
  airy( t, &ai, &aip, &bi, &bip );

  /* polynomials in expansion */
  u[0] = 1.0;
  zzi = 1.0/zz;
  u[1] = polevl( zzi, P1, 1 )/sz;
  u[2] = polevl( zzi, P2, 2 )/zz;
  u[3] = polevl( zzi, P3, 3 )/(sz*zz);
  pp = zz*zz;
  u[4] = polevl( zzi, P4, 4 )/pp;
  u[5] = polevl( zzi, P5, 5 )/(pp*sz);
  pp *= zz;
  u[6] = polevl( zzi, P6, 6 )/pp;
  u[7] = polevl( zzi, P7, 7 )/(pp*sz);

#if DEBUG
  for( k=0; k<=7; k++ )
    printf( "u[%d] = %.5E\n", k, u[k] );
#endif

  pp = 0.0;
  qq = 0.0;
  np = 1.0;
  /* flags to stop when terms get larger */
  doa = 1;
  dob = 1;
  akl = MAXNUM;
  bkl = MAXNUM;

  for( k=0; k<=3; k++ )
  {
    tk = 2 * k;
    tkp1 = tk + 1;
    zp = 1.0;
    ak = 0.0;
    bk = 0.0;
    for( s=0; s<=tk; s++ )
    {
      if( doa )
      {
        if( (s & 3) > 1 )
          sign = nflg;
        else
          sign = 1;
        ak += sign * mu[s] * zp * u[tk-s];
      }

      if( dob )
      {
        m = tkp1 - s;
        if( ((m+1) & 3) > 1 )
          sign = nflg;
        else
          sign = 1;
        bk += sign * lambda[s] * zp * u[m];
      }
      zp *= z32i;
    }

    if( doa )
    {
      ak *= np;
      t = fabs(ak);
      if( t < akl )
      {
        akl = t;
        pp += ak;
      }
      else
        doa = 0;
    }

    if( dob )
    {
      bk += lambda[tkp1] * zp * u[0];
      bk *= -np/sqz;
      t = fabs(bk);
      if( t < bkl )
      {
        bkl = t;
        qq += bk;
      }
      else
        dob = 0;
    }
#if DEBUG
    printf("a[%d] %.5E, b[%d] %.5E\n", k, ak, k, bk );
#endif
    if( np < MACHEP )
      break;
    np /= n*n;
  }

  /* normalizing factor ( 4*zeta/(1 - z**2) )**1/4	*/
  t = 4.0 * zeta/zz;
  t = std::sqrt( std::sqrt(t) );

  t *= ai*pp/cbrt(n)  +  aip*qq/(n23*n);
  return(t);
}

/* Asymptotic expansion for transition region,
* n large and x close to n.
* AMS55 #9.3.23.
*/

static double PF2[] = {
  -9.0000000000000000000e-2,
  8.5714285714285714286e-2
};
static double PF3[] = {
  1.3671428571428571429e-1,
  -5.4920634920634920635e-2,
  -4.4444444444444444444e-3
};
static double PF4[] = {
  1.3500000000000000000e-3,
  -1.6036054421768707483e-1,
  4.2590187590187590188e-2,
  2.7330447330447330447e-3
};
static double PG1[] = {
  -2.4285714285714285714e-1,
  1.4285714285714285714e-2
};
static double PG2[] = {
  -9.0000000000000000000e-3,
  1.9396825396825396825e-1,
  -1.1746031746031746032e-2
};
static double PG3[] = {
  1.9607142857142857143e-2,
  -1.5983694083694083694e-1,
  6.3838383838383838384e-3
};


static double jnt(double n, double x)
{
  double z, zz, z3;
  double cbn, n23, cbtwo;
  double ai, aip, bi, bip;	/* Airy functions */
  double nk, fk, gk, pp, qq;
  double F[5], G[4];
  int k;

  cbn = cbrt(n);
  z = (x - n)/cbn;
  cbtwo = cbrt( 2.0 );

  /* Airy function */
  zz = -cbtwo * z;
  airy( zz, &ai, &aip, &bi, &bip );

  /* polynomials in expansion */
  zz = z * z;
  z3 = zz * z;
  F[0] = 1.0;
  F[1] = -z/5.0;
  F[2] = polevl( z3, PF2, 1 ) * zz;
  F[3] = polevl( z3, PF3, 2 );
  F[4] = polevl( z3, PF4, 3 ) * z;
  G[0] = 0.3 * zz;
  G[1] = polevl( z3, PG1, 1 );
  G[2] = polevl( z3, PG2, 2 ) * z;
  G[3] = polevl( z3, PG3, 2 ) * zz;
#if DEBUG
  for( k=0; k<=4; k++ )
    printf( "F[%d] = %.5E\n", k, F[k] );
  for( k=0; k<=3; k++ )
    printf( "G[%d] = %.5E\n", k, G[k] );
#endif
  pp = 0.0;
  qq = 0.0;
  nk = 1.0;
  n23 = cbrt( n * n );

  for( k=0; k<=4; k++ )
  {
    fk = F[k]*nk;
    pp += fk;
    if( k != 4 )
    {
      gk = G[k]*nk;
      qq += gk;
    }
#if DEBUG
    printf("fk[%d] %.5E, gk[%d] %.5E\n", k, fk, k, gk );
#endif
    nk /= n23;
  }

  fk = cbtwo * ai * pp/cbn  +  cbrt(4.0) * aip * qq/n;
  return(fk);
}


//////////////////////////// end of    jv.c ///////////////////////


/***************************************************
***  END of the Bessel and Gamma functions.      ***
****************************************************/

