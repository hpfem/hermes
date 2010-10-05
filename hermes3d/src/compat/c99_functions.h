#ifndef __H2D_C99_FUNCTIONS_H
#define __H2D_C99_FUNCTIONS_H

#ifdef IMPLELENT_C99

/* Definitions of C99 specification. Used in a case of MSVC 2008 and
 * below because MSVC follows C++ rather than C
 */

// Not-a-number constant.
extern H2D_API const double NAN;

// functions
extern H2D_API double exp2(double x); ///< exp 2
extern H2D_API double log2(double x); ///< log 2
extern H2D_API double cbrt(double x); ///< cubic root

#endif /* IMPLELENT_C99 */

#endif
