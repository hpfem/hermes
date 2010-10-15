#ifndef __H3D_C99_FUNCTIONS_H
#define __H3D_C99_FUNCTIONS_H

#ifdef IMPLEMENT_C99

/* Definitions of C99 specification. Used in a case of MSVC 2008 and
 * below because MSVC follows C++ rather than C
 */

// Not-a-number constant.
extern H3D_API const double NAN;

// functions
extern H3D_API double exp2(double x); ///< exp 2
extern H3D_API double log2(double x); ///< log 2
extern H3D_API double cbrt(double x); ///< cubic root

#endif /* IMPLEMENT_C99 */

#endif
