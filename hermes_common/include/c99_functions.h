#ifndef __HERMES_COMMON_C99_FUNCTIONS_H
#define __HERMES_COMMON_C99_FUNCTIONS_H

#ifdef IMPLEMENT_C99

/* Definitions of C99 specification. Used in a case of MSVC 2008 and
 * below because MSVC follows C++ rather than C
 */

// Not-a-number constant.
#define NAN 0x7fffffffffffffffL;

// functions
extern HERMES_API double exp2(double x); ///< exp 2
extern HERMES_API double log2(double x); ///< log 2
extern HERMES_API double cbrt(double x); ///< cubic root

#endif /* IMPLEMENT_C99 */

#endif
