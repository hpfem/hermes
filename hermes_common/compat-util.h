#ifndef COMPAT_UTIL_H
#define COMPAT_UTIL_H

#ifdef __GNUC__
#define NORETURN __attribute__((noreturn))
#else
#define NORETURN
#ifndef __attribute__
#define __attribute__(x)
#endif
#endif


#endif
