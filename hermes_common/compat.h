/*
   This header file contains platform compatibility layer.

   It is included from common.h, so it is automatically included in all hermes
   sources. The implementation of the functions in this file is in the
   src/compat directory.
*/

#ifndef __HERMES_COMPAT_H
#define __HERMES_COMPAT_H
#include <stdio.h>

#ifndef HAVE_FMEMOPEN
/// Implementation of GNU fmemopen. Intended to be used if the current platform does not support it.
FILE *fmemopen (void *buf, size_t size, const char *opentype);
#endif

// Windows DLL export/import definitions
// Visual Studio 10 throws an error : 
// Compiler Error C2252 - cannot explicitly instantiate template in current scope
// if attempted to instantiate the templates => only for older versions.
#if (defined(WIN32) || defined(_WINDOWS)) && _MSC_VER < 1600
  #if defined(EXPORT_HERMES_DLL)
  // when building DLL (target project defines this macro)
    #define HERMES_API __declspec(dllexport)
    #define HERMES_API_USED_TEMPLATE(__implementation) template class HERMES_API __implementation
  #else  
  // when using the DLL by a client project
    #define HERMES_API __declspec(dllimport)
    #define HERMES_API_USED_TEMPLATE(__implementation)
    //#define HERMES_API_USED_TEMPLATE(__implementation) extern template class HERMES_API __implementation
  #endif
  #define HERMES_API_USED_STL_VECTOR(__type) HERMES_API_USED_TEMPLATE(std::allocator<__type>); HERMES_API_USED_TEMPLATE(std::vector<__type>)
  
#else 

  #define HERMES_API
  #define HERMES_API_USED_TEMPLATE(__implementation)
  #define HERMES_API_USED_STL_VECTOR(__type)

#endif

#ifndef HAVE_STRCASECMP
#define strcasecmp strcmp
#endif

//C99 functions
#include "compat/c99_functions.h"

#ifdef __GNUC__
#define NORETURN __attribute__((noreturn))
#else
#define NORETURN
#ifndef __attribute__
#define __attribute__(x)
#endif
#endif

#endif
