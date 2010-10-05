/*
   This header file contains platform compatibility layer.

   It is included from common.h, so it is automatically included in all hermes
   sources. The implementation of the functions in this file is in the
   src/compat directory.
*/

#ifndef __H2D_COMPAT_H
#define __H2D_COMPAT_H


#ifndef HAVE_FMEMOPEN
/// Implementation of GNU fmemopen. Intended to be used if the current platform does not support it.
FILE *fmemopen (void *buf, size_t size, const char *opentype);
#endif

//Windows DLL export/import definitions
#if defined(WIN32) || defined(_WINDOWS)
# if defined(_HERMESDLL)
#   define H2D_API __declspec(dllexport)
#   define H2D_API_USED_TEMPLATE(__implementation) template class H2D_API __implementation
#   define H2D_API_USED_STL_VECTOR(__type) H2D_API_USED_TEMPLATE(std::allocator<__type>); H2D_API_USED_TEMPLATE(std::vector<__type>)
# else
#   define H2D_API __declspec(dllimport)
#   define H2D_API_USED_TEMPLATE(__implementation)
//#   define H2D_API_USED_TEMPLATE(__implementation) extern template class H2D_API __implementation
#   define H2D_API_USED_STL_VECTOR(__type) H2D_API_USED_TEMPLATE(std::allocator<__type>); H2D_API_USED_TEMPLATE(std::vector<__type>)
# endif
#else
# define H2D_API
# define H2D_API_USED_TEMPLATE(__implementation)
# define H2D_API_USED_STL_VECTOR(__type)
#endif

//C99 functions
#include "compat/c99_functions.h"

#endif
