// This file is part of HermesCommon
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
/*! \file compat.h
\brief File containing platform compatibility layer, especially for Win / MSVC.
*/

#ifndef __HERMES_COMMON_COMPAT_H
#define __HERMES_COMMON_COMPAT_H
#include <stdio.h>
#include "stddef.h"
#include "config.h"
#include <cstddef>

#ifndef HAVE_FMEMOPEN
/// Implementation of GNU fmemopen. Intended to be used if the current platform does not support it.
FILE *fmemopen (void *buf, size_t size, const char *opentype);
#endif

// Windows DLL export/import definitions.
#if defined (HERMES_STATIC_LIBS)
  #define HERMES_API
#else
  #if defined (AGROS)
    #if defined (HERMES_FOR_AGROS)
      #define HERMES_API __declspec(dllexport)
    #else
      #define HERMES_API __declspec(dllimport)
    #endif
  #else
    #if defined(WIN32) || defined(_WINDOWS)
      // Visual Studio 2010.
      #if defined(EXPORT_HERMES_DLL)
      // when building DLL (target project defines this macro)
        #define HERMES_API __declspec(dllexport)
      #else  
      // when using the DLL by a client project
      #define HERMES_API __declspec(dllimport)
      #endif
    #else 
      #define HERMES_API
    #endif
  #endif
#endif

#if defined(WIN32) || defined(_WINDOWS)
#ifdef HERMES_COMMON
  #define HERMES_COMMON_API __declspec(dllexport)
#else
  #define HERMES_COMMON_API __declspec(dllimport)
#endif
#else
  #define HERMES_COMMON_API 
#endif

#ifndef HAVE_STRCASECMP
#define strcasecmp strcmp
#endif

#ifdef __GNUC__
  #define NORETURN __attribute__((noreturn))
#else
  #define NORETURN
  #ifndef __attribute__
    #define __attribute__(x)
  #endif
#endif

// If C++ 11 is not supported
namespace std
{
#if defined(__GNUC__) && ((__GNUC__ == 4 && __GNUC_MINOR__ >= 6) || (__GNUC__ >= 5))
# define HACK_GCC_ITS_CPP0X 1
#endif
#if defined(nullptr_t) || (__cplusplus >= 199711L) || defined(HACK_GCC_ITS_CPP0X)
#else
#define nullptr NULL
#endif
}

#endif
