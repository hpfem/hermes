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
/*! \file c99_functions.h
    \brief File containing definitions from the C99 standard that are missing in MSVC.
*/
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
