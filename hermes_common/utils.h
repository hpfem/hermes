// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#ifndef _UTILS_H_
#define _UTILS_H_
#include "compat.h"

//
// Miscelaneous utils
//

int HERMES_API maxn(int count, ...);

inline int max(int a, int b) {
	return a > b ? a : b;
}

//void hermes_fwrite(const void* ptr, size_t size, size_t nitems, FILE* stream);
//void hermes_fread(void* ptr, size_t size, size_t nitems, FILE* stream);

#define countof(a) 								(sizeof(a)/sizeof(a[0]))

#endif
