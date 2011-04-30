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

#include <stdarg.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include "utils.h"
#include "error.h"

/// Find the largest value
///
/// @return The largest value among numbers passed as arguments
/// @param[in] count - the number of values passed into this function (all of them are searched)
///
int maxn(int count, ...) {
	va_list ap;
	va_start(ap, count);
	int mx = INT_MIN;

	for (int i = 0; i < count; i++) {
		int num = va_arg(ap, int);
		if (num > mx)
			mx = num;
	}

	return mx;
}

void hermes_fwrite(const void* ptr, size_t size, size_t nitems, FILE* stream) {
	if (fwrite(ptr, size, nitems, stream) != nitems || ferror(stream))
		EXIT("Error writing to file: %s", strerror(ferror(stream)));
}

void hermes_fread(void* ptr, size_t size, size_t nitems, FILE* stream) {
	if (fread(ptr, size, nitems, stream) != nitems || ferror(stream))
		EXIT("Error reading file: %s", strerror(ferror(stream)));
}
