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

#include "error.h"
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include "callstack.h"

static void report(const char *prefix, const char *err, va_list params)
{
	char msg[1024];
	vsnprintf(msg, sizeof(msg), err, params);
	fprintf(stderr, "%s%s\n", prefix, msg);
}

static void report_w_loc(const char *prefix, const char *file, int line, const char *func,
                         const char *err, va_list params)
{
	char msg[2048];
	vsnprintf(msg, sizeof(msg), err, params);
	fprintf(stderr, "%s%s:%d: %s: %s\n", prefix, file, line, func, msg);
}

// errors
void h_exit(int line, const char *func, const char *file, char const *fmt, ...)
{
	va_list params;

	va_start(params, fmt);
	report_w_loc("FATAL: ", file, line, func, fmt, params);
	va_end(params);
	get_callstack().dump();
	exit(128);
}

// FIXME: this function should be removed (use the one from logging.cpp):
void error_function(const char *err, ...)
{
	va_list params;

	va_start(params, err);
	report("FATAL ERROR: ", err, params);
	va_end(params);
	exit(128);
}

void warning(const char *warn, ...)
{
	va_list params;

	va_start(params, warn);
	report("WARNING: ", warn, params);
	va_end(params);
}

void h_mem_check(int line, const char *func, const char *file, void *var)
{
        //va_list params;

	if (var == NULL) {
	  //report_w_loc("FATAL: ", file, line, func, "Out of memory.", params);
	  report_w_loc("FATAL: ", file, line, func, "Out of memory.", NULL);
		get_callstack().dump();
		exit(EXIT_FAILURE);
	}
}
