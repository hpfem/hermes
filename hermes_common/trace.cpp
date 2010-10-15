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

#include <cstdio>
#include <cstdarg>
#include <ctime>

#include "trace.h"

// VERBOSE

static int verbose_level = 0;

void set_verbose_level(int lvl) {
	verbose_level = lvl;
}

void verbose_printf(int level, char const *fmt, ...) {
	if (level <= verbose_level) {
		va_list ap;
		va_start(ap, fmt);
		vprintf(fmt, ap);
		va_end(ap);
	}
}


// DEBUGGING stuff

static int debug_level = 0;

void debug_output_on() {
	debug_level++;
}

void debug_output_off() {
	if (debug_level > 0)
		debug_level--;
}

void debug_printf(char const *fmt, ...) {
	if (debug_level > 0) {
		va_list ap;
		va_start(ap, fmt);
		vfprintf(stderr, fmt, ap);
		va_end(ap);
		fflush(stderr);
	}
}

// TRACING stuff

FILE *trace_file = NULL;
static int trace_level = 0;
static clock_t trace_time = clock();

void trace_on() {
	trace_level++;
}

void trace_off() {
	if (trace_level > 0)
		trace_level--;
}

// Initialize tracing to the file
void trace_start(const char *file_name) {
	trace_file = fopen(file_name, "w");
	if (trace_file == NULL) {
		perror("trace_start");
	}

	trace_on();
}

//
void trace_end() {
	fclose(trace_file);
	trace_file = NULL;
}

void trace(int line, const char *func, const char *file, char const *fmt, ...) {
	static bool warn_trace = false;		// we need to warn the user only once

	if (trace_level > 0) {
		if (trace_file != NULL) {
			va_list ap;
			va_start(ap, fmt);
			vfprintf(trace_file, fmt, ap);
			va_end(ap);
			fflush(stderr);
		}
		else {
			if (!warn_trace)
				fprintf(stderr, "Warning: TRACING not initialized (call TRACE_START at the beginning of your program).\n");
			warn_trace = true;
		}
	}
}
