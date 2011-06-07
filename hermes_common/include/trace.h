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
/*! \file trace.h
    \brief Library output.
*/
#ifndef __HERMES_COMMON_TRACE_H_
#define __HERMES_COMMON_TRACE_H_

//
// This is used for library output
//
// Use VERBOSE macro for library outputs (behaves like printf, but the first parameter is level).
// The output is printed if the output level is smaller than current verbose level.
// Use SET_VERBOSE_LEVEL to change the verbose level. This macro is not intended for the use within
// the library. Inspite of that, VERBOSE macro should be used only inside the library.
//
// There are the following levels of verboseness:
//   0 - nothing is printed (silent mode)
//   1 - only basic things are displayed
//   3 - quite verbose
//   5 - really verbose
//   7 - psychotic output (enjoy reading the outputs :-)
//
//  Levels 2, 4, 6 are reserved (for now)
//
#define VERBOSE(lvl, fmt, ...) verbose_printf(lvl, fmt, ## __VA_ARGS__);
#define SET_VERBOSE_LEVEL(lvl) set_verbose_level(lvl);

void verbose_printf(int level, char const *fmt, ...);
void set_verbose_level(int lvl);

//
// Use this for debugging purposes
//
// DEBUG_PRINT behaves like printf and prints on stderr
// Debugging output can be disabled by DEBUG_OUTPUT_OFF and re-enabled by DEBUG_OUTPUT_ON
//
#ifdef DEBUG
	#define DEBUG_PRINT(fmt, ...) debug_printf(fmt, ## __VA_ARGS__)
	#define DEBUG_OUTPUT_ON debug_output_on();
	#define DEBUG_OUTPUT_OFF debug_output_off();

	void debug_output_on();
	void debug_output_off();
	void debug_printf(char const *fmt, ...);
#else
	#define DEBUG_PRINT(fmt, ...)
	#define DEBUG_OUTPUT_ON
	#define DEBUG_OUTPUT_OFF
#endif

//
// Use this for tracing (for long term runs to see what was happening in the code)
//
// Tracing is done into a file.
// First, call TRACE_START with the name of a file where you want the output.
// Use TRACE_ON and TRACE_OFF to turn on or off the tracing output
// Use TRACE to output something (behaves like printf but into a file specified by TRACE_START)
// At the end, call TRACE_END to close the tracing file
//
// TRACE_START and TRACE_END should be called outside the library (they should not be used inside the library)
//
#ifdef TRACING
	#define TRACE_START(fn) trace_start(fn)
	#define TRACE_END trace_end()

	#define TRACE_ON trace_on()
	#define TRACE_OFF trace_off()

	#define FUNC_BEGIN(...) func_begin(__LINE__, __PRETTY_FUNCTION__, __FILE__, ## __VA_ARGS__)
	#define FUNC_END(...) func_end(__LINE__, __PRETTY_FUNCTION__, __FILE__, ## __VA_ARGS__)
	#define TRACE(...) trace(__LINE__, __PRETTY_FUNCTION__, __FILE__, ## __VA_ARGS__)

	void func_begin(int line, const char *func, const char *file, char const *fmt, ...);
	void func_end(int line, const char *func, const char *file, char const *fmt, ...);
	void trace(int line, const char *func, const char *file, char const *fmt, ...);
	void trace_start(const char *file_name);
	void trace_end();

	void trace_on();
	void trace_off();
#else
	#define TRACE_START(fn)
	#define TRACE_END

	#define TRACE_ON
	#define TRACE_OFF

	#define FUNC_BEGIN(...)
	#define FUNC_END(...)
	#define TRACE(...)
#endif

// TODO: implement this
#define INFO(...)

#endif
