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
/*! \file callstack.h
    \brief File containing functionality for investigating call stack.
*/
#ifndef __HERMES_COMMON_CALLSTACK_H_
#define __HERMES_COMMON_CALLSTACK_H_

#include <stdio.h>
#include "compat.h"

// __PRETTY_FUNCTION__ missing on MSVC
#ifndef __GNUC__
  #define _F_ CallStackObj __call_stack_obj(__LINE__, __FUNCTION__, __FILE__);
#else
  #define _F_ CallStackObj __call_stack_obj(__LINE__, __PRETTY_FUNCTION__, __FILE__);
#endif

/// Holds data for one call stack object.
class HERMES_API CallStackObj
{
public:
	CallStackObj(int ln, const char *func, const char *file);
	~CallStackObj();

	int line;					// line number in the file
	const char *file;			// file
	const char *func;			// function name
};

/// Call stack class.
class HERMES_API CallStack
{
public:
	CallStack(int max_size = 32);
	~CallStack();

	// dump the call stack objects to standard error
	void dump();
  const char * getLastFunc();

protected:
	CallStackObj **stack;
	int size;
	int max_size;

	friend class CallStackObj;
};

/// Global call stack instance.
extern CallStack callstack;

#endif
