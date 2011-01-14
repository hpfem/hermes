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

#include "callstack.h"
#include "third_party_codes/trilinos-teuchos/Teuchos_stacktrace.hpp"
#include <signal.h>
#include <stdlib.h>

// global instance of the call stack object
static CallStack callstack;

// Call Stack Object ////

CallStackObj::CallStackObj(int ln, const char *func, const char *file) {
	this->line = ln;
	this->func = func;
	this->file = file;

	// add this object to the call stack
	if (callstack.size < callstack.max_size) {
		callstack.stack[callstack.size] = this;
		callstack.size++;
	}
}

CallStackObj::~CallStackObj() {
	// remove the object only if it is on the top of the call stack
	if (callstack.size > 0 && callstack.stack[callstack.size - 1] == this) {
		callstack.size--;
		callstack.stack[callstack.size] = NULL;
	}
}

// Signals ////

static
void sighandler(int signo) {
	const char *sig_name[64];

	sig_name[SIGABRT] = "Abort";
	sig_name[SIGSEGV] = "Segmentation violation";

	fprintf(stderr, "Caught signal %d (%s)\n", signo, sig_name[signo]);
	callstack.dump();
	exit(EXIT_FAILURE);
}

// Comment this out stop using Teuchos stacktrace (in that case the stacktrace
// code originally in h3d will be used). Teuchos stacktrace not used for WIN32 (execinfo.h and cxxabi.h absent).
#ifndef _WIN32
  #define HERMES_USE_TEUCHOS_STACKTRACE
#endif

void callstack_initialize() {
	// install our signal handlers
#ifdef HERMES_USE_TEUCHOS_STACKTRACE
    Teuchos::print_stack_on_segfault();
#else
	signal(SIGSEGV, sighandler);
	signal(SIGABRT, sighandler);
#endif
}

void callstack_finalize() {
}

// Call Stack ////

CallStack &get_callstack() { return callstack; }

CallStack::CallStack(int max_size) {
	this->max_size = max_size;
	this->size = 0;
	this->stack = new CallStackObj *[max_size];

	// initialize signals
	callstack_initialize();
}

CallStack::~CallStack() {
	delete [] stack;
}

void CallStack::dump() {
	if (size > 0) {
		fprintf(stderr, "Call stack:\n");
		for (int i = size - 1; i >= 0; i--)
			fprintf(stderr, "  %s:%d: %s\n", stack[i]->file, stack[i]->line, stack[i]->func);
	}
	else {
		fprintf(stderr, "No call stack available.\n");
	}
}
