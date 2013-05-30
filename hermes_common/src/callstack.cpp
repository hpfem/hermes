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
/*! \file callstack.cpp
    \brief File containing functionality for investigating call stack.
*/
#include "callstack.h"
#include <signal.h>
#include <stdlib.h>

// Basically GNU stuff
#ifdef EXECINFO_FOUND
  #include <execinfo.h>
  void handler(int sig)
  {
    void *array[20];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, 20);

    // print out all the frames to stderr
    printf("Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, 2);
    exit(1);
  }

#endif

void CallStack::dump(int signalCode)
{
// Basically WIN stuff
#ifdef WITH_STACKTRACE
  #ifdef _WINDOWS
    MyStackWalker sw;
    sw.ShowCallstack();
  #else
    #ifdef EXECINFO_FOUND
      handler(signalCode);
    #endif
  #endif
#endif
}