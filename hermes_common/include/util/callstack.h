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

#include "util/compat.h"
#include "common.h"
#ifdef WITH_STACKTRACE
  #ifdef _WINDOWS
    #include <windows.h>
    #include "StackWalker.h"

    class MyStackWalker : public StackWalker
    {
    public:
      MyStackWalker() : StackWalker(StackWalker::RetrieveSymbol & StackWalker::RetrieveLine, NULL, GetCurrentProcessId(), GetCurrentProcess())
      {
        int a = 1;
      }
    protected:
      virtual void OnOutput(LPCSTR szText)
        { printf(szText);}
    };

  #endif
#endif

/// Call stack class.
class HERMES_API CallStack
{
public:
  // dump the call stack objects to standard error
  static void dump(int signalCode);
};
#endif