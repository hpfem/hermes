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
#include "third_party_codes/trilinos-teuchos/Teuchos_stacktrace.hpp"
#include <signal.h>
#include <stdlib.h>

/// Definition of the global CallStack instance.
CallStack callstack;

CallStackObj::CallStackObj(int ln, const char *func, const char *file)
{
}

CallStackObj::~CallStackObj()
{
}

static void sighandler(int signo)
{
}

void callstack_initialize()
{
}

void callstack_finalize()
{
}

CallStack::CallStack(int max_size)
{
}

CallStack::~CallStack()
{
}

void CallStack::dump()
{
}
const char * CallStack::getLastFunc()
{
  if (size>0)
    return stack[size-1]->func;
  else
    return NULL;
}
