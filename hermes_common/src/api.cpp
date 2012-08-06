// This file is part of Hermes2D
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

#include "api.h"
#include <utility>
#include "callstack.h"
#include "common.h"
#include "exceptions.h"
#include "matrix.h"

namespace Hermes
{
  Api::Parameter::Parameter(int defaultVal)
  {
    this->defaultVal = defaultVal;
    this->userSet = false;
  }

  Api::Api()
  {
    signal(SIGABRT, CallStack::dump);
    signal(SIGFPE, CallStack::dump);
    signal(SIGILL, CallStack::dump);
    signal(SIGSEGV, CallStack::dump);
    signal(SIGTERM, CallStack::dump);

    this->parameters.insert(std::pair<HermesCommonApiParam, Parameter*> (Hermes::exceptionsPrintCallstack,new Parameter(0)));
    this->parameters.insert(std::pair<HermesCommonApiParam, Parameter*> (Hermes::matrixSolverType,new Parameter(SOLVER_UMFPACK)));
  }

  Api::~Api()
  {
    this->parameters.clear();
  }

  int Api::getParamValue(HermesCommonApiParam param)
  {
    if(this->parameters.find(param) == parameters.end())
      throw Hermes::Exceptions::Exception("Wrong Hermes::Api parameter name:%i", param);
    if(this->parameters.find(param)->second->userSet)
      return this->parameters.find(param)->second->userVal;
    else
      return this->parameters.find(param)->second->defaultVal;
  }

  void Api::setParamValue(HermesCommonApiParam param, int value)
  {
    if(this->parameters.find(param) == parameters.end())
      throw Hermes::Exceptions::Exception("Wrong Hermes::Api parameter name:%i", param);
    this->parameters.find(param)->second->userSet = true;
    this->parameters.find(param)->second->userVal = value;
  }

  Hermes::Api HermesCommonApi;
}