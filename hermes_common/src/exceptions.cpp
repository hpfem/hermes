// This file is part of HermesCommon
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http:// hpfem.org/.
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

#include "exceptions.h"
#include <string>
#include "api.h"
#include "callstack.h"

namespace Hermes
{
  namespace Exceptions
  {
    Exception::Exception() : std::exception(), message(new char[1000])
    {
    }

    Exception::Exception(const char * msg, ...) : std::exception(), message(new char[1000])
    {
      char text[1024];

      // print the message
      va_list arglist;
      va_start(arglist, msg);
      vsprintf(text, msg, arglist);
      va_end(arglist);

      strcpy(message, text);
    }

    void Exception::printMsg() const
    {
      if(message)
        printf("Exception: %s\n", message);
      else
        printf("Default exception\n");
      if(Hermes::HermesCommonApi.getParamValue(Hermes::exceptionsPrintCallstack) == 1)
        CallStack::dump(0);
    }

    Exception* Exception::clone()
    {
      return new Exception(*this);
    }

    const char * Exception::what() const throw()
    {
      return message;
    }

    NullException::NullException(int paramIdx) : Exception()
    {
      this->paramIdx = paramIdx;
      this->itemIdx = -1;
      char * msg = new char[27];
      sprintf(msg, "Parameter number %d is NULL", paramIdx);
      message = msg;
    }

    NullException::NullException(int paramIdx, int itemIdx) : Exception()
    {
      this->paramIdx = paramIdx;
      this->itemIdx = itemIdx;
      char * msg = new char[55];
      sprintf(msg, "Element number %d of parameter number %d is NULL", itemIdx, paramIdx);
      message = msg;
    }

    int NullException::getParamIdx() const
    {
      return paramIdx;
    }

    int NullException::getItemIdx() const
    {
      return itemIdx;
    }

    NullException::NullException(const NullException & e)
    {
      char * msg = new char[strlen(e.what())+1];
      strcpy(msg, e.what());
      message = msg;
      paramIdx = e.getParamIdx();
      itemIdx = e.getItemIdx();
    }

    Exception* NullException::clone()
    {
      return new NullException(*this);
    }

    LengthException::LengthException(int paramIdx, int wrong, int right) : Exception()
    {
      fstParamIdx = paramIdx;
      this->wrong = wrong;
      this->right = right;
      this->sndParamIdx = -1;
      char * msg = new char[60];
      sprintf(msg, "Parameter number %d have length %d and should have %d", fstParamIdx, wrong, right);
      message = msg;
    }

    LengthException::LengthException(int fstParamIdx, int sndParamIdx, int first, int second) : Exception()
    {
      this->fstParamIdx = fstParamIdx;
      this->sndParamIdx = sndParamIdx;
      this->wrong = first;
      this->right = second;
      char * msg = new char[110];
      sprintf(msg, "Parameter number %d have length %d and parameter number %d have length %d. The lengths should be same",
            fstParamIdx, wrong, sndParamIdx, right);
      message = msg;
    }

    int LengthException::getFirstParamIdx() const
    {
      return fstParamIdx;
    }

    int LengthException::getSecondParamIdx() const
    {
      return sndParamIdx;
    }

    int LengthException::getFirstLength() const
    {
      return wrong;
    }

    int LengthException::getExpectedLength() const
    {
      return right;
    }

    LengthException::LengthException(const LengthException&e) : Exception()
    {
      char * msg = new char[strlen(e.what())+1];
      strcpy(msg, e.what());
      message = msg;
      this->fstParamIdx = e.getFirstParamIdx();
      this->sndParamIdx = e.getSecondParamIdx();
      this->wrong = e.getFirstLength();
      this->right = e.getExpectedLength();
    }

    Exception* LengthException::clone()
    {
      return new LengthException(*this);
    }

    LinearMatrixSolverException::LinearMatrixSolverException() : Exception()
    {
      char * msg =  new char[22];
      sprintf(msg, "Linear solver failed.");
      message = msg;
    }

    LinearMatrixSolverException::LinearMatrixSolverException(const char * reason) : Exception()
    {
      char * msg =  new char[34 + strlen(reason)];
      sprintf(msg, "Linear solver failed because:\"%s\"", reason);
      message = msg;
    }

    LinearMatrixSolverException::LinearMatrixSolverException(const LinearMatrixSolverException&e) : Exception()
    {
      char * msg = new char[strlen(e.what())+1];
      strcpy(msg, e.what());
      message = msg;
    }

    Exception* LinearMatrixSolverException::clone()
    {
      return new LinearMatrixSolverException(*this);
    }

    ValueException::ValueException(const char * name, double value, double allowed) : Exception()
    {
      char * msg =  new char[55 + strlen(name)];
      if(value>allowed)
        sprintf(msg, "Variable %s is %f but maximum allowed value is %f", name, value, allowed);
      else
        sprintf(msg, "Variable %s is %f but minimum allowed value is %f", name, value, allowed);
      message = msg;
      this->value = value;
      this->allowed = allowed;
    }

    ValueException::ValueException(const char * name, double value, double min, double max) : Exception()
    {
      char * msg = new char[70+strlen(name)];
      sprintf(msg, "Variable %s is %f allowed range is %f -- %f", name, value, min, max);
      message = msg;
      this->value = value;
      if(value>min)
        this->allowed = max;
      else
        this->allowed = min;
    }

    ValueException::ValueException(const char * name, std::string passed) : Exception()
    {
      char * msg = new char[70+strlen(name)];
      sprintf(msg, "Variable %s does not support value %s.", name, passed.c_str());
      message = msg;
    }

    double ValueException::getValue() const
    {
      return value;
    }

    double ValueException::getAllowed() const
    {
      return allowed;
    }

    ValueException::ValueException(const ValueException&e)
    {
      char * msg = new char[strlen(e.what())+1];
      strcpy(msg, e.what());
      message = msg;
      this->value = e.getValue();
      this->allowed = e.getAllowed();
    }

    Exception* ValueException::clone()
    {
      return new ValueException(*this);
    }

    FunctionNotOverridenException::FunctionNotOverridenException(const char * name, ...) : Exception()
    {
      char text[1024];

      // print the message
      va_list arglist;
      va_start(arglist, name);
      vsprintf(text, name, arglist);
      va_end(arglist);

      message = text;
    }

    FunctionNotOverridenException::FunctionNotOverridenException(const FunctionNotOverridenException&e)
    {
      char * msg = new char[strlen(e.what())+1];
      strcpy(msg, e.what());
      message = msg;
    }

    Exception* FunctionNotOverridenException::clone()
    {
      return new FunctionNotOverridenException(*this);
    }

    MeshLoadFailureException::MeshLoadFailureException(const char * reason, ...) : Exception()
    {
      char * text = new char[strlen(reason)+1];

      // print the message
      va_list arglist;
      va_start(arglist, reason);
      vsprintf(text, reason, arglist);
      va_end(arglist);

      message = text;
    }

    MeshLoadFailureException::MeshLoadFailureException(const MeshLoadFailureException&e)
    {
      char * msg = new char[strlen(e.what())+1];
      strcpy(msg, e.what());
      message = msg;
    }

    Exception* MeshLoadFailureException::clone()
    {
      return new MeshLoadFailureException(*this);
    }

    SpaceLoadFailureException::SpaceLoadFailureException(const char * reason, ...) : Exception()
    {
      char * text = new char[strlen(reason)+1];

      // print the message
      va_list arglist;
      va_start(arglist, reason);
      vsprintf(text, reason, arglist);
      va_end(arglist);

      message = text;
    }

    SpaceLoadFailureException::SpaceLoadFailureException(const SpaceLoadFailureException&e)
    {
      char * msg = new char[strlen(e.what())+1];
      strcpy(msg, e.what());
      message = msg;
    }

    Exception* SpaceLoadFailureException::clone()
    {
      return new SpaceLoadFailureException(*this);
    }

    SolutionSaveFailureException::SolutionSaveFailureException(const char * reason, ...) : Exception()
    {
      char * text = new char[strlen(reason)+1];

      // print the message
      va_list arglist;
      va_start(arglist, reason);
      vsprintf(text, reason, arglist);
      va_end(arglist);

      message = text;
    }

    SolutionSaveFailureException::SolutionSaveFailureException(const SolutionSaveFailureException&e)
    {
      char * msg = new char[strlen(e.what())+1];
      strcpy(msg, e.what());
      message = msg;
    }

    Exception* SolutionSaveFailureException::clone()
    {
      return new SolutionSaveFailureException(*this);
    }

    SolutionLoadFailureException::SolutionLoadFailureException(const char * reason, ...) : Exception()
    {
      char * text = new char[strlen(reason)+1];

      // print the message
      va_list arglist;
      va_start(arglist, reason);
      vsprintf(text, reason, arglist);
      va_end(arglist);

      message = text;
    }

    SolutionLoadFailureException::SolutionLoadFailureException(const SolutionLoadFailureException&e)
    {
      char * msg = new char[strlen(e.what())+1];
      strcpy(msg, e.what());
      message = msg;
    }
    
    Exception* SolutionLoadFailureException::clone()
    {
      return new SolutionLoadFailureException(*this);
    }
  }
}