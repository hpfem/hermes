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

#include"exceptions.h"

namespace Hermes
{
  namespace Exceptions
  {

    Exception::Exception()
    {
      message = NULL;
      func = callstack.getLastFunc();
    }

    Exception::Exception(const char * msg)
    {
      message = msg;
      func = callstack.getLastFunc();
    }

    void Exception::printMsg() const
    {
      if (message)
        fprintf(stderr, "%s\n", message);
      else
        fprintf(stderr, "Default exception\n");
      if (func)
        fprintf(stderr, "in %s\n", func);
    }

    const char * Exception::getFuncName() const
    {
      return func;
    }

    const char * Exception::getMsg() const
    {
      return message;
    }

    NullException::NullException(int paramIdx)
    {
      this->paramIdx = paramIdx;
      this->itemIdx = -1;
      char * msg = new char[27];
      sprintf(msg, "Parameter number %d is NULL", paramIdx);
      message = msg;
    }

    NullException::NullException(int paramIdx, int itemIdx)
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

    NullException::~NullException()
    {
      delete[] message;
    }

    NullException::NullException(const NullException & e)
    {
      char * msg= new char[strlen(e.getMsg())+1];
      strcpy(msg, e.getMsg());
      message=msg;
      paramIdx=e.getParamIdx();
      itemIdx=e.getItemIdx();
    }

    LengthException::LengthException(int paramIdx, int wrong, int right)
    {
      fstParamIdx = paramIdx;
      this->wrong = wrong;
      this->right = right;
      this->sndParamIdx = -1;
      char * msg = new char[60];
      sprintf(msg, "Parameter number %d have length %d and should have %d", fstParamIdx, wrong, right);
      message = msg;
    }

    LengthException::LengthException(int fstParamIdx, int sndParamIdx, int first, int second)
    {
      this->fstParamIdx=fstParamIdx;
      this->sndParamIdx=sndParamIdx;
      this->wrong=first;
      this->right=second;
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

    LengthException::~LengthException()
    {
      delete[]message;
    }

    LengthException::LengthException(const LengthException&e)
    {
      char * msg= new char[strlen(e.getMsg())+1];
      strcpy(msg, e.getMsg());
      message=msg;
      this->fstParamIdx=e.getFirstParamIdx();
      this->sndParamIdx=e.getSecondParamIdx();
      this->wrong=e.getFirstLength();
      this->right=e.getExpectedLength();
    }

    LinearMatrixSolverException::LinearMatrixSolverException()
    {
      char * msg =  new char[22];
      sprintf(msg, "Linear solver failed.");
      message = msg;
    }

    LinearMatrixSolverException::LinearMatrixSolverException(const char * reason)
    {
      char * msg =  new char[34 + strlen(reason)];
      sprintf(msg, "Linear solver failed because:\"%s\"", reason);
      message = msg;
    }

    LinearMatrixSolverException::~LinearMatrixSolverException()
    {
      delete[] message;
    }

    LinearMatrixSolverException::LinearMatrixSolverException(const LinearMatrixSolverException&e)
    {
      char * msg= new char[strlen(e.getMsg())+1];
      strcpy(msg, e.getMsg());
      message=msg;
    }
    
    ValueException::ValueException(const char * name, double value, double allowed)
    {
      char * msg =  new char[55 + strlen(name)];
      if (value>allowed)
        sprintf(msg, "Variable %s is %f but maximum allowed value is %f", name, value, allowed);
      else
        sprintf(msg, "Variable %s is %f but minimum allowed value is %f", name, value, allowed);
      message = msg;
      this->value = value;
      this->allowed = allowed;
    }

    ValueException::ValueException(const char * name, double value, double min, double max)
    {
      char * msg= new char[70+strlen(name)];
      sprintf(msg, "Variable %s is %f allowed range is %f -- %f", name, value, min, max);
      message=msg;
      this->value=value;
      if (value>min)
        this->allowed = max;
      else
        this->allowed = min;
    }


    ValueException::ValueException(const char * name, std::string passed)
    {
      char * msg= new char[70+strlen(name)];
      sprintf(msg, "Variable %s does not support value %s.", name, passed.c_str());
      message=msg;
    }

    double ValueException::getValue() const
    {
      return value;
    }

    double ValueException::getAllowed() const
    {
      return allowed;
    }

    ValueException::~ValueException()
    {
      delete[] message;
    }

    ValueException::ValueException(const ValueException&e)
    {
      char * msg= new char[strlen(e.getMsg())+1];
      strcpy(msg, e.getMsg());
      message=msg;
      this->value=e.getValue();
      this->allowed=e.getAllowed();
    }

  }
}
