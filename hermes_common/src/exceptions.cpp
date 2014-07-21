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
#include "util/memory_handling.h"
#include "mixins.h"

namespace Hermes
{
  namespace Exceptions
  {
    Exception::Exception() : std::exception()
    {
    }

    Exception::Exception(const char * msg, ...) : std::exception()
    {
      char text[10000];

      // print the message
      va_list arglist;
      va_start(arglist, msg);
      vsprintf(text, msg, arglist);
      va_end(arglist);

      this->message << text;
    }

    Exception::Exception(const Exception& e)
    {
      this->message << e.message.str();
    }

    void Exception::print_msg() const
    {
      std::cout << "Exception: " << message.str() << std::endl;
    }

    const char * Exception::what() const throw()
    {
      sprintf(const_cast<char*>(this->message_char), "%s", this->message.str().c_str());
      return const_cast<const char*>(this->message_char);
    }

    std::string Exception::info() const
    {
      return this->message.str();
    }

    IOException::IOException(ReadWrite readWrite, const char* filename_) : Exception(), readWrite(readWrite)
    {
      std::stringstream ss;
      ss << filename;
      this->filename = ss.str();
    }

    IOException::IOException(ReadWrite readWrite, std::string filename_) : Exception(), readWrite(readWrite), filename(filename_)
    {
    }

    IOException::IOException(const IOException & e)
    {
      this->readWrite = e.readWrite;
      this->filename = e.filename;
    }

    NullException::NullException() : Exception()
    {
      this->message << "A pointer is invalid (nullptr)";
    }

    NullException::NullException(unsigned int param_idx) : Exception()
    {
      this->param_idx = param_idx;
      this->item_idx = -1;
      this->message.clear();
      this->message << "Parameter number " << param_idx << " is nullptr";
    }

    NullException::NullException(unsigned int param_idx, unsigned int item_idx) : Exception()
    {
      this->param_idx = param_idx;
      this->item_idx = item_idx;
      this->message.clear();
      this->message << "Element number " << item_idx << " of parameter number " << param_idx << " is nullptr";
    }

    unsigned int NullException::get_param_idx() const
    {
      return param_idx;
    }

    unsigned int NullException::get_item_idx() const
    {
      return item_idx;
    }

    NullException::NullException(const NullException & e)
    {
      this->message.clear();
      this->message << e.message.str();
      param_idx = e.get_param_idx();
      item_idx = e.get_item_idx();
    }

    LengthException::LengthException() : Exception()
    {
      this->message.clear();
      this->message << "Two instances do not have the same length and they should";
    }

    LengthException::LengthException(unsigned int wrong, unsigned int right) : Exception()
    {
      this->wrong = wrong;
      this->right = right;
      this->message.clear();
      this->message << "An instance has length " << wrong << " and should have " << right;
    }

    LengthException::LengthException(unsigned int param_idx, unsigned int wrong, unsigned int right) : Exception()
    {
      fst_param_idx = param_idx;
      this->wrong = wrong;
      this->right = right;
      this->snd_param_idx = -1;
      this->message.clear();
      this->message << "Parameter number " << fst_param_idx << " have length " << wrong << " and should have " << right;
    }

    LengthException::LengthException(unsigned int fst_param_idx, unsigned int snd_param_idx, unsigned int first, unsigned int second) : Exception()
    {
      this->fst_param_idx = fst_param_idx;
      this->snd_param_idx = snd_param_idx;
      this->wrong = first;
      this->right = second;
      this->message.clear();
      this->message << "Parameter number " << fst_param_idx << " have length " << wrong << " and parameter number " << snd_param_idx << " have length " << right << " The lengths should be same.";
    }

    unsigned int LengthException::get_first_param_idx() const
    {
      return fst_param_idx;
    }

    unsigned int LengthException::get_second_param_idx() const
    {
      return snd_param_idx;
    }

    unsigned int LengthException::get_first_length() const
    {
      return wrong;
    }

    unsigned int LengthException::get_expected_length() const
    {
      return right;
    }

    LengthException::LengthException(const LengthException&e) : Exception()
    {
      this->message << e.message.str();
      this->fst_param_idx = e.get_first_param_idx();
      this->snd_param_idx = e.get_second_param_idx();
      this->wrong = e.get_first_length();
      this->right = e.get_expected_length();
    }

    LinearMatrixSolverException::LinearMatrixSolverException() : Exception()
    {
      this->message << "Linear solver failed.";
    }

    LinearMatrixSolverException::LinearMatrixSolverException(const char * reason, ...) : Exception()
    {
      char text[10000];

      // print the message
      va_list arglist;
      va_start(arglist, reason);
      vsprintf(text, reason, arglist);
      va_end(arglist);

      this->message << "Linear solver failed because: " << text;
    }

    LinearMatrixSolverException::LinearMatrixSolverException(const LinearMatrixSolverException&e) : Exception()
    {
      this->message << e.message.str();
    }

    ValueException::ValueException(const char * name, double value, double allowed) : Exception()
    {
      if (value > allowed)
        this->message << "Variable " << name << " is " << value << " but maximum allowed value is " << allowed;
      else
        this->message << "Variable " << name << " is " << value << " but minimum allowed value is " << allowed;
      this->value = value;
      this->allowed = allowed;
    }

    ValueException::ValueException(const char * name, double value, double min, double max) : Exception()
    {
      this->message << "Variable " << name << " is " << value << " allowed range is " << min << " -- " << max;
      this->value = value;
      if (value > min)
        this->allowed = max;
      else
        this->allowed = min;
    }

    ValueException::ValueException(const char * name, std::string passed) : Exception()
    {
      this->message << "Variable " << name << " does not support value " << passed;
    }

    double ValueException::get_value() const
    {
      return value;
    }

    double ValueException::get_allowed() const
    {
      return allowed;
    }

    ValueException::ValueException(const ValueException&e)
    {
      this->message << e.message.str();
      this->value = e.get_value();
      this->allowed = e.get_allowed();
    }

    MethodNotOverridenException::MethodNotOverridenException(const char * name, ...) : Exception()
    {
      char text[10000];

      // print the message
      va_list arglist;
      va_start(arglist, name);
      vsprintf(text, name, arglist);
      va_end(arglist);

      this->message << "Method not overriden: " << text;
    }

    MethodNotOverridenException::MethodNotOverridenException(const MethodNotOverridenException&e)
    {
      this->message << e.message.str();
    }

    MethodNotImplementedException::MethodNotImplementedException(const char * name, ...) : Exception()
    {
      char text[10000];

      // print the message
      va_list arglist;
      va_start(arglist, name);
      vsprintf(text, name, arglist);
      va_end(arglist);

      this->message << "Sorry, method not implemented so far: " << text;
    }

    MethodNotImplementedException::MethodNotImplementedException(const MethodNotImplementedException&e)
    {
      this->message << e.message.str();
    }

    MeshLoadFailureException::MeshLoadFailureException(const char * reason, ...) : Exception()
    {
      char text[10000];

      // print the message
      va_list arglist;
      va_start(arglist, reason);
      vsprintf(text, reason, arglist);
      va_end(arglist);

      this->message << "Mesh loading failed: " << (std::string)text;
    }

    MeshLoadFailureException::MeshLoadFailureException(const MeshLoadFailureException&e)
    {
      this->message << e.message.str();
    }

    SpaceLoadFailureException::SpaceLoadFailureException(const char * reason, ...) : Exception()
    {
      char text[10000];

      // print the message
      va_list arglist;
      va_start(arglist, reason);
      vsprintf(text, reason, arglist);
      va_end(arglist);

      this->message << "Space loading failed: " << text;
    }

    SpaceLoadFailureException::SpaceLoadFailureException(const SpaceLoadFailureException&e)
    {
      this->message << e.message.str();
    }

    SolutionSaveFailureException::SolutionSaveFailureException(const char * reason, ...) : Exception()
    {
      char text[10000];

      // print the message
      va_list arglist;
      va_start(arglist, reason);
      vsprintf(text, reason, arglist);
      va_end(arglist);

      this->message << "Solution saving failed: " << text;
    }

    SolutionSaveFailureException::SolutionSaveFailureException(const SolutionSaveFailureException&e)
    {
      this->message << e.message.str();
    }

    SolutionLoadFailureException::SolutionLoadFailureException(const char * reason, ...) : Exception()
    {
      char text[10000];

      // print the message
      va_list arglist;
      va_start(arglist, reason);
      vsprintf(text, reason, arglist);
      va_end(arglist);

      this->message << "Solution loading failed: " << text;
    }

    SolutionLoadFailureException::SolutionLoadFailureException(const SolutionLoadFailureException&e)
    {
      this->message << e.message.str();
    }
  }
}