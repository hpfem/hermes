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
/*! \file exceptions.h
    \brief File containing definition of exceptions classes
*/
#ifndef __HERMES_COMMON_EXCEPTIONS_H_
#define __HERMES_COMMON_EXCEPTIONS_H_

#include<stdio.h>
#include<string>
#include"callstack.h"
#include<string.h>

namespace Hermes
{
  namespace Exceptions
  {
    /// \brief Exception interface
    class HERMES_API Exception
    {
      public:
        /// \brief Init exception with default message.
        Exception();
        /// Init exception with message.
        /// \param[in] msg message
        Exception(const char * msg);
        /// \brief print error message to stderr
        void printMsg() const;
        /// \brief get pointer to error message
        const char * getMsg() const;
        /// \return name of function where exception was created.
        const char * getFuncName() const;
        virtual ~Exception(){};
      protected:
        const char * message;
      private:
        const char * func;
    };

    /// \brief Null parameter exception.
    /// Exception occurs when some parameter is Null or empty and it shouldn't be.
    class HERMES_API NullException : public Exception
    {
      public:
        /// Constructor
        /// \param[in] paramnIdx index of null parameter.
        NullException(int paramIdx);
        /// Null item is passed in vector or array.
        /// \param[in] paramnIdx index of parameter.
        /// \param[in] elementIdx index of null item in array parameter.
        NullException(int paramIdx, int itemIdx);
        /// \return index of null parameter.
        int getParamIdx() const;
        /// \return index of null item in array parameter. Returns -1 if bad parrameter is not array with null item.
        int getItemIdx() const;
        ~NullException();
        NullException(const NullException & e);
      private:
        int paramIdx, itemIdx;
    };

    /// \brief Parameter length parameter exception.
    /// Exception occurs when some parameter has wrong length.
    class HERMES_API LengthException : public Exception
    {
      public:
        /// One parameter has wrong length.
        /// \param[in] paramnIdx index wrong parameter.
        /// \param[in] wrong actual length of parameter.
        /// \param[in] right right length of parameter.
        LengthException(int paramIdx, int wrong, int right);
        /// Two parameters should have same length and they dont have.
        /// \param[in] fstParamnIdx index first parameter.
        /// \param[in] sndParamnIdx index second parameter.
        /// \param[in] first actual length of first parameter.
        /// \param[in] second actual length of second parameter.
        LengthException(int fstParamIdx, int sndParamIdx, int first, int second);
        /// \return index of first wrong parameter.
        int getFirstParamIdx() const;
        /// \return index of second wrong parameter. Returns -1 when only one parameter is wrong.
        int getSecondParamIdx() const;
        /// \return length of first parameter.
        int getFirstLength() const;
        /// \return expected length of first parameter.
        int getExpectedLength() const;
        ~LengthException();
        LengthException(const LengthException & e);
      private:
        int fstParamIdx, sndParamIdx, wrong, right;
    };

    /// \brief Linear solver failed.
    class HERMES_API LinearMatrixSolverException : public Exception
    {
      public:
        /// \brief Linear solver failed from unknown reason.
        LinearMatrixSolverException();
        /// Linear solver failed from spevific reason.
        /// \param[in] reasen specification of solver fail.
        LinearMatrixSolverException(const char * reason);
        ~LinearMatrixSolverException();
        LinearMatrixSolverException(const LinearMatrixSolverException & e);
    };

    /// \brief Value is out of allowed range
    class HERMES_API ValueException : public Exception
    {
      public:
        /// Value is greather or lower than allowed.
        /// \param[in] name name of variable
        /// \param[in] value value of variable
        /// \param[in] allowed allowed value (maximum or minimum)
        ValueException(const char * name, double value, double allowed);
        /// Value is out of range.
        /// \param[in] name name of variable.
        /// \param[in] value value of variable.
        /// \param[in] min minimum allowed value.
        /// \param[in] max minimum allowed value.
        ValueException(const char * name, double value, double min, double max);
        /// String value is not supported.
        ValueException(const char * name, std::string passed);
        /// \return bad value of variable.
        double getValue() const;
        /// return allowed value of variable.
        double getAllowed() const;
        ~ValueException();
        ValueException(const ValueException & e);
      private:
        double value, allowed;
    };
  }
}
#endif
