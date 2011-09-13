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

namespace Hermes
{
  namespace Exceptions
  {
    /// \brief Exception interface
    class Exception
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
      protected:
        const char * message;

    };

    /// \brief Null parameter exception.
    /// Exception occurs when some parameter is Null or empty and it shouldn't be.
    class NullException : public Exception
    {
      public:
        /// Constructor
        /// \param[in] paramnIdx index null parameter.
        NullException(int paramIdx);
        /// \return index of null parameter.
        int getParamIdx() const;
        ~NullException();
      private:
        int paramIdx;
    };

    /// \brief Parameter length parameter exception.
    /// Exception occurs when some parameter has wrong length.
    class LengthException : public Exception
    {
      public:
        /// One parameter has wrong length.
        /// \param[in] paramnIdx index wrong parameter.
        /// \param[in] wrong actual length of parameter.
        /// \param[in] right right length of parameter.
        LengthException(int paramIdx,int wrong, int right);
        /// Two parameters should have same length and they dont have.
        /// \param[in] fstParamnIdx index first parameter.
        /// \param[in] sndParamnIdx index second parameter.
        /// \param[in] first actual length of first parameter.
        /// \param[in] second actual length of second parameter.
        LengthException(int fstParamIdx, int sndParmIdx, int first, int second);
        /// \return index of first wrong parameter.
        int getFirstParamIdx() const;
        /// \return index of second wrong parameter. Returns -1 when only one parameter is wrong.
        int getSecondParamIdx() const;
        /// \return length of first parameter.
        int getFirstLength() const;
        /// \return expected length of first parameter.
        int getExpectedLength() const;
        ~LengthException();
      private:
        int fstParamIdx,sndParamIdx,wrong,right;
    };
  }
}
#endif
