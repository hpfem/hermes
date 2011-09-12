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
        /// init exception with message
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
    /// Exception occurs when some parameter is Null and it shouldn't be.
    class NullException : public Exception
    {
      public:
        /// Constructor
        /// \param[in] paramnIdx index null parameter.
        NullException(int paramIdx);
        /// \return index of null parameter.
        int getParamIdx() const;
      private:
        int paramIdx;
    };
  }
}
#endif
