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

#include "common.h"
#include "compat.h"

namespace Hermes
{
  namespace Exceptions
  {
    /// \brief Exception interface
    class HERMES_API Exception : public std::exception
    {
      public:
        /// \brief Init exception with default message.
        Exception();
        /// Init exception with message.
        /// \param[in] msg message
        Exception(const char * msg, ...);
        /// \brief print error message to stderr
        void printMsg() const;
        /// \brief get pointer to error message
        virtual const char * what() const throw();
        /// \return name of function where exception was created.
        const char * getFuncName() const;
        virtual ~Exception() throw() {};

        virtual Exception* clone();
      protected:
        char * message;
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
        ~NullException() throw() {};
        NullException(const NullException & e);
        virtual Exception* clone();
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
        ~LengthException() throw() {};
        LengthException(const LengthException & e);
        virtual Exception* clone();
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
        /// \param[in] reason specification of solver fail.
        LinearMatrixSolverException(const char * reason);
        ~LinearMatrixSolverException() throw() {};
        LinearMatrixSolverException(const LinearMatrixSolverException & e);
        virtual Exception* clone();
    };

    /// \brief Numeric value is out of allowed range
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
        ~ValueException() throw() {};
        ValueException(const ValueException & e);
        virtual Exception* clone();
      private:
        double value, allowed;
    };

    /// \brief Linear solver failed.
    class HERMES_API FunctionNotOverridenException : public Exception
    {
      public:
        /// Constructor
        /// \param[in] name Name of the function.
        FunctionNotOverridenException(const char * msg, ...);
        ~FunctionNotOverridenException() throw() {};
        FunctionNotOverridenException(const FunctionNotOverridenException & e);
        virtual Exception* clone();
    };

    /// \brief Linear solver failed.
    class HERMES_API MeshLoadFailureException : public Exception
    {
      public:
        /// Constructor
        /// \param[in] name Name of the function.
        MeshLoadFailureException(const char * msg, ...);
        ~MeshLoadFailureException() throw() {};
        MeshLoadFailureException(const MeshLoadFailureException & e);
        virtual Exception* clone();
    };

    /// \brief Linear solver failed.
    class HERMES_API SpaceLoadFailureException : public Exception
    {
      public:
        /// Constructor
        /// \param[in] name Name of the function.
        SpaceLoadFailureException(const char * msg, ...);
        ~SpaceLoadFailureException() throw() {};
        SpaceLoadFailureException(const SpaceLoadFailureException & e);
        virtual Exception* clone();
    };

    /// \brief Linear solver failed.
    class HERMES_API SolutionSaveFailureException : public Exception
    {
      public:
        /// Constructor
        /// \param[in] name Name of the function.
        SolutionSaveFailureException(const char * msg, ...);
        ~SolutionSaveFailureException() throw() {};
        SolutionSaveFailureException(const SolutionSaveFailureException & e);
        virtual Exception* clone();
    };

    /// \brief Linear solver failed.
    class HERMES_API SolutionLoadFailureException : public Exception
    {
      public:
        /// Constructor
        /// \param[in] name Name of the function.
        SolutionLoadFailureException(const char * msg, ...);
        ~SolutionLoadFailureException() throw() {};
        SolutionLoadFailureException(const SolutionLoadFailureException & e);
        virtual Exception* clone();
    };
  }
}
#endif