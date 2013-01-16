// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public Licenserix
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __H2D_ELEMWISE_PARAMETER_H
#define __H2D_ELEMWISE_PARAMETER_H

#include "config.h"
#include "compat.h"
#include "common.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// Enum for the basic types. Upon the attribute of this enum's type, the assembling decides how to handle the particular parameter.
    enum ElemwiseParameterType
    {
      ElemwiseParameterTypeFunc,
      ElemwiseParameterTypeNonlinear
    };

    /// General parameter.
    template<typename Scalar>
    class HERMES_API ElemwiseParameter
    {
    public:
      /// Constructor.
      ElemwiseParameter();

      /// Destructor.
      virtual ~ElemwiseParameter();

      /// Get the type of this instance.
      virtual ElemwiseParameterType get_type() = 0;
    protected:
      template<typename T> friend class DiscreteProblem;
      template<typename T> friend class DiscreteProblemLinear;
    };
  }
}
#endif
