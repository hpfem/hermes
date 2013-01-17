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

#ifndef __H2D_PARAMETER_ELEMWISE_H
#define __H2D_PARAMETER_ELEMWISE_H

#include "parameter.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// Enum for basic Element-wise Parameter types.
    enum ParameterElemwiseValueType
    {
      ParameterElemwiseFunctionValue = 1,
      ParameterElemwiseMeshFunctionValue = 2
    };

    /// General parameter.
    template<typename Scalar>
    class HERMES_API ParameterElemwise : public Parameter<Scalar>
    {
    public:
      /// Constructor.
      ParameterElemwise();

      /// Destructor.
      virtual ~ParameterElemwise();

    protected:
      virtual ParameterElemwiseValueType get_type() = 0;
      virtual Scalar get_value(double x, double y) = 0;

      template<typename T> friend class DiscreteProblem;
      template<typename T> friend class DiscreteProblemLinear;
    };
  }
}
#endif
