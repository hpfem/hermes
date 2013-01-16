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

#ifndef __H2D_ELEMWISE_PARAMETER_NONLINEAR_H
#define __H2D_ELEMWISE_PARAMETER_NONLINEAR_H

#include "elemwise_parameter.h"
#include "../spline.h"
#include "../mixins2d.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// Nonlinear element-wise parameter.
    template<typename Scalar>
    class HERMES_API ElemwiseParameterNonlinear : public ElemwiseParameter<Scalar>
    {
    public:
      /// Constructor.
      ElemwiseParameterNonlinear();

      /// Destructor.
      virtual ~ElemwiseParameterNonlinear();

      /// Get the type of this instance.
      virtual ElemwiseParameterType get_type();
    protected:
      virtual Scalar get_value(Scalar u_ext_value) = 0;

      template<typename T> friend class DiscreteProblem;
      template<typename T> friend class DiscreteProblemLinear;
    };

    /// Spline-represented nonlinear parameter.
    class HERMES_API ElemwiseParameterSpline : public ElemwiseParameterNonlinear<double>, public Hermes2D::Mixins::StateQueryable
    {
    public:
      /// Constructor.
      ElemwiseParameterSpline(CubicSpline* spline);
      ElemwiseParameterSpline();

      /// Destructor.
      virtual ~ElemwiseParameterSpline();

      /// Ask if the instance is fine.
      virtual bool isOkay() const;

        /// Get class name, for the purpose of messaging.
      virtual std::string getClassName() const;
    private:
      virtual double get_value(double u_ext_value);

      CubicSpline* spline;

      template<typename T> friend class DiscreteProblem;
      template<typename T> friend class DiscreteProblemLinear;
    };
  }
}
#endif
