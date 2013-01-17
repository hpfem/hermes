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
#include "../function/solution.h"
#include "../mixins2d.h"
#include "../mesh/mesh.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// Nonlinear parameter value type.
    enum ElemwiseParameterNonlinearValueType
    {
      ElemwiseParameterNonlinearValue,
      ElemwiseParameterNonlinearDerivative
    };

    /// Nonlinear element-wise parameter.
    template<typename Scalar>
    class HERMES_API ElemwiseParameterNonlinear : public ElemwiseParameter<Scalar>
    {
    public:
      /// Constructor.
      ElemwiseParameterNonlinear(ElemwiseParameterNonlinearValueType value_type = ElemwiseParameterNonlinearValue);

      /// Destructor.
      virtual ~ElemwiseParameterNonlinear();

      /// Get the type of this instance.
      virtual ElemwiseParameterType get_type();
    protected:
      virtual Scalar get_value(Solution<Scalar>* u_ext, Element* e) = 0;

      ElemwiseParameterNonlinearValueType value_type;

      template<typename T> friend class DiscreteProblem;
      template<typename T> friend class DiscreteProblemLinear;
    };

    /// Nonlinear element-wise parameter.
    template<typename Scalar>
    class HERMES_API ElemwiseParameterNonlinearHermesFunc : public ElemwiseParameterNonlinear<Scalar>
    {
    public:
      /// Constructor.
      ElemwiseParameterNonlinearHermesFunc(Hermes1DFunction<Scalar>* function, ElemwiseParameterNonlinearValueType value_type = ElemwiseParameterNonlinearValue);
      ElemwiseParameterNonlinearHermesFunc(ElemwiseParameterNonlinearValueType value_type = ElemwiseParameterNonlinearValue);

      /// Destructor.
      virtual ~ElemwiseParameterNonlinearHermesFunc();

      /// Ask if the instance is fine.
      virtual bool isOkay() const;

        /// Get class name, for the purpose of messaging.
      virtual std::string getClassName() const;

    protected:
      virtual Scalar get_value(Solution<Scalar>* u_ext, Element* e);

      Hermes1DFunction<Scalar>* function;

      template<typename T> friend class DiscreteProblem;
      template<typename T> friend class DiscreteProblemLinear;
    };

    /// Spline-represented nonlinear parameter.
    class HERMES_API ElemwiseParameterSpline : public ElemwiseParameterNonlinearHermesFunc<double>, public Hermes2D::Mixins::StateQueryable
    {
    public:
      /// Constructor.
      ElemwiseParameterSpline(CubicSpline* spline, ElemwiseParameterNonlinearValueType value_type = ElemwiseParameterNonlinearValue);
      ElemwiseParameterSpline(ElemwiseParameterNonlinearValueType value_type = ElemwiseParameterNonlinearValue);

      /// Destructor.
      virtual ~ElemwiseParameterSpline();
    private:
      template<typename T> friend class DiscreteProblem;
      template<typename T> friend class DiscreteProblemLinear;
    };
  }
}
#endif
