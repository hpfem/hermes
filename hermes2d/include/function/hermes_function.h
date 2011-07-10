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
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __H2D_HERMES_FUNCTION_H
#define __H2D_HERMES_FUNCTION_H

#include "mesh_function.h"
#include "../mesh/refmap.h"
#include "../forms.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// Generic class for functions of one and two variables.
    template<typename Scalar>
    class HERMES_API HermesFunction
    {
    public:
      /// Constructor.
      HermesFunction();

      /// Constructor for the constant case.
      HermesFunction(Scalar value);

      /// One-dimensional function value.
      virtual Scalar value(Scalar x) const;

      /// One-dimensional function integration order.
      virtual Hermes::Ord value_ord(Hermes::Ord x) const;

      /// Two-dimensional function value.
      virtual Scalar value(Scalar x, Scalar y) const;

      /// Two-dimensional function integration order.
      virtual Hermes::Ord value_ord(Hermes::Ord x, Hermes::Ord y) const;

      /// One-dimensional function derivative value.
      virtual Scalar derivative(Scalar x) const;

      /// One-dimensional function derivative integration order.
      virtual Hermes::Ord derivative_ord(Hermes::Ord x) const;

      /// Two-dimensional function derivative value.
      virtual Scalar derivative(Scalar x, Scalar y) const;

      /// Two-dimensional function derivative integration order.
      virtual Hermes::Ord derivative_ord(Hermes::Ord x, Hermes::Ord y) const;

      /// The function is constant.
      /// Returns the value of is_const.
      bool is_constant() const;

    protected:
      /// The function is constant.
      bool is_const;
      /// If the function is constant, this is the value.
      Scalar const_value;
    };
  }
}

#endif