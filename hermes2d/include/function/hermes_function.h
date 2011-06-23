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

#include "mesh_function.h"
#include "../mesh/refmap.h"
#include "../form/forms.h"
namespace Hermes
{
  namespace Hermes2D
  {
    // Generic class for functions of one and two variables.
    template<typename Scalar>
    class HERMES_API HermesFunction
    {
    public:
      HermesFunction();

      HermesFunction(Scalar value);

      virtual Scalar value(Scalar x) const;

      virtual Hermes::Ord value_ord(Hermes::Ord x) const;

      virtual Scalar value(Scalar x, Scalar y) const;

      virtual Hermes::Ord value_ord(Hermes::Ord x, Hermes::Ord y) const;

      virtual Scalar derivative(Scalar x) const;

      virtual Hermes::Ord derivative_ord(Hermes::Ord x) const;

      virtual Scalar derivative(Scalar x, Scalar y) const;

      virtual Hermes::Ord derivative_ord(Hermes::Ord x, Hermes::Ord y) const;

      bool is_constant() const;

    protected:
      bool is_const;
      Scalar const_value;
    };
  }
}