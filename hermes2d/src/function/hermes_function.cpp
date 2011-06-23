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
#include "refmap.h"
#include "hermes_function.h"
namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    HermesFunction<Scalar>::HermesFunction()
    {
      this->is_const = true;
      this->const_value = -9999.0;
    };

    template<typename Scalar>
    HermesFunction<Scalar>::HermesFunction(Scalar value)
    {
      this->is_const = true;
      this->const_value = value;
    };

    template<typename Scalar>
    Scalar HermesFunction<Scalar>::value(Scalar x) const
    {
      return const_value;
    };

    template<typename Scalar>
    Ord HermesFunction<Scalar>::value_ord(Ord x) const
    {
      return Ord(1);
    };

    template<typename Scalar>
    Scalar HermesFunction<Scalar>::value(Scalar x, Scalar y) const
    {
      return const_value;
    };

    template<typename Scalar>
    Ord HermesFunction<Scalar>::value_ord(Ord x, Ord y) const
    {
      return Ord(1);
    };

    template<>
    double HermesFunction<double>::derivative(double x) const
    {
      return 0.0;
    };
    template<>
    std::complex<double> HermesFunction<std::complex<double> >::derivative(std::complex<double> x) const
    {
      return std::complex<double>(0.0, 0.0);
    };

    template<typename Scalar>
    Ord HermesFunction<Scalar>::derivative_ord(Ord x) const
    {
      return Ord(1);
    };

    template<>
    double HermesFunction<double>::derivative(double x, double y) const
    {
      return 0.0;
    };
    template<>
    std::complex<double> HermesFunction<std::complex<double> >::derivative(std::complex<double> x, std::complex<double> y) const
    {
      return 0.0;
    };

    template<typename Scalar>
    Ord HermesFunction<Scalar>::derivative_ord(Ord x, Ord y) const
    {
      return Ord(1);
    };

    template<typename Scalar>
    bool HermesFunction<Scalar>::is_constant() const
    {
      return is_const;
    };

    template class HERMES_API HermesFunction<double>;
    template class HERMES_API HermesFunction<std::complex<double> >;
  }
}