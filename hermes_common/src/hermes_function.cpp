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

#include "hermes_function.h"
#include "hermes_logging.h"

namespace Hermes
{
  template<typename Scalar>
  bool Hermes1DFunction<Scalar>::is_constant() const
  {
    return is_const;
  };

  template<typename Scalar>
  bool Hermes2DFunction<Scalar>::is_constant() const
  {
    return is_const;
  };

  template<typename Scalar>
  bool Hermes3DFunction<Scalar>::is_constant() const
  {
    return is_const;
  };


  template<typename Scalar>
  Hermes1DFunction<Scalar>::Hermes1DFunction()
  {
    this->is_const = true;
  };

  template<typename Scalar>
  Hermes1DFunction<Scalar>::Hermes1DFunction(Scalar value)
  {
    this->is_const = true;
    this->const_value = value;
  };

  template<>
  double Hermes1DFunction<double>::value(double x) const
  {
    if(this->is_const)
      return const_value;
    else
    {
      error("Abstract function call.");
      return 0.0;
    }
  };
  template<>
  std::complex<double> Hermes1DFunction<std::complex<double> >::value(std::complex<double> x) const
  {
    if(this->is_const)
      return const_value;
    else
    {
      error("Abstract function call.");
      return std::complex<double>(0.0, 0.0);
    }
  };

  template<typename Scalar>
  Ord Hermes1DFunction<Scalar>::value_ord(Ord x) const
  {
    if(this->is_const)
      return Ord(0);
    else
    {
      error("Abstract function call.");
      return Ord(99);
    }
  };

  template<>
  double Hermes1DFunction<double>::derivative(double x) const
  {
    if(this->is_const)
      return 0.0;
    else
    {
      error("Abstract function call.");
      return 0.0;
    }
  };
  template<>
  std::complex<double> Hermes1DFunction<std::complex<double> >::derivative(std::complex<double> x) const
  {
    if(this->is_const)
      return std::complex<double>(0.0, 0.0);
    else
    {
      error("Abstract function call.");
      return std::complex<double>(0.0, 0.0);
    }
  };

  template<typename Scalar>
  Ord Hermes1DFunction<Scalar>::derivative_ord(Ord x) const
  {
    if(this->is_const)
      return Ord(0);
    else
    {
      error("Abstract function call.");
      return Ord(99);
    }
  };

    
  template<typename Scalar>
  Hermes2DFunction<Scalar>::Hermes2DFunction()
  {
    this->is_const = true;
  };

  template<typename Scalar>
  Hermes2DFunction<Scalar>::Hermes2DFunction(Scalar value)
  {
    this->is_const = true;
    this->const_value = value;
  };

  template<>
  double Hermes2DFunction<double>::value(double x, double y) const
  {
    if(this->is_const)
      return const_value;
    else
    {
      error("Abstract function call.");
      return 0.0;
    }
  };

  template<>
  std::complex<double> Hermes2DFunction<std::complex<double> >::value(std::complex<double> x, std::complex<double> y) const
  {
    if(this->is_const)
      return const_value;
    else
    {
      error("Abstract function call.");
      return std::complex<double>(0.0, 0.0);
    }
  };

  template<typename Scalar>
  Ord Hermes2DFunction<Scalar>::value_ord(Ord x, Ord y) const
  {
    if(this->is_const)
      return Ord(0);
    else
    {
      error("Abstract function call.");
      return Ord(99);
    }
  };

  template<>
  double Hermes2DFunction<double>::derivative(double x, double y) const
  {
    if(this->is_const)
      return 0.0;
    else
    {
      error("Abstract function call.");
      return 0.0;
    }
  };
  template<>
  std::complex<double> Hermes2DFunction<std::complex<double> >::derivative(std::complex<double> x, std::complex<double> y) const
  {
    if(this->is_const)
      return std::complex<double>(0.0, 0.0);
    else
    {
      error("Abstract function call.");
      return std::complex<double>(0.0, 0.0);
    }
  };

  template<typename Scalar>
  Ord Hermes2DFunction<Scalar>::derivative_ord(Ord x, Ord y) const
  {
    if(this->is_const)
      return Ord(0);
    else
    {
      error("Abstract function call.");
      return Ord(99);
    }
  };

  template<typename Scalar>
  Hermes3DFunction<Scalar>::Hermes3DFunction()
  {
    this->is_const = true;
  };

  template<typename Scalar>
  Hermes3DFunction<Scalar>::Hermes3DFunction(Scalar value)
  {
    this->is_const = true;
    this->const_value = value;
  };

  template<>
  double Hermes3DFunction<double>::value(double x, double y, double z) const
  {
    if(this->is_const)
      return 0.0;
    else
    {
      error("Abstract function call.");
      return 0.0;
    }
  };
  template<>
  std::complex<double> Hermes3DFunction<std::complex<double> >::value(std::complex<double> x, std::complex<double> y, std::complex<double> z) const
  {
    if(this->is_const)
      return std::complex<double>(0.0, 0.0);
    else
    {
      error("Abstract function call.");
      return std::complex<double>(0.0, 0.0);
    }
  };


  template<typename Scalar>
  Ord Hermes3DFunction<Scalar>::value_ord(Ord x, Ord y, Ord z) const
  {
    if(this->is_const)
      return Ord(0);
    else
    {
      error("Abstract function call.");
      return Ord(99);
    }
  };

  template<>
  double Hermes3DFunction<double>::derivative(double x, double y, double z) const
  {
    if(this->is_const)
      return 0.0;
    else
    {
      error("Abstract function call.");
      return 0.0;
    }
  };
  template<>
  std::complex<double> Hermes3DFunction<std::complex<double> >::derivative(std::complex<double> x, std::complex<double> y, std::complex<double> z) const
  {
    if(this->is_const)
      return std::complex<double>(0.0, 0.0);
    else
    {
      error("Abstract function call.");
      return std::complex<double>(0.0, 0.0);
    }
  };

  template<typename Scalar>
  Ord Hermes3DFunction<Scalar>::derivative_ord(Ord x, Ord y, Ord z) const
  {
    if(this->is_const)
      return Ord(0);
    else
    {
      error("Abstract function call.");
      return Ord(99);
    }
  };
    
  template class HERMES_API Hermes1DFunction<double>;
  template class HERMES_API Hermes1DFunction<std::complex<double> >;
  template class HERMES_API Hermes2DFunction<double>;
  template class HERMES_API Hermes2DFunction<std::complex<double> >;
}