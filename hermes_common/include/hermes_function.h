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

#include "util/compat.h"
#include "ord.h"
#include "exceptions.h"
#include "mixins.h"

namespace Hermes
{
  struct HERMES_API SplineCoeff
  {
    double a, b, c, d;    // four coefficients of a cubic spline.
  };

  /// Generic class for functions of one variable.
  /// Typical usage: user creates a descendant:
  /// class MyRealSinSquared : public Hermes1DFunction<double>
  /// {
  ///&nbsp;/// Constructor.
  ///&nbsp;MyRealSinSquared() : Hermes1DFunction<double>() {}
  ///&nbsp;double value(double x) const { return std::sin(x) * std::sin(x); }
  ///&nbsp;Ord ord(Hermes::Ord x) const { return std::sin(x) * std::sin(x); // Hermes::Ord arithmetics will take care of everything, the resulting order will be as high as possible, see other examples in this file for another possibility how to go about it. }
  ///&nbsp;....
  /// }
  template<typename Scalar>
  class HERMES_API Hermes1DFunction : public Hermes::Mixins::Loggable
  {
  public:
    /// Constructor.
    Hermes1DFunction();

    /// Constructor for the constant case.
    Hermes1DFunction(Scalar value);

    /// One-dimensional function value.
    virtual Scalar value(Scalar x) const;

    /// One-dimensional function integration order.
    virtual Hermes::Ord value(Hermes::Ord x) const;

    /// One-dimensional function derivative value.
    virtual Scalar derivative(Scalar x) const;

    /// One-dimensional function derivative integration order.
    virtual Hermes::Ord derivative(Hermes::Ord x) const;

    /// The function is constant.
    /// Returns the value of is_const.
    bool is_constant() const;

  protected:
    /// The function is constant.
    bool is_const;
    /// If the function is constant, this is the value.
    Scalar const_value;
  };

  /// Generic class for functions of two variables.
  /// Typical usage: user creates a descendant:
  /// class MyComplexASquaredMinusBSquared : public Hermes2DFunction<std::complex<double> >
  /// {
  ///&nbsp;/// Constructor.
  ///&nbsp;MyComplexASquaredMinusBSquared() : Hermes2DFunction<std::complex<double> >() {}
  ///&nbsp;std::complex<double> value(std::complex<double> x, std::complex<double> y) const { return x.real * y.real() - x.imag() * y.imag(); }
  ///&nbsp;Hermes::Ord value(Hermes::Ord x, Hermes::Ord y) const { return Hermes::Ord(10); // Just say that this particular function is of polynomial order 10. }
  ///&nbsp;....
  /// }
  template<typename Scalar>
  class HERMES_API Hermes2DFunction : public Hermes::Mixins::Loggable
  {
  public:
    /// Constructor.
    Hermes2DFunction();

    /// Constructor for the constant case.
    Hermes2DFunction(Scalar value);

    /// Two-dimensional function value.
    virtual Scalar value(Scalar x, Scalar y) const;

    /// Two-dimensional function integration order.
    virtual Hermes::Ord value(Hermes::Ord x, Hermes::Ord y) const;

    /// Two-dimensional function derivative value.
    virtual Scalar derivative_x(Scalar x, Scalar y) const;
    virtual Scalar derivative_y(Scalar x, Scalar y) const;

    /// Two-dimensional function derivative integration order.
    virtual Hermes::Ord derivative_x(Hermes::Ord x, Hermes::Ord y) const;
    virtual Hermes::Ord derivative_y(Hermes::Ord x, Hermes::Ord y) const;

    /// The function is constant.
    /// Returns the value of is_const.
    bool is_constant() const;

    /// The function is constant.
    bool is_const;
    /// If the function is constant, this is the value.
    Scalar const_value;
  };

  /// Generic class for functions of three variables.
  /// Typical usage: user creates a descendant:
  /// class MyReal3DSum : public Hermes3DFunction<double>
  /// {
  ///&nbsp;/// Constructor.
  ///&nbsp;MyReal3DSum() : Hermes3DFunction<double>() {}
  ///&nbsp;virtual double value(double x, double y, double z) const { return x + y + z; }
  ///&nbsp;....
  /// }
  template<typename Scalar>
  class HERMES_API Hermes3DFunction : public Hermes::Mixins::Loggable
  {
  public:
    /// Constructor.
    Hermes3DFunction();

    /// Constructor for the constant case.
    Hermes3DFunction(Scalar value);

    /// Two-dimensional function value.
    virtual Scalar value(Scalar x, Scalar y, Scalar z) const;

    /// Two-dimensional function integration order.
    virtual Hermes::Ord value(Hermes::Ord x, Hermes::Ord y, Hermes::Ord z) const;

    /// Two-dimensional function derivative value.
    virtual Scalar derivative_x(Scalar x, Scalar y, Scalar z) const;
    virtual Scalar derivative_y(Scalar x, Scalar y, Scalar z) const;
    virtual Scalar derivative_z(Scalar x, Scalar y, Scalar z) const;

    /// Two-dimensional function derivative integration order.
    virtual Hermes::Ord derivative_x(Hermes::Ord x, Hermes::Ord y, Hermes::Ord z) const;
    virtual Hermes::Ord derivative_y(Hermes::Ord x, Hermes::Ord y, Hermes::Ord z) const;
    virtual Hermes::Ord derivative_z(Hermes::Ord x, Hermes::Ord y, Hermes::Ord z) const;

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

#endif