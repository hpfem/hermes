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

#include "hermes2d.h"

#ifndef __H2D_SPLINE_H
#define __H2D_SPLINE_H
namespace Hermes
{
  namespace Hermes2D
  {
    class HERMES_API CubicSpline
    {
    public:
      /// Constructor (general case).
      CubicSpline(Hermes::vector<double> points, Hermes::vector<double> values, 
        double bc_left, double bc_right, 
        bool first_der_left = true, bool first_der_right = true,
        bool extend_der_left = true, bool extend_der_right = true);
      /// Constructor (trivial constant case).
      CubicSpline(double const_value) : is_constant(true), const_value(const_value) { };

      /// Destructor.
      ~CubicSpline() { 
        if (is_constant == false) {
          if (coeffs != NULL) delete [] coeffs;
          points.clear();
          values.clear();
        }
      };

      /// Calculates coefficients.
      bool calculate_coeffs();

      /// \brief Get the value at a given point.
      double get_value(double x_in);

      /// \brief Get the value at a given point.
      double get_value(std::complex<double> x_in) {
        return get_value(x_in.real());
      };

      /// \brief For order calculation in Hermes.
      Hermes::Ord get_value(Hermes::Ord x_in) {return Hermes::Ord(3);};
      /// \brief Get first derivative at a given point.
      double get_derivative(double x_in);
      /// \brief Get first derivative at a given point.
      double get_derivative(std::complex<double> x_in) {
        return get_derivative(x_in.real());
      };

      /// For order calculation in Hermes.
      Hermes::Ord get_derivative(Hermes::Ord x_in) {return Hermes::Ord(2);};

      /// Plots the spline in format for Pylab (just pairs 
      /// x-coordinate and value per line). The interval of definition 
      /// of the spline will be extended by "extension" both to the left 
      /// and to the right. This allows the user to see how the code will
      /// handle the spline if used for points that lie outside of its
      /// interval of definition. If plot_value == false, derivative is plotted.
      void plot(const char* filename, double extension, bool plot_derivative = false, int subdiv = 50);

    protected:
      /// Flag indicating whether the spline is constant (simplified evaluation).
      bool is_constant;

      /// Value for constant spline.
      double const_value;
      /// Uses a bisection method to locale interval where a given point lies.
      /// Returns false if point lies outside.
      bool find_interval(double x_in, int& m);

      /// Extrapolate the value of the spline outside of its interval of definition.
      double extrapolate_value(double point_end, double value_end, double derivative_end, double x_in);
      /// Grid points, ordered.
      Hermes::vector<double> points;

      /// Values at the grid points.
      Hermes::vector<double> values;

      /// Boundary conditions.
      double bc_left, bc_right;

      /// Flags that determine the meaning of the boundary constants.
      /// first_der_left == true means that the left BC is the first derivative.
      /// first_der_left == false means that the left BC is the second derivative.
      /// (same on the right)
      bool first_der_left, first_der_right;
      /// If extend_der_left == true then the spline is extended to the left of the 
      /// interval of definition as a linear function whose slope matches the derivative 
      /// at the left-most point. Otherwise the spline is extended as a constant value
      /// that matches the spline value at the left-most point. Analogously for the 
      /// flag extend_der_right.
      bool extrapolate_der_left, extrapolate_der_right;

      /// Values and derivatives at end points for extrapolation purposes.
      double point_left, value_left, derivative_left;
      double point_right, value_right, derivative_right;

      /// A set of four coefficients a, b, c, d for an elementary cubic spline.
      SplineCoeff* coeffs;

      /// Gets derivative at a point that lies in interval 'm'.
      double get_derivative_from_interval(double x_in, int m);

      /// Gets value at a point that lies in interval 'm'.
      double get_value_from_interval(double x_in, int m);
    };
  }
}
#endif
