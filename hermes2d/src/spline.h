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

class HERMES_API CubicSpline
{
public:
  /// Constructor.
  CubicSpline(std::vector<double> points, std::vector<double> values, 
              double bc_left, double bc_right, 
              bool first_der_left = true, bool first_der_right = true);

  /// Destructor.
  ~CubicSpline() { 
    delete coeffs;
    points.clear();
    values.clear();
  };

  /// Calculate the value at a given point. Return true if point is
  /// in range, otherwise return false.
  bool get_value(double x_in, double& val_out);

  /// Calculate first derivative at a given point. Return true if point is
  /// in range, otherwise return false.
  bool get_derivative(double x_in, double& der_out);

protected:
  /// Calculates coefficients.
  bool calculate_coeffs();

  /// Uses a bisection method to locale interval where a given point lies.
  /// Returns false if point lies outside.
  bool find_interval(double x_in, int& m);

  /// Gets value at a point that lies in interval 'm'.
  double get_value_from_interval(double x_in, int m);

  /// Gets derivative at a point that lies in interval 'm'.
  double get_derivative_from_interval(double x_in, int m);

  /// Plots the spline in format for Pylab (just pairs 
  /// x-coordinate and value per line).
  void plot(const char* filename, int subdiv = 20);

  /// Grid points, ordered.
  std::vector<double> points;

  /// Values at the grid points.
  std::vector<double> values;

  /// Boundary conditions.
  double bc_left, bc_right;

  /// Flags that determine the meaning of the boundary constants.
  /// first_der_left == true means that the left BC is the first derivative.
  /// first_der_left == false means that the left BC is the second derivative.
  /// (same on the right)
  bool first_der_left, first_der_right;

  /// A set of four coefficients a, b, c, d for an elementary cubic spline.
  SplineCoeff* coeffs;
};

#endif
