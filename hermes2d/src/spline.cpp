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

#include <typeinfo>
#include "spline.h"

using namespace Hermes::Algebra::DenseMatrixOperations;
namespace Hermes
{
  namespace Hermes2D
  {
    CubicSpline::CubicSpline(Hermes::vector<double> points, Hermes::vector<double> values, 
      double bc_left, double bc_right, 
      bool first_der_left, bool first_der_right,
      bool extrapolate_der_left, bool extrapolate_der_right) 
      : is_constant(false), const_value(-9999), points(points), values(values), 
      bc_left(bc_left), bc_right(bc_right), first_der_left(first_der_left), 
      first_der_right(first_der_right), extrapolate_der_left(extrapolate_der_left), 
      extrapolate_der_right(extrapolate_der_right) 
    {
      bool success = this->calculate_coeffs(); 
      if (!success) error("There was a problem constructing a cubic spline.");
    }

    double CubicSpline::get_value(double x_in) 
    {
      // For constant case.
      if (this->is_constant) 
      {
        return this->const_value;
      }

      // For general case.
      int m = -1;
      if (!this->find_interval(x_in, m)) 
      {
        // Point lies on the left of interval of definition.
        if (x_in <= point_left) 
        {
          // Spline should be extrapolated by constant function 
          // matching the value at the end.
          if (extrapolate_der_left == false) return value_left;
          // Spline should be extrapolated as a linear function 
          // matching the derivative at the end.
          else return extrapolate_value(point_left, value_left, derivative_left, x_in);
        }
        // Point lies on the right of interval of definition.
        else 
        {
          // Spline should be extrapolated by constant function 
          // matching the value at the end.
          if (extrapolate_der_right == false) return value_right;
          // Spline should be extrapolated as a linear function 
          // matching the derivative at the end.
          else return extrapolate_value(point_right, value_right, derivative_right, x_in);
        }
      }

      return get_value_from_interval(x_in, m); 
    }

    double CubicSpline::extrapolate_value(double point_end, double value_end, 
      double derivative_end, double x_in)
    {
      return value_end + derivative_end * (x_in - point_end);
    }

    double CubicSpline::get_value_from_interval(double x_in, int m) 
    {
      double x2 = x_in * x_in;
      double x3 = x2 * x_in;
      return   this->coeffs[m].a + this->coeffs[m].b * x_in + this->coeffs[m].c * x2
        + this->coeffs[m].d * x3;
    }

    double CubicSpline::get_derivative(double x_in) 
    {
      // For constant case.
      if (this->is_constant) return 0.0;

      // For general case.
      int m = -1;
      if (!this->find_interval(x_in, m)) 
      {
        // Point lies on the left of interval of definition.
        if (x_in <= point_left) 
        {
          // Spline should be extrapolated by constant function 
          // matching the value at the end.
          if (extrapolate_der_left == false) return 0;
          // Spline should be extrapolated as a linear function 
          // matching the derivative at the end.
          else return derivative_left;
        }
        // Point lies on the right of interval of definition.
        else 
        {
          // Spline should be extrapolated by constant function 
          // matching the value at the end.
          if (extrapolate_der_right == false) return 0;
          // Spline should be extrapolated as a linear function 
          // matching the derivative at the end.
          else return derivative_right;
        }
      }

      return get_derivative_from_interval(x_in, m); 
    };

    double CubicSpline::get_derivative_from_interval(double x_in, int m) 
    {
      double x2 = x_in * x_in;
      return this->coeffs[m].b + 2 * this->coeffs[m].c * x_in
        + 3 * this->coeffs[m].d * x2;
    }

    bool CubicSpline::find_interval(double x_in, int &m) 
    {
      int i_left = 0;
      int i_right = points.size() - 1;

      if (x_in < points[i_left]) return false;
      if (x_in > points[i_right]) return false;

      while (i_left + 1 < i_right) 
      {
        int i_mid = (i_left + i_right) / 2;
        if (points[i_mid] < x_in) i_left = i_mid;   
        else i_right = i_mid;  
      } 

      m = i_left;
      return true;
    };

    void CubicSpline::plot(const char* filename, double extension, bool plot_derivative, int subdiv) 
    {
      FILE *f = fopen(filename, "wb");
      if (f == NULL) error("Could not open a spline file for writing."); 

      // Plotting on the left of the area of definition.
      double x_left = point_left - extension;
      double h = extension / subdiv;
      for (int j = 0; j < subdiv; j++) 
      {
        double x = x_left + j * h;
        double val;
        if (!plot_derivative) val = get_value(x);
        else val = get_derivative(x); 
        fprintf(f, "%g %g\n", x, val);
      }
      double x_last = point_left;
      double val_last;
      if (!plot_derivative) val_last = get_value(x_last);
      else val_last = get_derivative(x_last); 
      fprintf(f, "%g %g\n", x_last, val_last);

      // Plotting inside the interval of definition.
      for (unsigned int i = 0; i < points.size() - 1; i++) 
      {
        double h = (points[i+1] - points[i]) / subdiv;
        for (int j = 0; j < subdiv; j++) 
        {
          double x = points[i] + j * h;
          double val;
          // Do not use get_value_from_interval() here.
          if (!plot_derivative) val = this->get_value(x); 
          else val = get_derivative(x); 
          fprintf(f, "%g %g\n", x, val);
        }
      }
      x_last = points[points.size() - 1];
      // Do not use get_value_from_interval() here.
      if (!plot_derivative) val_last = get_value(x_last); 
      else val_last = get_derivative(x_last);
      fprintf(f, "%g %g\n", x_last, val_last);

      // Plotting on the right of the area of definition.
      double x_right = point_right + extension;
      h = extension / subdiv;
      for (int j = 0; j < subdiv; j++) 
      {
        double x = point_right + j * h;
        double val;
        if (!plot_derivative) val = get_value(x); 
        else val = get_derivative(x); 
        fprintf(f, "%g %g\n", x, val);
      }
      x_last = x_right;
      if (!plot_derivative) val_last = get_value(x_last);
      else val_last = get_derivative(x_last);
      fprintf(f, "%g %g\n", x_last, val_last);

      fclose(f);
    }

    bool CubicSpline::calculate_coeffs() 
    {
      int nelem = points.size() - 1;

      // Basic sanity checks.
      if (points.empty() || values.empty()) 
      {
        warn("Empty points or values vector in CubicSpline, returning false.");
        return false;
      }
      if (points.size() < 2 || values.size() < 2) 
      {
        warn("At least two points and values required in CubicSpline, returning false.");
        return false;
      }
      if (points.size() != values.size()) 
      {
        warn("Mismatched number fo points and values in CubicSpline, returning false.");
        return false;
      }

      // Check for improperly ordered or duplicated points.
      double eps = 1e-8;
      for (int i = 0; i < nelem; i++) 
      {
        if (points[i+1] < points[i] + eps) 
        {
          warn("Duplicated or improperly ordered points in CubicSpline detected, returning false.");
          return false;
        }
      }

      /* START COMPUTATION */

      // Initializing coefficient array.
      coeffs = new SplineCoeff[nelem];

      // Allocate matrix and rhs.
      const int n = 4 * nelem;
      double** matrix = new_matrix<double>(n, n);
      for (int i = 0; i < n; i++) 
      {
        for (int j = 0; j < n; j++) 
        {
          matrix[i][j] = 0;
        }
      } 
      double* rhs = new double[n];
      for (int j = 0; j < n; j++) 
      {
        rhs[j] = 0;
      }

      // Fill the rhs vector.
      for (int i=0; i < nelem; i++) 
      {
        rhs[2*i] = values[i];
        rhs[2*i+1] = values[i+1];
      }

      // Fill the matrix. Step 1 - match values at interval endpoints.
      // This will generate the first 2*nelem rows.
      for (int i = 0; i < nelem; i++) { // Loop over elements.
        double xx = points[i];
        double xx2 = xx*xx;
        double xx3 = xx2 * xx;
        matrix[2*i][4*i + 0] = 1;
        matrix[2*i][4*i + 1] = xx;
        matrix[2*i][4*i + 2] = xx2;
        matrix[2*i][4*i + 3] = xx3;
        xx = points[i+1];
        xx2 = xx*xx;
        xx3 = xx2 * xx;
        matrix[2*i + 1][4*i + 0] = 1.0;
        matrix[2*i + 1][4*i + 1] = xx;
        matrix[2*i + 1][4*i + 2] = xx2;
        matrix[2*i + 1][4*i + 3] = xx3;
      }

      // Step 2: match first derivatives at all interior points. 
      // This will generate additional n_elem-1 rows in the matrix.
      int offset = 2*nelem;
      for (int i = 1; i < nelem; i++) { // Loop over internal points.
        double xx = points[i];
        double xx2 = xx*xx;
        matrix[offset + i-1][4*(i-1) + 1] = 1;
        matrix[offset + i-1][4*(i-1) + 2] = 2*xx;
        matrix[offset + i-1][4*(i-1) + 3] = 3*xx2;
        matrix[offset + i-1][4*(i-1) + 5] = -1;
        matrix[offset + i-1][4*(i-1) + 6] = -2*xx;
        matrix[offset + i-1][4*(i-1) + 7] = -3*xx2;
      }

      // Step 3: match second derivatives at all interior points.
      // This will generate additional n_elem-1 rows in the matrix.
      offset = 2*nelem + nelem - 1;
      for (int i = 1; i < nelem; i++) { // Loop over internal points.
        double xx = points[i];
        matrix[offset + i-1][4*(i-1) + 2] = 2;
        matrix[offset + i-1][4*(i-1) + 3] = 6*xx;
        matrix[offset + i-1][4*(i-1) + 6] = -2;
        matrix[offset + i-1][4*(i-1) + 7] = -6*xx;
      }

      // Step 4: Additional two conditions are needed to define 
      // a cubic spline. This will generate the last two rows in 
      // the matrix. Setting the second derivative (curvature) at both
      // endpoints equal to zero will result into "natural cubic spline",
      // but you can also prescribe non-zero values, or decide to 
      // prescribe first derivative (slope).  
      // Choose just one of the following two variables to be True,
      // and state the corresponding value for the derivative.
      offset = 2*nelem + 2 * (nelem - 1);
      double xx = points[0]; // Left end-point.
      if (first_der_left == false) 
      {
        matrix[offset + 0][2] = 2;
        matrix[offset + 0][3] = 6*xx;
        rhs[n-2] = bc_left; // Value of the second derivative.
      } 
      else 
      {
        matrix[offset + 0][1] = 1;
        matrix[offset + 0][2] = 2*xx;
        matrix[offset + 0][3] = 3*xx*xx;
        rhs[n-2] = bc_left; // Value of the first derivative.
      }
      xx = points[nelem]; // Right end-point.
      if (first_der_right == false) 
      {
        matrix[offset + 1][n-2] = 2;
        matrix[offset + 1][n-1] = 6*xx;
        rhs[n-1] = bc_right; // Value of the second derivative.
      }
      else 
      {
        matrix[offset + 1][n-3] = 1;
        matrix[offset + 1][n-2] = 2*xx;
        matrix[offset + 1][n-1] = 3*xx*xx;
        rhs[n-1] = bc_right; // Value of the first derivative.
      }

      // Solve the matrix problem.
      double d;
      int* perm = new int[n];
      ludcmp(matrix, n, perm, &d);
      lubksb<double>(matrix, n, perm, rhs);

      // Copy the solution into the coeffs array.
      for (int i = 0; i < nelem; i++) 
      {
        coeffs[i].a = rhs[4*i + 0];
        coeffs[i].b = rhs[4*i + 1];
        coeffs[i].c = rhs[4*i + 2];
        coeffs[i].d = rhs[4*i + 3];
      }

      // Define end point values and derivatives so that 
      // the points[] and values[] arrays are no longer
      // needed.
      point_left = points[0];
      value_left = values[0];
      derivative_left = get_derivative_from_interval(point_left, 0);
      point_right = points[points.size() - 1];
      value_right = values[values.size() - 1];
      derivative_right = get_derivative_from_interval(point_right, points.size() - 2);

      // Free the matrix and rhs vector.
      delete [] matrix;
      delete [] rhs;

      return true;
    }
  }
}