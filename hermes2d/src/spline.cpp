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

CubicSpline::CubicSpline(std::vector<double> points, std::vector<double> values, 
                         double bc_left, double bc_right, 
                         bool first_der_left, bool first_der_right) 
  : points(points), values(values), bc_left(bc_left), bc_right(bc_right), 
    first_der_left(first_der_left), first_der_right(first_der_right)
{
  // Sanity check.
  if (points.empty() || values.empty()) error("Supply both points and values when initializing a spline.");
  if (points.size() != values.size()) error("Mismatched number of spline points and values.");
}

bool CubicSpline::get_value(double x_in, double& val_out) 
{  
  int m = -1;
  if (!this->find_interval(x_in, m)) return false;

  val_out = get_value_from_interval(x_in, m); 

  return true;
}

double CubicSpline::get_value_from_interval(double x_in, int m) 
{  
  double x2 = x_in * x_in;
  double x3 = x2 * x_in;
  return   this->coeffs[m].a + this->coeffs[m].b * x_in + this->coeffs[m].c * x2
         + this->coeffs[m].d * x3;
}

bool CubicSpline::get_derivative(double x_in, double& der_out) 
{
  int m = -1;
  if (!this->find_interval(x_in, m)) return false;

  der_out = get_derivative_from_interval(x_in, m); 

  return true;
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

  while (i_left + 1 < i_right) {
    int i_mid = (i_left + i_right) / 2;
    if (points[i_mid] < x_in) i_left = i_mid;   
    else i_right = i_mid;  
  } 

  m = i_left;
  return true;
};

void CubicSpline::plot(const char* filename, int subdiv) 
{
  FILE *f = fopen(filename, "wb");
  if (f == NULL) error("Could not open a spline file for writing."); 
  
  for (unsigned int i = 0; i < points.size() - 1; i++) {
    double h = (points[i+1] - points[i]) / subdiv;
    for (int j = 0; j < subdiv; j++) {
      double x = points[i] + j * h;
      double val = get_value_from_interval(x, i); 
      fprintf(f, "%g %g\n", x, val);
    }
  }
  double x_last = points[points.size() - 1];
  double val_last = get_value_from_interval(x_last, points.size() - 1);
  fprintf(f, "%g %g\n", x_last, val_last);
  fclose(f);
}
