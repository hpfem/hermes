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

  // Initializing coefficient array.
  int nelem = points.size() - 1;
  coeffs = new SplineCoeff[nelem];
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

bool CubicSpline::calculate_coeffs() 
{
  int nelem = points.size() - 1;
  const int n = 4 * nelem;
  double** matrix = new_matrix<double>(n, n);
  memset(matrix, 0, n*n*sizeof(double));
  double* rhs = new double[n];
  memset(rhs, 0, n*sizeof(double));

  // Fill the rhs vector.
  for (int i=0; i < nelem; i++) {
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
  if (first_der_left == false) { 
    matrix[offset + 0][2] = 2;
    matrix[offset + 0][3] = 6*xx;
    rhs[n-2] = bc_left; // Value of the second derivative.
  } 
  else {
    matrix[offset + 0][1] = 1;
    matrix[offset + 0][2] = 2*xx;
    matrix[offset + 0][3] = 3*xx*xx;
    rhs[n-2] = bc_left; // Value of the first derivative.
  }
  xx = points[n-1]; // Right end-point.
  if (first_der_right == false) { 
    matrix[offset + 1][n-2] = 2;
    matrix[offset + 1][n-1] = 6*xx;
    rhs[n-1] = bc_right; // Value of the second derivative.
  }
  else { 
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
  for (int i = 0; i < nelem; i++) {
    coeffs[i].a = rhs[4*i + 0];
    coeffs[i].b = rhs[4*i + 1];
    coeffs[i].c = rhs[4*i + 2];
    coeffs[i].d = rhs[4*i + 3];
  }

  // Free the matrix and rhs vector.
  delete [] matrix;
  delete [] rhs;

  return true;
}
