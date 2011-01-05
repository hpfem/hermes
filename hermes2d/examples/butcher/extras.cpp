// This file is part of Hermes2D
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, see <http://www.gnu.prg/licenses/>.

#ifndef _BUTCHER_TABLE_
#define _BUTCHER_TABLE_

/// Butcher's table for Runge-Kutta methods.
class HERMES_API ButcherTable 
{
 public:
  ButcherTable(int size);
  int get_size() {return this->size;}; 
  double get_A(int i, int j) 
  {
    if (i < 0 || i > size || j < 0 || j > size) error("Invalid access to a Butcher's table.");
    return this->A[i][j];
  }
  double get_B(int i) 
  {
    if (i < 0 || i > size) error("Invalid access to a Butcher's table.");
    return this->B[i];
  }
  double get_C(int i) 
  {
    if (i < 0 || i > size) error("Invalid access to a Butcher's table.");
    return this->C[i];
  }
  void set_A(int i, int j, double val) 
  {
    if (i < 0 || i > size || j < 0 || j > size) error("Invalid access to a Butcher's table.");
    this->A[i][j] = val;
  }
  void set_B(int i, double val) 
  {
    if (i < 0 || i > size) error("Invalid access to a Butcher's table.");
    this->B[i] = val;
  }
  void set_C(int i, double val) 
  {
    if (i < 0 || i > size) error("Invalid access to a Butcher's table.");
    this->C[i] = val;
  }

 protected:
  int size;
  double** A;
  double* B;
  double* C;
};

ButcherTable::ButcherTable(int size) 
{
  this->size = size;
  // A array.
  this->A = new double*[size];
  for (int i=0; i<size; i++) {
    this->A[i] = new double[size];
    for (int j=0; j<size; j++) this->A[i][j] = 0;
  }
  // B array.
  this->B = new double[size];
  for (int j=0; j<size; j++) this->B[j] = 0;
  // C array.
  this->C = new double[size];
  for (int j=0; j<size; j++) this->C[j] = 0;
}


#endif
