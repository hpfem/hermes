// This file is part of Hermes
//
// Hermes is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// Hermes is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes; if not, see <http://www.gnu.prg/licenses/>.

#ifndef _TABLES_H_
#define _TABLES_H_

#include "common.h"

// General square table of real numbers.
class HERMES_API Table 
{
 public:
  Table();
  Table(int size);
  virtual void alloc(int size); 
  int get_size(); 
  double get_A(int i, int j);
  void set_A(int i, int j, double val);

 protected:
  int size;
  double** A;
};

/// Butcher's tables for Runge-Kutta methods.
class HERMES_API ButcherTable: public Table
{
 public:
 ButcherTable();
 ButcherTable(int size);
  virtual void alloc(int size);
  double get_B(int i);
  double get_C(int i);
  void set_B(int i, double val);
  void set_C(int i, double val); 

 protected:
  double* B;
  double* C;
};

#endif
