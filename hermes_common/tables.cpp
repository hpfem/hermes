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

#include "tables.h"

Table::Table(int size) 
{
  // Size.
  this->size = size;
  // A array.
  this->A = new double*[size];
  for (int i=0; i<size; i++) {
    this->A[i] = new double[size];
    for (int j=0; j<size; j++) this->A[i][j] = 0;
  }
}

int Table::get_size() 
{
  return this->size;
}

double Table::get_A(int i, int j) 
{
  if (i < 0 || i > size || j < 0 || j > size) error("Invalid access to a Butcher's table.");
  return this->A[i][j];
}

void Table::set_A(int i, int j, double val) 
{
  if (i < 0 || i > size || j < 0 || j > size) error("Invalid access to a Butcher's table.");
  this->A[i][j] = val;
}

ButcherTable::ButcherTable(int size) : Table(size)
{
  // B array.
  this->B = new double[size];
  for (int j=0; j<size; j++) this->B[j] = 0;
  // C array.
  this->C = new double[size];
  for (int j=0; j<size; j++) this->C[j] = 0;

  this->time_step = -1;
}

double ButcherTable::get_B(int i) 
{
  if (i < 0 || i > size) error("Invalid access to a Butcher's table.");
  return this->B[i];
}

double ButcherTable::get_C(int i) 
{
  if (i < 0 || i > size) error("Invalid access to a Butcher's table.");
  return this->C[i];
}

void ButcherTable::set_B(int i, double val) 
{
  if (i < 0 || i > size) error("Invalid access to a Butcher's table.");
  this->B[i] = val;
}

void ButcherTable::set_C(int i, double val) 
{
  if (i < 0 || i > size) error("Invalid access to a Butcher's table.");
  this->C[i] = val;
}

void ButcherTable::set_time_step(double tau) 
{
  this->time_step = tau;
}

double ButcherTable::get_time_step() 
{
  return this->time_step;
}

