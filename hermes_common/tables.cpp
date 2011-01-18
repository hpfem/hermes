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
#include "matrix.h"

Table::Table() 
{
  this->size = -1;
  this->A = NULL;
}

Table::Table(int size) 
{
  // Size.
  this->size = size;
  // A array.
  this->A = new_matrix<double>(size, size);
  for (int i=0; i<size; i++) {
    for (int j=0; j<size; j++) this->A[i][j] = 0;
  }
}

void Table::alloc(int size) 
{
  // Size.
  this->size = size;
  // A array.
  this->A = new_matrix<double>(size, size);
  for (int i=0; i<size; i++) {
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

ButcherTable::ButcherTable() : Table() 
{
  this->B = NULL;
  this->C = NULL;
}

ButcherTable::ButcherTable(int size) : Table(size)
{
  // B array.
  this->B = new double[size];
  for (int j=0; j<size; j++) this->B[j] = 0;
  // C array.
  this->C = new double[size];
  for (int j=0; j<size; j++) this->C[j] = 0;
}



ButcherTable::ButcherTable(ButcherTableType butcher_table)
{
  double gamma = 1./sqrt(2.);

  switch (butcher_table) {
    case Explicit_RK_1: // Explicit Euler.
      this->alloc(1);
      this->set_B(0, 1.);
    break;

    case Implicit_RK_1: // Implicit Euler.
      this->alloc(1);
      this->set_A(0, 0, 1.);
      this->set_B(0, 1.);
      this->set_C(0, 1.);
    break;

    case Explicit_RK_2: // Explicit RK-2.
      this->alloc(2);
      this->set_A(1, 0, 2./3.);
      this->set_A(1, 1, 0.);
      this->set_B(0, 1./4.);
      this->set_B(1, 3./4.);
      this->set_C(0, 2./3.);
    break;

    case Implicit_Crank_Nicolson_2: // Implicit Crank Nicolson.
      this->alloc(2);
      this->set_A(0, 0, 1./2.);
      this->set_A(0, 1, 1./2.);
      this->set_B(0, 1./2.);
      this->set_B(1, 1./2.);
      this->set_C(0, 1.);
    break;

    case Implicit_SDIRK_2: // Implicit SDIRK-2 (second-order).
      this->alloc(2);
      this->set_A(0, 0, 1. - gamma);
      this->set_A(0, 1, 0.);
      this->set_A(1, 0, gamma);
      this->set_A(1, 1, 1. - gamma);
      this->set_B(0, gamma);
      this->set_B(1, 1. - gamma);
      this->set_C(0, 1. - gamma);
      this->set_C(1, 1.);  
    break;

    case Implicit_Lobatto_IIIA_2: // Implicit Lobatto IIIA (second-order).
      this->alloc(2);
      this->set_A(1, 0, 1./2.);
      this->set_A(1, 1, 1./2.);
      this->set_B(0, 1./2.);
      this->set_B(1, 1./2.);
      this->set_C(1, 1.);
    break;

    case Implicit_Lobatto_IIIB_2: // Implicit Lobatto IIIB (second-order).
      this->alloc(2);
      this->set_A(0, 0, 1./2.);
      this->set_A(0, 1, 1./2.);
      this->set_B(0, 1./2.);
      this->set_B(1, 1./2.);
      this->set_C(0, 1./2.);
      this->set_C(1, 1./2.);
    break;

    case Implicit_Lobatto_IIIC_2: // Implicit Lobatto IIIC (second-order).
      this->alloc(2);
      this->set_A(0, 0, 1./2.);
      this->set_A(0, 1, -1./2.);
      this->set_A(1, 0, 1./2.);
      this->set_A(1, 1, 1./2.);
      this->set_B(0, 1./2.);
      this->set_B(1, 1./2.);
      this->set_C(0, 0.);
      this->set_C(1, 1.);
    break;

    case Explicit_RK_3: // Explicit RK-3.
      this->alloc(3);
      this->set_A(1, 0, 1./2.);
      this->set_A(2, 0, -1.);
      this->set_A(2, 1, 2.);
      this->set_B(0, 1./6.);
      this->set_B(1, 2./3.);
      this->set_B(2, 1./6.);
      this->set_C(1, 1./2.);
      this->set_C(2, 1.);
    break;

    case Explicit_RK_4: // Explicit RK-4.
      this->alloc(4);
      this->set_A(1, 0, 1./2.);
      this->set_A(2, 1, 1./2.);
      this->set_A(3, 2, 1.);
      this->set_B(0, 1./6.);
      this->set_B(1, 1./3.);
      this->set_B(2, 1./3.);
      this->set_B(3, 1./6.);
      this->set_C(1, 1./2.);
      this->set_C(2, 1./2.);
      this->set_C(3, 1.);
    break;

    case Implicit_Lobatto_IIIA_4: // Implicit Lobatto IIIA (fourth-order).
      this->alloc(3);
      this->set_A(1, 0, 5./24.);
      this->set_A(2, 0, 1./6.);
      this->set_A(1, 1, 1./3.);
      this->set_A(2, 1, 2./3.);
      this->set_A(1, 2, -1./24.);
      this->set_A(2, 2, 1./6.);
      this->set_B(0, 1./6.);
      this->set_B(1, 2./3.);
      this->set_B(2, 1./6.);
      this->set_C(1, 1./2.);
      this->set_C(2, 1.);
    break;

    case Implicit_Lobatto_IIIB_4: // Implicit Lobatto IIIB (fourth-order).
      this->alloc(3);
      this->set_A(0, 0, 1./6.);
      this->set_A(1, 0, 1./6.);
      this->set_A(2, 0, 1./6.);
      this->set_A(0, 1, -1./6.);
      this->set_A(1, 1, 1./3.);
      this->set_A(2, 1, 5./6.);
      this->set_B(0, 1./6.);
      this->set_B(1, 2./3.);
      this->set_B(2, 1./6.);
      this->set_C(1, 1./2.);
      this->set_C(2, 1.);
    break;

    case Implicit_Lobatto_IIIC_4: // Implicit Lobatto IIIC (fourth-order).
      this->alloc(3);
      this->set_A(0, 0, 1./6.);
      this->set_A(1, 0, 1./6.);
      this->set_A(2, 0, 1./6.);
      this->set_A(0, 1, -1./3.);
      this->set_A(1, 1, 5./12.);
      this->set_A(2, 1, 2./3.);
      this->set_A(0, 2, 1./6.);
      this->set_A(1, 2, -1./12.);
      this->set_A(2, 2, 1./6.);
      this->set_B(0, 1./6.);
      this->set_B(1, 2./3.);
      this->set_B(2, 1./6.);
      this->set_C(1, 1./2.);
      this->set_C(2, 1.);
    break;

    default: error("Unknown Butcher's table.");
  }
}

void ButcherTable::alloc(int size) 
{
  // Size.
  this->size = size;
  // A array.
  this->A = new_matrix<double>(size, size);
  for (int i=0; i<size; i++) {
    for (int j=0; j<size; j++) this->A[i][j] = 0;
  }
  // B array.
  this->B = new double[size];
  for (int j=0; j<size; j++) this->B[j] = 0;
  // C array.
  this->C = new double[size];
  for (int j=0; j<size; j++) this->C[j] = 0;
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

