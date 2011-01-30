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
  this->size = 0;
  this->A = NULL;
}

Table::Table(unsigned int size) 
{
  // Size.
  this->size = size;
  // A array.
  this->A = new_matrix<double>(size, size);
  for (unsigned int i=0; i<size; i++) {
    for (unsigned int j=0; j<size; j++) this->A[i][j] = 0;
  }
}

void Table::alloc(unsigned int size) 
{
  // Size.
  this->size = size;
  // A array.
  this->A = new_matrix<double>(size, size);
  for (unsigned int i=0; i<size; i++) {
    for (unsigned int j=0; j<size; j++) this->A[i][j] = 0;
  }
}

unsigned int Table::get_size() 
{
  return this->size;
}

double Table::get_A(unsigned int i, unsigned int j) 
{
  if (i > size || j > size) error("Invalid access to a Butcher's table.");
  return this->A[i][j];
}

void Table::set_A(unsigned int i, unsigned int j, double val) 
{
  if (i > size || j > size) error("Invalid access to a Butcher's table.");
  this->A[i][j] = val;
}

ButcherTable::ButcherTable() : Table() 
{
  this->B = NULL;
  this->B2 = NULL;
  this->C = NULL;
}

ButcherTable::ButcherTable(unsigned int size) : Table(size)
{
  // B array.
  this->B = new double[size];
  for (unsigned int j=0; j<size; j++) this->B[j] = 0;
  // B2 array.
  this->B2 = new double[size];
  for (unsigned int j=0; j<size; j++) this->B2[j] = 0;
  // C array.
  this->C = new double[size];
  for (unsigned int j=0; j<size; j++) this->C[j] = 0;
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

    case Implicit_Crank_Nicolson_2_2: // Implicit Crank Nicolson.
      this->alloc(2);
      this->set_A(0, 0, 1./2.);
      this->set_A(0, 1, 1./2.);
      this->set_B(0, 1./2.);
      this->set_B(1, 1./2.);
      this->set_C(0, 1.);
    break;

    case Implicit_SIRK_2_2: // Implicit SIRK-22 (second-order).
      this->alloc(2);
      this->set_A(0, 0, (5. - 3.*sqrt(2.))/4.);
      this->set_A(0, 1, (7. - 5.*sqrt(2.))/4.);
      this->set_A(1, 0, (1. + 1.*sqrt(2.))/4.);
      this->set_A(1, 1, (3. - 1.*sqrt(2.))/4.);
      this->set_B(0, (1. + 1.*sqrt(2.))/4.);
      this->set_B(1, (3. - 1.*sqrt(2.))/4.);
      this->set_C(0, 3. - 2.*sqrt(2.));
      this->set_C(1, 1.);  
    break;

    case Implicit_ESIRK_2_2: // Implicit ESIRK-22 (second-order).
      this->alloc(2);
      this->set_A(0, 0, (9. - 6.*sqrt(2.))/4.);
      this->set_A(0, 1, (-3. + 2.*sqrt(2.))/4.);
      this->set_A(1, 0, (11. - 6.*sqrt(2.))/4.);
      this->set_A(1, 1, (-1. + 2.*sqrt(2.))/4.);
      this->set_B(0, 2. - sqrt(2.0));
      this->set_B(1, -1. + sqrt(2.0));
      this->set_C(1, 1.);  
    break;

    case Implicit_SDIRK_2_2: // Implicit SDIRK-22 (second-order).
      this->alloc(2);
      this->set_A(0, 0, 1. - gamma);
      this->set_A(1, 0, gamma);
      this->set_A(1, 1, 1. - gamma);
      this->set_B(0, gamma);
      this->set_B(1, 1. - gamma);
      this->set_C(0, 1. - gamma);
      this->set_C(1, 1.);  
    break;

    case Implicit_Lobatto_IIIA_2_2: // Implicit Lobatto IIIA (second-order).
      this->alloc(2);
      this->set_A(1, 0, 1./2.);
      this->set_A(1, 1, 1./2.);
      this->set_B(0, 1./2.);
      this->set_B(1, 1./2.);
      this->set_C(1, 1.);
    break;

    case Implicit_Lobatto_IIIB_2_2: // Implicit Lobatto IIIB (second-order).
      this->alloc(2);
      this->set_A(0, 0, 1./2.);
      this->set_A(0, 1, 1./2.);
      this->set_B(0, 1./2.);
      this->set_B(1, 1./2.);
      this->set_C(0, 1./2.);
      this->set_C(1, 1./2.);
    break;

    case Implicit_Lobatto_IIIC_2_2: // Implicit Lobatto IIIC (second-order).
      this->alloc(2);
      this->set_A(0, 0, 1./2.);
      this->set_A(0, 1, -1./2.);
      this->set_A(1, 0, 1./2.);
      this->set_A(1, 1, 1./2.);
      this->set_B(0, 1./2.);
      this->set_B(1, 1./2.);
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

    case Implicit_Lobatto_IIIA_3_4: // Implicit Lobatto IIIA (fourth-order).
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

    case Implicit_Lobatto_IIIB_3_4: // Implicit Lobatto IIIB (fourth-order).
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

    case Implicit_Lobatto_IIIC_3_4: // Implicit Lobatto IIIC (fourth-order).
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

    case Implicit_Radau_IIA_3_5: // Implicit Radau IIA (fifth-order).
      this->alloc(3);
      this->set_A(0, 0, (88. - 7*sqrt((double)6)) / 360 );
      this->set_A(1, 0, (296 + 169 * sqrt((double)6)) / 1800 );
      this->set_A(2, 0, (16. - sqrt((double)6)) / 36 );
      this->set_A(0, 1, (296 - 169 * sqrt((double)6)) / 1800);
      this->set_A(1, 1, (88. + 7*sqrt((double)6)) / 360);
      this->set_A(2, 1, (16. + sqrt((double)6)) / 36); 
      this->set_A(0, 2, (-2. + 3 * sqrt((double)6)) / 225 );
      this->set_A(1, 2, (-2. - 3 * sqrt((double)6)) / 225 );
      this->set_A(2, 2, 1./9.);
      this->set_B(0, (16. - sqrt((double)6)) / 36 );
      this->set_B(1, (16. + sqrt((double)6)) / 36 );
      this->set_B(2, 1./9.);
      this->set_C(0, (4. - sqrt((double)6)) / 10 );
      this->set_C(1, (4. + sqrt((double)6)) / 10 );
      this->set_C(2, 1.);
    break;

    case Implicit_SDIRK_4_5: // Implicit SDIRK-45 (fifth-order).
      this->alloc(5);
      this->set_A(0, 0, 1./4.);
      this->set_A(1, 0, 1./2.);
      this->set_A(1, 1, 1./4.);
      this->set_A(2, 0, 17./50.);
      this->set_A(2, 1, -1./25.);
      this->set_A(2, 2, 1./4.);
      this->set_A(3, 0, 371./1360.);
      this->set_A(3, 1, -137./2720.);
      this->set_A(3, 2, 15./544.);
      this->set_A(3, 3, 1./4.);
      this->set_A(4, 0, 25./24.);
      this->set_A(4, 1, -49./48.);
      this->set_A(4, 2, 125./16.);
      this->set_A(4, 3, -85./12.);
      this->set_A(4, 4, 1./4.);
      this->set_B(0, 25./24.);
      this->set_B(1, -49./48.);
      this->set_B(2, 125./16.);
      this->set_B(3, -85./12.);
      this->set_B(4, 1./4.);
      this->set_C(0, 1./4.);
      this->set_C(1, 3./4.);  
      this->set_C(2, 11./20.);  
      this->set_C(3, 1./2.);  
      this->set_C(4, 1.);  
    break;

    case Implicit_DIRK_7_45_embedded: // Implicit embedded DIRK method with orders 4 and 5.
      this->alloc(7);
      this->set_A(0, 0, 0);
      this->set_A(1, 0, 0.28589);
      this->set_A(2, 0, 0.142945);
      this->set_A(3, 0, 0.16803599);
      this->set_A(4, 0, 0.182315);
      this->set_A(5, 0, 0.24756392);
      this->set_A(6, 0, 0.13001804);
      this->set_A(1, 1, 0.28589);
      this->set_A(2, 1, 0.924011005);
      this->set_A(3, 1, -0.049416510);
      this->set_A(4, 1, -0.112951603);
      this->set_A(5, 1, -0.425378071);
      this->set_A(6, 1, 0);
      this->set_A(2, 2, 0.28589);
      this->set_A(3, 2, -0.004509476);
      this->set_A(4, 2, -0.0277933233);
      this->set_A(5, 2, -0.107036282);
      this->set_A(6, 2, -0.019290177);
      this->set_A(3, 3, 0.28589);
      this->set_A(4, 3, 0.422539833);
      this->set_A(5, 3, 0.395700134);
      this->set_A(6, 3, 0.535386266);
      this->set_A(4, 4, 0.28589);
      this->set_A(5, 4, 0.503260302);
      this->set_A(6, 4, 0.234313169);
      this->set_A(5, 5, 0.28589);
      this->set_A(6, 5, -0.166317293);
      this->set_A(6, 6, 0.28589);
      this->set_B(0, 0.13001804);
      this->set_B(1, 0);
      this->set_B(2, -0.019290177);
      this->set_B(3, 0.535386266);
      this->set_B(4, 0.234313169);
      this->set_B(5, -0.166317293);
      this->set_B(6, 0.28589);
      this->set_B2(0, -0.094388662);
      this->set_B2(1, 0);
      this->set_B2(2, -0.039782614);
      this->set_B2(3, 0.745608552);
      this->set_B2(4, -0.505129807);
      this->set_B2(5, 0.704915206);
      this->set_B2(6, 0.28589);
      this->set_C(0, 0);
      this->set_C(1, 0.57178);
      this->set_C(2, 1.352846);
      this->set_C(3, 0.4);
      this->set_C(4, 0.75);
      this->set_C(5, 0.9);
      this->set_C(6, 1.0);
    break;

    default: error("Unknown Butcher's table.");
  }
}

void ButcherTable::alloc(unsigned int size) 
{
  // Size.
  this->size = size;
  // A array.
  this->A = new_matrix<double>(size, size);
  for (unsigned int i=0; i<size; i++) {
    for (unsigned int j=0; j<size; j++) this->A[i][j] = 0;
  }
  // B array.
  this->B = new double[size];
  for (unsigned int j=0; j<size; j++) this->B[j] = 0;
  // B2 array.
  this->B2 = new double[size];
  for (unsigned int j=0; j<size; j++) this->B2[j] = 0;
  // C array.
  this->C = new double[size];
  for (unsigned int j=0; j<size; j++) this->C[j] = 0;
}

double ButcherTable::get_B(unsigned int i) 
{
  if (i > size) error("Invalid access to a Butcher's table.");
  return this->B[i];
}

double ButcherTable::get_B2(unsigned int i) 
{
  if (i > size) error("Invalid access to a Butcher's table.");
  return this->B2[i];
}

double ButcherTable::get_C(unsigned int i) 
{
  if (i > size) error("Invalid access to a Butcher's table.");
  return this->C[i];
}

void ButcherTable::set_B(unsigned int i, double val) 
{
  if (i > size) error("Invalid access to a Butcher's table.");
  this->B[i] = val;
}

void ButcherTable::set_B2(unsigned int i, double val) 
{
  if (i > size) error("Invalid access to a Butcher's table.");
  this->B2[i] = val;
}

void ButcherTable::set_C(unsigned int i, double val) 
{
  if (i > size) error("Invalid access to a Butcher's table.");
  this->C[i] = val;
}

bool ButcherTable::is_explicit()
{
  bool result = true;
  for (int i=0; i<size; i++) {
    for (int j=0; j<size; j++) {
      double val_ij = get_A(i, j);
      if (j >= i && fabs(val_ij) > 1e-12) result = false;
    }
  }

  return result;
} 

bool ButcherTable::is_diagonally_implicit()
{
  bool result = true;
  for (int i=0; i<size; i++) {
    for (int j=0; j<size; j++) {
      double val_ij = get_A(i, j);
      if (j > i && fabs(val_ij) > 1e-12) result = false;
    }
  }

  return result;
} 

bool ButcherTable::is_fully_implicit()
{
  bool result = false;
  for (int i=0; i<size; i++) {
    for (int j=0; j<size; j++) {
      double val_ij = get_A(i, j);
      if (j > i && fabs(val_ij) > 1e-12) result = true;
    }
  }

  return result;
} 
