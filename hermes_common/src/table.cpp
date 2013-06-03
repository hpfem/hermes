// This file is part of HermesCommon
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
/*! \file tables.cpp
\brief Butcher tables. Including the class Table and enum ButcherTableType.
*/
#include "table.h"
#include "matrix.h"

using namespace Hermes::Algebra::DenseMatrixOperations;

namespace Hermes
{
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
    for (unsigned int i = 0; i < size; i++)
    {
      for (unsigned int j = 0; j < size; j++) this->A[i][j] = 0;
    }
  }

  void Table::alloc(unsigned int size)
  {
    // Size.
    this->size = size;
    // A array.
    this->A = new_matrix<double>(size, size);
    for (unsigned int i = 0; i < size; i++)
    {
      for (unsigned int j = 0; j < size; j++) this->A[i][j] = 0;
    }
  }

  unsigned int Table::get_size()
  {
    return this->size;
  }

  double Table::get_A(unsigned int i, unsigned int j)
  {
    if(i > size || j > size) throw Hermes::Exceptions::Exception("Invalid access to a Butcher's table.");
    return this->A[i][j];
  }

  void Table::set_A(unsigned int i, unsigned int j, double val)
  {
    if(i > size || j > size) throw Hermes::Exceptions::Exception("Invalid access to a Butcher's table.");
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
    for (unsigned int j = 0; j < size; j++) this->B[j] = 0;
    // B2 array.
    this->B2 = new double[size];
    for (unsigned int j = 0; j < size; j++) this->B2[j] = 0;
    // C array.
    this->C = new double[size];
    for (unsigned int j = 0; j < size; j++) this->C[j] = 0;
  }

  ButcherTable::ButcherTable(ButcherTableType butcher_table)
  {
    switch (butcher_table)
    {
      /* EXPLICIT METHODS */

      // Explicit Euler.
    case Explicit_RK_1:
      this->alloc(1);
      this->set_B(0, 1.);
      break;

      // Explicit RK-2.
    case Explicit_RK_2:
      this->alloc(2);
      this->set_A(1, 0, 2./3.);
      this->set_A(1, 1, 0.);
      this->set_B(0, 1./4.);
      this->set_B(1, 3./4.);
      this->set_C(0, 2./3.);
      break;

      // Explicit RK-3.
    case Explicit_RK_3:
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

      // Explicit RK-4.
    case Explicit_RK_4:
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

      /* IMPLICIT METHODS */

      // Implicit Euler.
    case Implicit_RK_1:
      this->alloc(1);
      this->set_A(0, 0, 1.);
      this->set_B(0, 1.);
      this->set_C(0, 1.);
      break;

      // Implicit Crank Nicolson.
    case Implicit_Crank_Nicolson_2_2:
      this->alloc(2);
      this->set_A(0, 0, 1./2.);
      this->set_A(0, 1, 1./2.);
      this->set_B(0, 1./2.);
      this->set_B(1, 1./2.);
      this->set_C(0, 1.);
      break;

      // Implicit SIRK-22 (second-order).
    case Implicit_SIRK_2_2:
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

      // Implicit ESIRK-22 (second-order).
    case Implicit_ESIRK_2_2:
      this->alloc(2);
      this->set_A(0, 0, (9. - 6.*sqrt(2.))/4.);
      this->set_A(0, 1, (-3. + 2.*sqrt(2.))/4.);
      this->set_A(1, 0, (11. - 6.*sqrt(2.))/4.);
      this->set_A(1, 1, (-1. + 2.*sqrt(2.))/4.);
      this->set_B(0, 2. - sqrt(2.0));
      this->set_B(1, -1. + sqrt(2.0));
      this->set_C(1, 1.);
      break;

      // Implicit SDIRK-2-2 (second-order).
    case Implicit_SDIRK_2_2:
      this->alloc(2);
      this->set_A(0, 0, 1. - 1./sqrt(2.));
      this->set_A(1, 0, 1./sqrt(2.));
      this->set_A(1, 1, 1. - 1./sqrt(2.));
      this->set_B(0, 1./sqrt(2.));
      this->set_B(1, 1. - 1./sqrt(2.));
      this->set_C(0, 1. - 1./sqrt(2.));
      this->set_C(1, 1.);
      break;

      // Implicit Lobatto IIIA (second-order).
    case Implicit_Lobatto_IIIA_2_2:
      this->alloc(2);
      this->set_A(1, 0, 1./2.);
      this->set_A(1, 1, 1./2.);
      this->set_B(0, 1./2.);
      this->set_B(1, 1./2.);
      this->set_C(1, 1.);
      break;

      // Implicit Lobatto IIIB (second-order).
    case Implicit_Lobatto_IIIB_2_2:
      this->alloc(2);
      this->set_A(0, 0, 1./2.);
      this->set_A(0, 1, 1./2.);
      this->set_B(0, 1./2.);
      this->set_B(1, 1./2.);
      this->set_C(0, 1./2.);
      this->set_C(1, 1./2.);
      break;

      // Implicit Lobatto IIIC (second-order).
    case Implicit_Lobatto_IIIC_2_2:
      this->alloc(2);
      this->set_A(0, 0, 1./2.);
      this->set_A(0, 1, -1./2.);
      this->set_A(1, 0, 1./2.);
      this->set_A(1, 1, 1./2.);
      this->set_B(0, 1./2.);
      this->set_B(1, 1./2.);
      this->set_C(1, 1.);
      break;

      // Implicit Lobatto IIIA (fourth-order).
    case Implicit_Lobatto_IIIA_3_4:
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

      // Implicit Lobatto IIIB (fourth-order).
    case Implicit_Lobatto_IIIB_3_4:
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

      // Implicit Lobatto IIIC (fourth-order).
    case Implicit_Lobatto_IIIC_3_4:
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

      // Implicit Radau IIA (fifth-order).
    case Implicit_Radau_IIA_3_5:
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

      // Implicit SDIRK-5-4 (fourth-order).
    case Implicit_SDIRK_5_4:
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

      /* EMBEDDED EXPLICIT METHODS */

      // Explicit Heun-Euler.
    case Explicit_HEUN_EULER_2_12_embedded:
      this->alloc(2);
      this->set_A(1, 0, 1.0);
      this->set_B(0, 0.5);
      this->set_B(1, 0.5);
      this->set_B2(0, 1.0);
      this->set_B2(1, 0.0);
      this->set_C(0, 0.0);
      this->set_C(1, 1.0);
      break;

      // Explicit Bogacki-Shampine.
    case Explicit_BOGACKI_SHAMPINE_4_23_embedded:
      this->alloc(4);
      this->set_A(1, 0, 1./2.);
      this->set_A(3, 0, 2./9.);
      this->set_A(2, 1, 3./4.);
      this->set_A(3, 1, 1./3.);
      this->set_A(3, 2, 4./9.);
      this->set_B(0, 2./9.);
      this->set_B(1, 1./3.);
      this->set_B(2, 4./9.);
      this->set_B(3, 0.0);
      this->set_B2(0, 7./24.);
      this->set_B2(1, 1./4.);
      this->set_B2(2, 1./3.);
      this->set_B2(3, 1./8.);
      this->set_C(1, 1./2.);
      this->set_C(2, 3./4.);
      this->set_C(3, 1.0);
      break;

      // Explicit Fehlberg.
    case Explicit_FEHLBERG_6_45_embedded:
      this->alloc(6);
      this->set_A(1, 0, 1./4.);
      this->set_A(2, 0, 3./32.);
      this->set_A(3, 0, 1932./2197.);
      this->set_A(4, 0, 439./216.);
      this->set_A(5, 0, -8./27.);
      this->set_A(2, 1, 9./32.);
      this->set_A(3, 1, -7200./2197.);
      this->set_A(4, 1, -8.);
      this->set_A(5, 1, 2.);
      this->set_A(3, 2, 7296./2197.);
      this->set_A(4, 2, 3680./513.);
      this->set_A(5, 2, -3544./2565.);
      this->set_A(4, 3, -845./4104.);
      this->set_A(5, 3, 1859./4104.);
      this->set_A(5, 4, -11./40.);
      this->set_B(0, 16./135.);
      this->set_B(1, 0);
      this->set_B(2, 6656./12825.);
      this->set_B(3, 28561./56430.);
      this->set_B(4, -9./50.);
      this->set_B(5, 2./55.);
      this->set_B2(0, 25./216.);
      this->set_B2(1, 0);
      this->set_B2(2, 1408./2565.);
      this->set_B2(3, 2197./4104.);
      this->set_B2(4, -1./5.);
      this->set_C(1, 1./4.);
      this->set_C(2, 3./8.);
      this->set_C(3, 12./13.);
      this->set_C(4, 1.);
      this->set_C(5, 1./2.);
      break;

      // Explicit Cash-Karp.
    case Explicit_CASH_KARP_6_45_embedded:
      this->alloc(6);
      this->set_A(1, 0, 1./5.);
      this->set_A(2, 0, 3./40.);
      this->set_A(3, 0, 3./10.);
      this->set_A(4, 0, -11./54.);
      this->set_A(5, 0, 1631./55296.);
      this->set_A(2, 1, 9./40.);
      this->set_A(3, 1, -9./10.);
      this->set_A(4, 1, 5./2.);
      this->set_A(5, 1, 175./512.);
      this->set_A(3, 2, 6./5.);
      this->set_A(4, 2, -70./27.);
      this->set_A(5, 2, 575./13824.);
      this->set_A(4, 3, 35./27.);
      this->set_A(5, 3, 44275./110592.);
      this->set_A(5, 4, 253./4096.);
      this->set_B(0, 37./378.);
      this->set_B(1, 0);
      this->set_B(2, 250./621.);
      this->set_B(3, 125./594.);
      this->set_B(4, 0);
      this->set_B(5, 512./1771.);
      this->set_B2(0, 2825./27648.);
      this->set_B2(1, 0);
      this->set_B2(2, 18575./48384.);
      this->set_B2(3, 13525./55296.);
      this->set_B2(4, 277./14336.);
      this->set_B2(5, 1./4.);
      this->set_C(1, 1./5.);
      this->set_C(2, 3./10.);
      this->set_C(3, 3./5.);
      this->set_C(4, 1.);
      this->set_C(5, 7./8.);
      break;

      // Explicit Dormand-Prince.
    case Explicit_DORMAND_PRINCE_7_45_embedded:
      this->alloc(7);
      this->set_A(1, 0, 1./5.);
      this->set_A(2, 0, 3./40.);
      this->set_A(3, 0, 44./45.);
      this->set_A(4, 0, 19372./6561.);
      this->set_A(5, 0, 9017./3168.);
      this->set_A(6, 0, 35./384.);
      this->set_A(2, 1, 9./40.);
      this->set_A(3, 1, -56./15.);
      this->set_A(4, 1, -25360./2187.);
      this->set_A(5, 1, -355./33.);
      this->set_A(6, 1, 0.);
      this->set_A(3, 2, 32./9.);
      this->set_A(4, 2, 64448./6561.);
      this->set_A(5, 2, 46732./5247.);
      this->set_A(6, 2, 500./1113.);
      this->set_A(4, 3, -212./729.);
      this->set_A(5, 3, 49./176.);
      this->set_A(6, 3, 125./192.);
      this->set_A(5, 4, -5103./18656.);
      this->set_A(6, 4, -2187./6784.);
      this->set_A(6, 5, 11./84.);
      this->set_B(0, 35./384.);
      this->set_B(1, 0.);
      this->set_B(2, 500./1113.);
      this->set_B(3, 125./192.);
      this->set_B(4, -2187./6784.);
      this->set_B(5, 11./84.);
      this->set_B(6, 0.);
      this->set_B2(0, 5179./57600.);
      this->set_B2(1, 0.);
      this->set_B2(2, 7571./16695.);
      this->set_B2(3, 393./640.);
      this->set_B2(4, -92097./339200.);
      this->set_B2(5, 187./2100.);
      this->set_B2(6, 1./40.);
      this->set_C(1, 1./5.);
      this->set_C(2, 3./10.);
      this->set_C(3, 4./5.);
      this->set_C(4, 8./9.);
      this->set_C(5, 1.);
      this->set_C(6, 1.);
      break;

      /* EMBEDDED IMPLICIT METHODS */

    case Implicit_ESDIRK_TRBDF2_3_23_embedded:
      this->alloc(3);
      this->set_A(1, 0, (2 - sqrt((double)2)) / 2.);
      this->set_A(2, 0, sqrt((double)2) / 4.);
      this->set_A(1, 1, (2 - sqrt((double)2)) / 2.);
      this->set_A(2, 1, sqrt((double)2) / 4.);
      this->set_A(2, 2, (2 - sqrt((double)2)) / 2.);
      this->set_B(0, sqrt((double)2) / 4.);
      this->set_B(1, sqrt((double)2) / 4.);
      this->set_B(2, (2 - sqrt((double)2)) / 2.);
      this->set_B2(0, (1 - sqrt((double)2)/4.) / 3.);
      this->set_B2(1, (3 * sqrt((double)2) / 4. + 1.) / 3.);
      this->set_B2(2, (2 - sqrt((double)2)) / 6.);
      this->set_C(1, 2 - sqrt((double)2));
      this->set_C(2, 1.0);
      break;

    case Implicit_ESDIRK_TRX2_3_23_embedded:
      this->alloc(3);
      this->set_A(1, 0, 0.25);
      this->set_A(2, 0, 0.25);
      this->set_A(1, 1, 0.25);
      this->set_A(2, 1, 0.5);
      this->set_A(2, 2, 0.25);
      this->set_B(0, 0.25);
      this->set_B(1, 0.5);
      this->set_B(2, 0.25);
      this->set_B2(0, 1./6.);
      this->set_B2(1, 2./3.);
      this->set_B2(2, 1./6.);
      this->set_C(1, 0.5);
      this->set_C(2, 1.0);
      break;

    case Implicit_SDIRK_CASH_3_23_embedded:
      this->alloc(3);
      this->set_A(0, 0, 0.435866521508);
      this->set_A(1, 0, 0.2820667320);
      this->set_A(2, 0, 1.208496649);
      this->set_A(1, 1, 0.435866521508);
      this->set_A(2, 1, -0.6443632015);
      this->set_A(2, 2, 0.435866521508);
      this->set_B(0, 1.208496649);
      this->set_B(1, -0.6443632015);
      this->set_B(2, 0.435866521508);
      this->set_B2(0, 0.77263013745746);
      this->set_B2(1, 0.22736986254254);
      this->set_C(0, 0.435866521508);
      this->set_C(1, 0.717933260755);
      this->set_C(2, 1.0);
      break;

    case Implicit_SDIRK_BILLINGTON_3_23_embedded:
      this->alloc(3);
      this->set_A(0, 0, 0.292893218813);
      this->set_A(1, 0, 0.798989873223);
      this->set_A(2, 0, 0.740789228841);
      this->set_A(1, 1, 0.292893218813);
      this->set_A(2, 1, 0.259210771159);
      this->set_A(2, 2, 0.292893218813);
      this->set_B(0, 0.691665115992);
      this->set_B(1, 0.503597029883);
      this->set_B(2, -0.195262145876);
      this->set_B2(0, 0.740789228840);
      this->set_B2(1, 0.259210771159);
      this->set_C(0, 0.292893218813);
      this->set_C(1, 1.091883092037);
      this->set_C(2, 1.292893218813);
      break;

    case Implicit_SDIRK_CASH_5_24_embedded:
      this->alloc(5);
      this->set_A(0, 0, 0.435866521508);
      this->set_A(1, 0, -1.13586652150);
      this->set_A(2, 0, 1.08543330679);
      this->set_A(3, 0, 0.416349501547);
      this->set_A(4, 0, 0.896869652944);
      this->set_A(1, 1, 0.435866521508);
      this->set_A(2, 1, -0.721299828287);
      this->set_A(3, 1, 0.190984004184);
      this->set_A(4, 1, 0.0182725272734);
      this->set_A(2, 2, 0.435866521508);
      this->set_A(3, 2, -0.118643265417);
      this->set_A(4, 2, -0.0845900310706);
      this->set_A(3, 3, 0.435866521508);
      this->set_A(4, 3, -0.266418670647);
      this->set_A(4, 4, 0.435866521508);
      this->set_B(0, 0.896869652944);
      this->set_B(1, 0.0182725272734);
      this->set_B(2, -0.0845900310706);
      this->set_B(3, -0.266418670647);
      this->set_B(4, 0.435866521508);
      this->set_B2(0, (-0.7 - 0.5) / (-0.7 - 0.435866521508));
      this->set_B2(1, (0.5 - 0.435866521508) / (-0.7 - 0.435866521508));
      this->set_B2(2, 0.0);
      this->set_B2(3, 0.0);
      this->set_B2(4, 0.0);
      this->set_C(0, 0.435866521508);
      this->set_C(1, -0.7);
      this->set_C(2, 0.8);
      this->set_C(3, 0.924556761814);
      this->set_C(4, 1.0);
      break;

    case Implicit_SDIRK_CASH_5_34_embedded:
      this->alloc(5);
      this->set_A(0, 0, 0.435866521508);
      this->set_A(1, 0, -1.13586652150);
      this->set_A(2, 0, 1.08543330679);
      this->set_A(3, 0, 0.416349501547);
      this->set_A(4, 0, 0.896869652944);
      this->set_A(1, 1, 0.435866521508);
      this->set_A(2, 1, -0.721299828287);
      this->set_A(3, 1, 0.190984004184);
      this->set_A(4, 1, 0.0182725272734);
      this->set_A(2, 2, 0.435866521508);
      this->set_A(3, 2, -0.118643265417);
      this->set_A(4, 2, -0.0845900310706);
      this->set_A(3, 3, 0.435866521508);
      this->set_A(4, 3, -0.266418670647);
      this->set_A(4, 4, 0.435866521508);
      this->set_B(0, 0.896869652944);
      this->set_B(1, 0.0182725272734);
      this->set_B(2, -0.0845900310706);
      this->set_B(3, -0.266418670647);
      this->set_B(4, 0.435866521508);
      this->set_B2(0, 0.776691932910);
      this->set_B2(1, 0.0297472791484);
      this->set_B2(2, -0.0267440239074);
      this->set_B2(3, 0.220304811849);
      this->set_B2(4, 0.0);
      this->set_C(0, 0.435866521508);
      this->set_C(1, -0.7);
      this->set_C(2, 0.8);
      this->set_C(3, 0.924556761814);
      this->set_C(4, 1.0);
      break;

    case Implicit_DIRK_ISMAIL_7_45_embedded: // Implicit embedded DIRK method with orders 4 and 5.
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
      this->set_A(4, 2, -0.027793233);
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
      this->set_B2(0, 0.094388663);
      this->set_B2(1, 0);
      this->set_B2(2, -0.039782614 );
      this->set_B2(3, 0.745608552);
      this->set_B2(4, -0.505129807);
      this->set_B2(5, 0.704915206);
      this->set_B2(6, 0);
      this->set_C(0, 0);
      this->set_C(1, 0.57178);
      this->set_C(2, 1.352846);
      this->set_C(3, 0.4);
      this->set_C(4, 0.75);
      this->set_C(5, 0.9);
      this->set_C(6, 1.0);
      break;

    default: throw Hermes::Exceptions::Exception("Unknown Butcher's table.");
    }
  }

  void ButcherTable::alloc(unsigned int size)
  {
    // Size.
    this->size = size;
    // A array.
    this->A = new_matrix<double>(size, size);
    for (unsigned int i = 0; i < size; i++)
    {
      for (unsigned int j = 0; j < size; j++) this->A[i][j] = 0;
    }
    // B array.
    this->B = new double[size];
    for (unsigned int j = 0; j < size; j++) this->B[j] = 0;
    // B2 array.
    this->B2 = new double[size];
    for (unsigned int j = 0; j < size; j++) this->B2[j] = 0;
    // C array.
    this->C = new double[size];
    for (unsigned int j = 0; j < size; j++) this->C[j] = 0;
  }

  double ButcherTable::get_B(unsigned int i)
  {
    if(i > size) throw Hermes::Exceptions::Exception("Invalid access to a Butcher's table.");
    return this->B[i];
  }

  double ButcherTable::get_B2(unsigned int i)
  {
    if(i > size) throw Hermes::Exceptions::Exception("Invalid access to a Butcher's table.");
    return this->B2[i];
  }

  double ButcherTable::get_C(unsigned int i)
  {
    if(i > size) throw Hermes::Exceptions::Exception("Invalid access to a Butcher's table.");
    return this->C[i];
  }

  void ButcherTable::set_B(unsigned int i, double val)
  {
    if(i > size) throw Hermes::Exceptions::Exception("Invalid access to a Butcher's table.");
    this->B[i] = val;
  }

  void ButcherTable::set_B2(unsigned int i, double val)
  {
    if(i > size) throw Hermes::Exceptions::Exception("Invalid access to a Butcher's table.");
    this->B2[i] = val;
  }

  void ButcherTable::set_C(unsigned int i, double val)
  {
    if(i > size) throw Hermes::Exceptions::Exception("Invalid access to a Butcher's table.");
    this->C[i] = val;
  }

  bool ButcherTable::is_explicit()
  {
    bool result = true;
    for (unsigned int i = 0; i<size; i++)
    {
      for (unsigned int j = 0; j<size; j++)
      {
        double val_ij = get_A(i, j);
        if(j >= i && fabs(val_ij) > Hermes::epsilon) result = false;
      }
    }

    return result;
  }

  bool ButcherTable::is_diagonally_implicit()
  {
    bool result = true;
    for (unsigned int i = 0; i < size; i++)
    {
      for (unsigned int j = 0; j < size; j++)
      {
        double val_ij = get_A(i, j);
        if(j > i && fabs(val_ij) > Hermes::epsilon) result = false;
      }
    }

    return result;
  }

  bool ButcherTable::is_fully_implicit()
  {
    bool result = false;
    for (unsigned int i = 0; i < size; i++)
    {
      for (unsigned int j = 0; j < size; j++)
      {
        double val_ij = get_A(i, j);
        if(j > i && fabs(val_ij) > Hermes::epsilon) result = true;
      }
    }

    return result;
  }

  bool ButcherTable::is_embedded()
  {
    // Test whether B2 row is not zero.
    double sum = 0;
    for (unsigned  int i = 0; i < size; i++) sum += fabs(B2[i]);
    if(sum < Hermes::epsilon) return false;
    else return true;
  }

  void ButcherTable::switch_B_rows()
  {
    // Test whether nonzero B2 row exists.
    if(this->is_embedded() == false)
      throw Hermes::Exceptions::Exception("ButcherTable::switch_B_rows(): Zero B2 row detected.");

    // Switch B rows.
    for (unsigned int i = 0; i < size; i++)
    {
      double tmp = B[i];
      B[i] = B2[i];
      B2[i] = tmp;
    }
  }
}