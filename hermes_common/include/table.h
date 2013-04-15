// This file is part of HermesCommon
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
/*! \file tables.h
\brief Butcher tables. Including the class Table and enum ButcherTableType.
*/
#ifndef __HERMES_COMMON_TABLES_H_
#define __HERMES_COMMON_TABLES_H_

#include "compat.h"

namespace Hermes
{
  /// \brief Butcher's tables type.
  /// The last number in the name always means order,
  /// the one before last (if provided) is the number of stages.
  enum ButcherTableType
  {
    /* EXPLICIT METHODS */

    Explicit_RK_1,               ///< Explicit Runge-Kutta RK-1 (explicit Euler).
    Explicit_RK_2,               ///< Explicit Runge-Kutta RK-2 method.
    Explicit_RK_3,               ///< Explicit Runge-Kutta RK-3 method.
    Explicit_RK_4,               ///< Explicit Runge-Kutta RK-4 method.

    /* IMPLICIT METHODS */

    Implicit_RK_1,               ///< Implicit Runge-Kutta RK-1 (implicit Euler).
    Implicit_Crank_Nicolson_2_2, ///< Implicit Crank_Nicolson method.
    Implicit_SIRK_2_2,           ///< Implicit SIRK-2-2 method.
    Implicit_ESIRK_2_2,          ///< Implicit ESIRK-2-2 method.
    Implicit_SDIRK_2_2,          ///< Implicit SDIRK-2-2 method.
    Implicit_Lobatto_IIIA_2_2,   ///< Implicit Lobatto IIIA-2 method.
    Implicit_Lobatto_IIIB_2_2,   ///< Implicit Lobatto IIIB-2 method.
    Implicit_Lobatto_IIIC_2_2,   ///< Implicit Lobatto IIIB-2 method.
    Implicit_Lobatto_IIIA_3_4,   ///< Implicit Lobatto IIIA-4 method.
    Implicit_Lobatto_IIIB_3_4,   ///< Implicit Lobatto IIIB-4 method.
    Implicit_Lobatto_IIIC_3_4,   ///< Implicit Lobatto IIIB-4 method.
    Implicit_Radau_IIA_3_5,      ///< Implicit Radau IIA-5 method.
    Implicit_SDIRK_5_4,          ///< Implicit SDIRK-5-4 method.

    /* EMBEDDED EXPLICIT METHODS */

    Explicit_HEUN_EULER_2_12_embedded,          ///< Heun-Euler (orders 1 and 2),
    Explicit_BOGACKI_SHAMPINE_4_23_embedded,    ///< Bogacki-Shampine (orders 2 and 3),
    Explicit_FEHLBERG_6_45_embedded,            ///< Fehlberg (orders 4 and 5),
    Explicit_CASH_KARP_6_45_embedded,           ///< Cash-Karp (orders 4 and 5),
    Explicit_DORMAND_PRINCE_7_45_embedded,      ///< Dormand-Prince (orders 4 and 5),

    /* EMBEDDED IMPLICIT METHODS */

    Implicit_ESDIRK_TRBDF2_3_23_embedded,    ///< From the paper Analysis and implementation
    ///< of TR-BDF2 by Hosea and Shampine.
    Implicit_ESDIRK_TRX2_3_23_embedded,      ///< From the paper Analysis and implementation
    ///< of TR-BDF2 by Hosea and Shampine.
    Implicit_SDIRK_CASH_3_23_embedded,       ///< From the paper Diagonally Implicit Runge-Kutta Formulae
    ///< with Error Estimates by J.R. Cash
    Implicit_SDIRK_BILLINGTON_3_23_embedded, ///< From the Master Thesis by S.R. Billington: Type Insensitive
    ///< Codes for the Solution of Stiff and Non-Stiff systems of ODE's
    Implicit_SDIRK_CASH_5_24_embedded,       ///< From the paper Diagonally Implicit Runge-Kutta Formulae
    ///< with Error Estimates by J.R. Cash
    Implicit_SDIRK_CASH_5_34_embedded,       ///< From the paper Diagonally Implicit Runge-Kutta Formulae
    ///< with Error Estimates by J.R. Cash
    Implicit_DIRK_ISMAIL_7_45_embedded       ///< Implicit embedded DIRK method pair of orders four in five (from the paper
    ///< Fudziah Ismail et all: Embedded Pair of Diagonally Implicit Runge-Kutta
    ///< Method for Solving Ordinary Differential Equations). The method has
    ///< 7 stages but the first one is explicit.
  };

  /// \brief General square table of real numbers.
  class HERMES_API Table
  {
  public:
    Table();
    Table(unsigned int size);
    virtual void alloc(unsigned int size);
    unsigned int get_size();
    double get_A(unsigned int i, unsigned int j);
    void set_A(unsigned int i, unsigned int j, double val);

  protected:
    unsigned int size;
    double** A;
  };

  /// \brief Butcher's tables for Runge-Kutta methods.
  class HERMES_API ButcherTable: public Table
  {
  public:
    ButcherTable();
    ButcherTable(unsigned int size);
    ButcherTable(ButcherTableType butcher_table);
    virtual void alloc(unsigned int size);
    double get_B(unsigned int i);
    double get_B2(unsigned int i);
    double get_C(unsigned int i);
    void set_B(unsigned int i, double val);
    void set_B2(unsigned int i, double val);
    void set_C(unsigned int i, double val);
    bool is_explicit();
    bool is_diagonally_implicit();
    bool is_fully_implicit();
    bool is_embedded();
    void switch_B_rows(); ///< For experimental purposes. Switches the B and B2 rows. B2 row
    ///< must be nonzero, otherwise error is thrown.

  protected:
    double* B;
    double* B2;  ///< This is the second B-row for adaptivity based
    ///< on embedded R-K methods.
    double* C;
  };
}
#endif