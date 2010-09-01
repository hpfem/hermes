// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#include "common.h"
#include "lobatto.h"

int lobatto_order_1d[] = { 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
int legendre_order_1d[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };

static double lobatto_fn_0(double x) { return l0(x); }
static double lobatto_fn_1(double x) { return l1(x); }
static double lobatto_fn_2(double x) { return l2(x); }
static double lobatto_fn_3(double x) { return l3(x); }
static double lobatto_fn_4(double x) { return l4(x); }
static double lobatto_fn_5(double x) { return l5(x); }
static double lobatto_fn_6(double x) { return l6(x); }
static double lobatto_fn_7(double x) { return l7(x); }
static double lobatto_fn_8(double x) { return l8(x); }
static double lobatto_fn_9(double x) { return l9(x); }
static double lobatto_fn_10(double x) { return l10(x); }
static double lobatto_fn_11(double x) { return l11(x); }

shape_fn_1d_t lobatto_fn_tab_1d[] = {
	lobatto_fn_0, lobatto_fn_1, lobatto_fn_2, lobatto_fn_3, lobatto_fn_4, lobatto_fn_5,
	lobatto_fn_6, lobatto_fn_7, lobatto_fn_8, lobatto_fn_9, lobatto_fn_10, lobatto_fn_11
};


static double lobatto_der_0(double x) { return dl0(x); }
static double lobatto_der_1(double x) { return dl1(x); }
static double lobatto_der_2(double x) { return dl2(x); }
static double lobatto_der_3(double x) { return dl3(x); }
static double lobatto_der_4(double x) { return dl4(x); }
static double lobatto_der_5(double x) { return dl5(x); }
static double lobatto_der_6(double x) { return dl6(x); }
static double lobatto_der_7(double x) { return dl7(x); }
static double lobatto_der_8(double x) { return dl8(x); }
static double lobatto_der_9(double x) { return dl9(x); }
static double lobatto_der_10(double x) { return dl10(x); }
static double lobatto_der_11(double x) { return dl11(x); }

shape_fn_1d_t lobatto_der_tab_1d[] = {
	lobatto_der_0, lobatto_der_1, lobatto_der_2, lobatto_der_3, lobatto_der_4, lobatto_der_5,
	lobatto_der_6, lobatto_der_7, lobatto_der_8, lobatto_der_9, lobatto_der_10, lobatto_der_11
};


static double legendre_fn_0(double x) { return  legendre0(x); }
static double legendre_fn_1(double x) { return  legendre1(x); }
static double legendre_fn_2(double x) { return  legendre2(x); }
static double legendre_fn_3(double x) { return  legendre3(x); }
static double legendre_fn_4(double x) { return  legendre4(x); }
static double legendre_fn_5(double x) { return  legendre5(x); }
static double legendre_fn_6(double x) { return  legendre6(x); }
static double legendre_fn_7(double x) { return  legendre7(x); }
static double legendre_fn_8(double x) { return  legendre8(x); }
static double legendre_fn_9(double x) { return  legendre9(x); }
static double legendre_fn_10(double x) { return  legendre10(x); }
static double legendre_fn_11(double x) { return  legendre11(x); }

shape_fn_1d_t legendre_fn_tab_1d[] = {
	legendre_fn_0, legendre_fn_1, legendre_fn_2, legendre_fn_3, legendre_fn_4, legendre_fn_5,
	legendre_fn_6, legendre_fn_7, legendre_fn_8, legendre_fn_9, legendre_fn_10, legendre_fn_11
};


static double legendre_der_0(double x) { return  legendre0x(x); }
static double legendre_der_1(double x) { return  legendre1x(x); }
static double legendre_der_2(double x) { return  legendre2x(x); }
static double legendre_der_3(double x) { return  legendre3x(x); }
static double legendre_der_4(double x) { return  legendre4x(x); }
static double legendre_der_5(double x) { return  legendre5x(x); }
static double legendre_der_6(double x) { return  legendre6x(x); }
static double legendre_der_7(double x) { return  legendre7x(x); }
static double legendre_der_8(double x) { return  legendre8x(x); }
static double legendre_der_9(double x) { return  legendre9x(x); }
static double legendre_der_10(double x) { return  legendre10x(x); }
static double legendre_der_11(double x) { return  legendre11x(x); }

shape_fn_1d_t legendre_der_tab_1d[] = {
	legendre_der_0, legendre_der_1, legendre_der_2, legendre_der_3, legendre_der_4, legendre_der_5,
	legendre_der_6, legendre_der_7, legendre_der_8, legendre_der_9, legendre_der_10, legendre_der_11
};
