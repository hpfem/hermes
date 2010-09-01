// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef SHAPESET_LOBATTO_H_
#define SHAPESET_LOBATTO_H_

#include <math.h>

#include "common.h"

void fill_lobatto_array_ref(double x, 
			double lobatto_array_val[MAX_P+1],
			double lobatto_array_der[MAX_P+1]);
double lobatto_val_ref(double x, int n);
double lobatto_der_ref(double x, int n);

// Poly orders of Lobatto functions
extern int lobatto_order_1d[];

// Precalculated values of Lobatto polynomials and their derivatives 
// at all Gauss quadrature rules on the reference interval (-1, 1). 
// The first index runs through Gauss quadrature 
// orders. The second index runs through the quadrature points of 
// the corresponding rule, and the third through the values of 
// Legendre polynomials at that point. 
extern double lobatto_val_ref_tab[MAX_QUAD_ORDER][MAX_QUAD_PTS_NUM][MAX_P + 1];
extern double lobatto_der_ref_tab[MAX_QUAD_ORDER][MAX_QUAD_PTS_NUM][MAX_P + 1];
extern void precalculate_lobatto_1d();

// The first index runs through Gauss quadrature 
// Precalculated values of Lobatto polynomials and their derivatives 
// defined in (-1, 1) at all Gauss quadrature rules which are 
// transformed to the interval (-1, 0).
// The second index runs through the quadrature points of 
// the corresponding rule, and the third through the values of 
// Legendre polynomials at that point. 
extern double lobatto_val_ref_tab_left[MAX_QUAD_ORDER][MAX_QUAD_PTS_NUM][MAX_P + 1];
extern double lobatto_der_ref_tab_left[MAX_QUAD_ORDER][MAX_QUAD_PTS_NUM][MAX_P + 1];
extern void precalculate_lobatto_1d_left();

// Precalculated values of Lobatto polynomials and their derivatives 
// defined in (-1, 1) at all Gauss quadrature rules which are 
// transformed to the interval (0, 1).
// The first index runs through Gauss quadrature 
// orders. The second index runs through the quadrature points of 
// the corresponding rule, and the third through the values of 
// Legendre polynomials at that point. 
extern double lobatto_val_ref_tab_right[MAX_QUAD_ORDER][MAX_QUAD_PTS_NUM][MAX_P + 1];
extern double lobatto_der_ref_tab_right[MAX_QUAD_ORDER][MAX_QUAD_PTS_NUM][MAX_P + 1];
extern void precalculate_lobatto_1d_right();

#endif /* SHAPESET_LOBATTO_H_ */
