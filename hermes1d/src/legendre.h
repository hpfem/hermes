// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef SHAPESET_LEGENDRE_H_
#define SHAPESET_LEGENDRE_H_

#include <math.h>

#include "common.h"
#include "quad_std.h"

extern double leg_norm_const_ref(int n);
extern void fill_legendre_array_ref(double x, 
                                double val_array[MAX_P+1],
                                double der_array[MAX_P+1]);
extern double legendre_val_ref(double x, int n);
extern double legendre_der_ref(double x, int n);

// Poly orders of Legendre polynomials
extern int legendre_order_1d[];

// Precalculated values of Legendre polynomials and their derivatives 
// at all Gauss quadrature rules on the reference
// interval (-1, 1). The first index runs through Gauss quadrature 
// orders. The second index runs through the quadrature points of 
// the corresponding rule, and the third through the values of 
// Lobatto polynomials at that point. 
extern double legendre_val_ref_tab[MAX_QUAD_ORDER][MAX_QUAD_PTS_NUM][MAX_P + 1];
extern double legendre_der_ref_tab[MAX_QUAD_ORDER][MAX_QUAD_PTS_NUM][MAX_P + 1];
extern void precalculate_legendre_1d();

// Precalculated values of Legendre polynomials and their derivatives 
// (defined in (-1, 1)) at all Gauss quadrature rules which are 
// transformed to the interval (-1, 0).
// The first index runs through Gauss quadrature 
// orders. The second index runs through the quadrature points of 
// the corresponding rule, and the third through the values of 
// Lobatto polynomials at that point. 
extern double legendre_val_ref_tab_left[MAX_QUAD_ORDER][MAX_QUAD_PTS_NUM][MAX_P + 1];
extern double legendre_der_ref_tab_left[MAX_QUAD_ORDER][MAX_QUAD_PTS_NUM][MAX_P + 1];
extern void precalculate_legendre_1d_left();

// Precalculated values of Legendre polynomials and their derivatives 
// (defined in (-1, 1)) at all Gauss quadrature rules which are 
// transformed to the interval (0, 1).
// The first index runs through Gauss quadrature 
// orders. The second index runs through the quadrature points of 
// the corresponding rule, and the third through the values of 
// Lobatto polynomials at that point. 
extern double legendre_val_ref_tab_right[MAX_QUAD_ORDER][MAX_QUAD_PTS_NUM][MAX_P + 1];
extern double legendre_der_ref_tab_right[MAX_QUAD_ORDER][MAX_QUAD_PTS_NUM][MAX_P + 1];
extern void precalculate_legendre_1d_right();

// transforms point 'x_phys' from element (x1, x2) to (-1, 1)
double inverse_map(double x1, double x2, double x_phys);

// returns values of normalized Legendre polynomials on (a, b), for 
// an arbitrary point 'x'
double legendre_val_phys_plot(int i, double a, double b, double x);

// Returns values of normalized Legendre polynomials on (a, b) in
// Gauss quadrature points of order 'quad_order'.
// flag == 0: entire polynomial defined in interval (a,b)
// flag == -1: only left half of polynomial defined in interval (a,b)
// flag == 1: only right half of polynomial defined in interval (a,b)
void legendre_val_phys_quad(int flag, int quad_order, int fns_num, 
                       double a, double b,  
                       double leg_pol_val[MAX_QUAD_PTS_NUM][MAX_P+1]); 

// returns derivatives of normalized Legendre polynomials on (a, b), for
// an arbitrary point 'x'
double legendre_der_phys_plot(int i, double a, double b, double x);

// Returns derivatives of normalized Legendre polynomials on (a, b) in
// Gauss quadrature points of order 'quad_order'
// flag == 0: entire polynomial defined in interval (a,b)
// flag == -1: only left half of polynomial defined in interval (a,b)
// flag == 1: only right half of polynomial defined in interval (a,b)
void legendre_der_phys_quad(int flag, int quad_order, int fns_num, 
                       double a, double b,  
                       double leg_pol_der[MAX_QUAD_PTS_NUM][MAX_P+1]);

#endif /* SHAPESET_LEGENDRE_H_ */
