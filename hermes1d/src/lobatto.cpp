// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "lobatto.h"
#include "legendre.h"

double lobatto_val_ref_tab[MAX_QUAD_ORDER][MAX_QUAD_PTS_NUM][MAX_P + 1];
double lobatto_der_ref_tab[MAX_QUAD_ORDER][MAX_QUAD_PTS_NUM][MAX_P + 1];
double lobatto_val_ref_tab_left[MAX_QUAD_ORDER][MAX_QUAD_PTS_NUM][MAX_P + 1];
double lobatto_der_ref_tab_left[MAX_QUAD_ORDER][MAX_QUAD_PTS_NUM][MAX_P + 1];
double lobatto_val_ref_tab_right[MAX_QUAD_ORDER][MAX_QUAD_PTS_NUM][MAX_P + 1];
double lobatto_der_ref_tab_right[MAX_QUAD_ORDER][MAX_QUAD_PTS_NUM][MAX_P + 1];

int lobatto_order_1d[] = {
1,
 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
11,12,13,14,15,16,17,18,19,20,
21,22,23,24,25,26,27,28,29,30,
31,32,33,34,35,36,37,38,39,40,
41,42,43,44,45,46,47,48,49,50,
51,52,53,54,55,56,57,58,59,60,
61,62,63,64,65,66,67,68,69,70,
71,72,73,74,75,76,77,78,79,80,
81,82,83,84,85,86,87,88,89,90,
91,92,93,94,95,96,97,98,99,100,
};

static double lobatto_fn_0(double x) {
    return 1.0/2.0 - x/2.;
}

static double lobatto_fn_1(double x) {
    return 1.0/2.0 + x/2.;
}

static double lobatto_der_0(double x) {
    return -1.0/2.0;
}

static double lobatto_der_1(double x) {
    return 1.0/2.0;
}

// Fills an array of length MAX_P + 1 with Lobatto shape 
// functions (integrated normalized Legendre polynomials) 
// at point 'x'. 
extern void fill_lobatto_array_ref(double x, 
                                   double lobatto_array_val[MAX_P+1],
                                   double lobatto_array_der[MAX_P+1]) {
    double legendre_array[MAX_P + 1];
    // calculating (non-normalized) Legendre polynomials
    legendre_array[0] = 1.;
    legendre_array[1] = x;
    for (int i=1; i < MAX_P; i++) {
      legendre_array[i+1]  = (2*i+1)*x*legendre_array[i] // last index is MAX_P
                             - i*legendre_array[i-1]; 
      legendre_array[i+1] /= i+1; 
    }
    // first fill the two linear Lobatto shape functions 
    lobatto_array_val[0] = lobatto_fn_0(x);
    lobatto_array_val[1] = lobatto_fn_1(x);
    lobatto_array_der[0] = lobatto_der_0(x);
    lobatto_array_der[1] = lobatto_der_1(x);
    // then fill the quadratic and higher which actually are 
    // the integrated Legendre polynomials
    for (int i=1; i < MAX_P; i++) {
      lobatto_array_val[i+1] =                           // last index is MAX_P
        (legendre_array[i+1] - legendre_array[i-1]) / (2.*i + 1.);
      lobatto_array_val[i+1] /= leg_norm_const_ref(i);
      lobatto_array_der[i+1] = legendre_array[i]; 
      lobatto_array_der[i+1] /= leg_norm_const_ref(i);
    }
}

// FIXME - this function is inefficient, it fills the
// whole array
extern double lobatto_val_ref(double x, int n) 
{
    double val_array[MAX_P + 1];
    double der_array[MAX_P + 1];
    fill_lobatto_array_ref(x, val_array, der_array);
    //printf("val_array[%d] = %g\n", n, val_array[n]);
    return val_array[n];
}

// FIXME - this function is inefficient, it fills the
// whole array
extern double lobatto_der_ref(double x, int n) 
{
    double val_array[MAX_P + 1];
    double der_array[MAX_P + 1];
    fill_lobatto_array_ref(x, val_array, der_array);
    return der_array[n];
}

// integrated Legendre polynomials in (-1, 1)
void precalculate_lobatto_1d() 
{
  // erasing
  for (int quad_order=0; quad_order < MAX_QUAD_ORDER; quad_order++) {
    for (int point_id=0; point_id < MAX_QUAD_PTS_NUM; point_id++) {
      for (int poly_deg=0; poly_deg < MAX_P + 1; poly_deg++) {
        lobatto_val_ref_tab[quad_order][point_id][poly_deg] = 0;
      }
    }
  }

  for (int quad_order=0; quad_order < MAX_QUAD_ORDER; quad_order++) {
    int pts_num = g_quad_1d_std.get_num_points(quad_order);
    double2 *ref_tab = g_quad_1d_std.get_points(quad_order);
    for (int point_id=0; point_id < pts_num; point_id++) {
      double x_ref = ref_tab[point_id][0];
      fill_lobatto_array_ref(x_ref, lobatto_val_ref_tab[quad_order][point_id],
                             lobatto_der_ref_tab[quad_order][point_id]);
    }
  }
}

// half-polynomials in (-1, 0)
void precalculate_lobatto_1d_left() 
{
  // erasing
  for (int quad_order=0; quad_order < MAX_QUAD_ORDER; quad_order++) {
    for (int point_id=0; point_id < MAX_QUAD_PTS_NUM; point_id++) {
      for (int poly_deg=0; poly_deg < MAX_P + 1; poly_deg++) {
        lobatto_val_ref_tab_left[quad_order][point_id][poly_deg] = 0;
        lobatto_der_ref_tab_left[quad_order][point_id][poly_deg] = 0;
      }
    }
  }

  for (int quad_order=0; quad_order < MAX_QUAD_ORDER; quad_order++) {
    int pts_num = g_quad_1d_std.get_num_points(quad_order);
    double2 *ref_tab = g_quad_1d_std.get_points(quad_order);
    for (int point_id=0; point_id < pts_num; point_id++) {
      double x_ref = (ref_tab[point_id][0] - 1.) / 2.;  // transf to (-1, 0)
      fill_lobatto_array_ref(x_ref, 
            lobatto_val_ref_tab_left[quad_order][point_id],
            lobatto_der_ref_tab_left[quad_order][point_id]);
    }
  }
}

// half-polynomials in (0, 1)
void precalculate_lobatto_1d_right() 
{
  // erasing
  for (int quad_order=0; quad_order < MAX_QUAD_ORDER; quad_order++) {
    for (int point_id=0; point_id < MAX_QUAD_PTS_NUM; point_id++) {
      for (int poly_deg=0; poly_deg < MAX_P + 1; poly_deg++) {
        lobatto_val_ref_tab_right[quad_order][point_id][poly_deg] = 0;
        lobatto_der_ref_tab_right[quad_order][point_id][poly_deg] = 0;
      }
    }
  }

  for (int quad_order=0; quad_order < MAX_QUAD_ORDER; quad_order++) {
    int pts_num = g_quad_1d_std.get_num_points(quad_order);
    double2 *ref_tab = g_quad_1d_std.get_points(quad_order);
    for (int point_id=0; point_id < pts_num; point_id++) {
      double x_ref = (ref_tab[point_id][0] + 1.) / 2.;  // transf to (0, 1)
      fill_lobatto_array_ref(x_ref, 
              lobatto_val_ref_tab_right[quad_order][point_id],
              lobatto_der_ref_tab_right[quad_order][point_id]);
    }
  }
}

