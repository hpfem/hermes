// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "quad_std.h"

Quad1DStd g_quad_1d_std;

// Gauss quadrature of order 'order' in (-1,1)
void create_ref_element_quadrature(int order, double *x_ref, 
                                   double *w_ref, int *pts_num) {
  double2 *ref_tab = g_quad_1d_std.get_points(order);
  *pts_num = g_quad_1d_std.get_num_points(order);
  for (int i=0;i<*pts_num;i++) {
    x_ref[i] = ref_tab[i][0]; 
    w_ref[i] = ref_tab[i][1];
  }
};

// Gauss quadrature of order 'order' in (a, b)
void create_phys_element_quadrature(double a, double b, 
                                    int order, double *x_phys, 
                                    double *w_phys, int *pts_num) {
  double2 *ref_tab = g_quad_1d_std.get_points(order);
  *pts_num = g_quad_1d_std.get_num_points(order);
  for (int i=0;i<*pts_num;i++) {
    //change points and weights from (-1, 1) to (a, b)
    x_phys[i] = (b-a)/2.*ref_tab[i][0] + (b+a)/2.; 
    w_phys[i] = ref_tab[i][1]*(b-a)/2.;
  }
};

Quad1DStd::Quad1DStd()
{
  tables = std_tables_1d;
  np = std_np_1d;
  ref_vert[0] = -1.0;
  ref_vert[1] = 1.0;
}

