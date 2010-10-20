// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef _PROJECTION_H_
#define _PROJECTION_H_

#include "common.h"
#include "mesh.h"
#include "matrix.h"

double L2_projection_biform(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data);

double L2_projection_liform(int num, double *x, double *weights,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double *v, double *dvdx, void *user_data);

double H1_projection_biform(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data);

double H1_projection_liform(int num, double *x, double *weights,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double *v, double *dvdx, void *user_data);

#define H1D_L2_ortho_global 0
#define H1D_H1_ortho_global 1

/*
    Returns the values of the function in f[] and derivatives in dfdx[].

    For all physical 'x' in x[]. If dfdx == NULL, don't return the derivative.

    This function is completely general and must work for any number of points
    'n'. Typically, this function is being called on each element, in which
    case the x[] are Gauss integratin points::

        ExactFunction f = <initialize it>;
        double x[10] = <initialize Gauss points>;
        double val[10], dfdx[10];
        f(10, x, val, dfdx);

    If you want to get just a value and derivative at a point, call it like::

        ExactFunction f = <initialize it>;
        double x = 3.15;
        double val, dfdx;
        f(1, &x, &val, &dfdx);

    or if you only need the value::

        ExactFunction f = <initialize it>;
        double x = 3.15;
        double val;
        f(1, &x, &val, NULL);
*/
typedef void(*ExactFunction)(int n, double x[], double f[], double dfdx[]);


void assemble_projection_matrix_rhs(Mesh *mesh, Matrix *A, double *rhs,
        ExactFunction fn, int projection_type=H1D_L2_ortho_global);

#endif
