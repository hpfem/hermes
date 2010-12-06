// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef _PROJECTION_H_
#define _PROJECTION_H_

#include "../../hermes_common/common.h"
#include "../../hermes_common/matrix.h"
#include "weakform.h"

/*
    Returns the values of the function in f[] and derivatives in dfdx[].

    This function type is used in calculation of the right hand side of the
    projection linear system.

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

class HERMES_API OGProjection
{
public:
  // Project the fine mesh solution (defined on space_ref) onto the coarse mesh (defined on space). Projects all solution components.
  static void project_global(Space *space, Space* space_ref, MatrixSolverType matrix_solver = SOLVER_UMFPACK, ProjNormType proj_norm = HERMES_H1_NORM);

  // Project the ExactFunction fun onto the Space space, and store it in [sln_to_save] component.
  static void project_global(Space *space, ExactFunction fun, int sln_to_save = 0, MatrixSolverType matrix_solver = SOLVER_UMFPACK, ProjNormType proj_norm = HERMES_H1_NORM);

protected:

  // Internal function facilitating the projection.
  static void project_internal(Space *space, MatrixSolverType matrix_solver, ProjNormType proj_norm, int sln_to_save = -1);

  // Standard projection forms.
  static double L2_projection_biform(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data);

  static double L2_projection_liform(int num, double *x, double *weights,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double *v, double *dvdx, void *user_data);

  static double H1_projection_biform(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data);

  static double H1_projection_liform(int num, double *x, double *weights,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double *v, double *dvdx, void *user_data);

  // Forms for different solution components.
  // Due to architecture of Space and its components, there are functions
  // for specific solution components.
  static double H1_projection_liform_0(int num, double *x, double *weights,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double *v, double *dvdx, void *user_data);
  static double H1_projection_liform_1(int num, double *x, double *weights,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double *v, double *dvdx, void *user_data);
  static double H1_projection_liform_2(int num, double *x, double *weights,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double *v, double *dvdx, void *user_data);
  static double H1_projection_liform_3(int num, double *x, double *weights,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double *v, double *dvdx, void *user_data);
  static double H1_projection_liform_4(int num, double *x, double *weights,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double *v, double *dvdx, void *user_data);
    
  static double L2_projection_liform_0(int num, double *x, double *weights,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double *v, double *dvdx, void *user_data);
  static double L2_projection_liform_1(int num, double *x, double *weights,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double *v, double *dvdx, void *user_data);
  static double L2_projection_liform_2(int num, double *x, double *weights,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double *v, double *dvdx, void *user_data);
  static double L2_projection_liform_3(int num, double *x, double *weights,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double *v, double *dvdx, void *user_data);
  static double L2_projection_liform_4(int num, double *x, double *weights,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double *v, double *dvdx, void *user_data);

  // Internal.
  // Solution component being changed by projection forms.
  static int sln;

  // Reference space used in the case of projecting solution defined on 
  // a reference mesh onto a coarse one.
  static Space* ref_space;

  // Current ExactFunction used for projections.
  static ExactFunction fn;
  
  // ExactFunction used in the case of projecting solution defined on 
  // a reference mesh onto a coarse one.
  static void ref_mesh_fn(int n, double x[], double f[], double dfdx[]);
};

#define H1D_L2_ortho_global 0
#define H1D_H1_ortho_global 1

void assemble_projection_matrix_rhs(Space *space, SparseMatrix *A, Vector *rhs,
        ExactFunction fn, int projection_type=H1D_L2_ortho_global);

#endif
