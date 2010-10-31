// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef _PROJECTION_H_
#define _PROJECTION_H_

#include "../../hermes_common/common.h"
#include "../../hermes_common/matrix.h"
#include "weakform.h"

typedef void(*ExactFunction)(int n, double x[], double f[], double dfdx[]);

class HERMES_API OGProjection
{
public:
  static void project_global(Space *space, Space* space_ref, MatrixSolverType matrix_solver = SOLVER_UMFPACK, ProjNormType proj_norm = HERMES_H1_NORM);

protected:
  static void project_global(int sln, Space *space, Space* space_ref, MatrixSolverType matrix_solver = SOLVER_UMFPACK, ProjNormType proj_norm = HERMES_H1_NORM);

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

  // Internal.
  static int sln;

};
#endif
