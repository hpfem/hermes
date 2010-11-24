// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef _PROJECTION_H_
#define _PROJECTION_H_

#include "../../hermes_common/common.h"
#include "../../hermes_common/matrix.h"
#include "weakform.h"

// Function type used in calculation of the right hand side of the projection linear system.
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
#endif
