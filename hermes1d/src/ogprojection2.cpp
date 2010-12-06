// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "ogprojection.h"
#include "math.h"
#include "discrete_problem.h"

namespace {

// The following functions are meant to be used from this file only, so we put
// them into the anonymous namespace. This file should be merged with
// ogprojections.cpp eventually.

double L2_projection_biform(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
    double val = 0;
    for (int i=0; i<num; i++) {
        val += weights[i] * (u[i] * v[i]);
    }
    return val;
}

double H1_projection_biform(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
    double val = 0;
    for (int i=0; i<num; i++) {
        val += weights[i] * (u[i] * v[i] + dudx[i] * dvdx[i]);
    }
    return val;
}

void f_sin(int n, double x[], double f[], double dfdx[])
{
    for (int i=0; i<n; i++) {
        f[i] = sin(x[i]);
        if (dfdx != NULL)
            dfdx[i] = cos(x[i]);
    }
}

void f_unimplemented(int n, double x[], double f[], double dfdx[])
{
    error("Internal error: you need to implement _f");
}


ExactFunction _f = f_unimplemented;

double L2_projection_liform(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                double *v, double *dvdx, void *user_data)
{
    // Value of the projected function at gauss points:
    double* f = new double[num];
    _f(num, x, f, NULL);
    double val = 0;
    for (int i=0; i<num; i++) {
        val += weights[i] * (f[i] * v[i]);
    }
    delete [] f;
    return val;
}

double H1_projection_liform(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                double *v, double *dvdx, void *user_data)
{
    // Value of the projected function at gauss points:
    double* f = new double[num];
    double* dfdx = new double[num];
    _f(num, x, f, dfdx);
    double val = 0;
    for (int i=0; i<num; i++) {
        val += weights[i] * (f[i] * v[i] + dfdx[i] * dvdx[i]);
    }
    delete [] f;
    delete [] dfdx;
    return val;
}

} // anonymous namespace

void assemble_projection_matrix_rhs(Space *space, SparseMatrix *A, Vector *rhs,
        ExactFunction fn, int projection_type)
{
  WeakForm* wf = new WeakForm;
  if (projection_type == H1D_L2_ortho_global) {
      wf->add_matrix_form(L2_projection_biform);
      wf->add_vector_form(L2_projection_liform);
  } else if (projection_type == H1D_H1_ortho_global) {
      wf->add_matrix_form(H1_projection_biform);
      wf->add_vector_form(H1_projection_liform);
  } else
      throw std::runtime_error("Unknown projection type");

  DiscreteProblem *dp1 = new DiscreteProblem(wf, space);

  _f = fn;
  int N_dof = space->assign_dofs();
  info("Assembling projection linear system. ndofs: %d", N_dof);
  dp1->assemble(NULL, A, rhs);
  info("  Done assembling.");
  delete dp1;
}
