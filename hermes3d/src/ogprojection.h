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

#ifndef __H3D_OGPROJECTION_H
#define __H3D_OGPROJECTION_H

#include "space/space.h"
#include "discrete_problem.h"

class HERMES_API OGProjection
{
public:
  
  // global orthogonal projection
  static void project_global(Hermes::vector<Space *> spaces, 
                              Hermes::vector<Solution*> sols_src, Hermes::vector<Solution*> sols_dest, 
                              MatrixSolverType matrix_solver = SOLVER_UMFPACK, 
                              Hermes::vector<ProjNormType> proj_norms = Hermes::vector<ProjNormType>());
  static void project_global(Hermes::vector<Space *> spaces, Hermes::vector<MeshFunction *> source_meshfns, 
                              scalar* target_vec, MatrixSolverType matrix_solver = SOLVER_UMFPACK, 
                              Hermes::vector<ProjNormType> proj_norms = Hermes::vector<ProjNormType>());  
  static void project_global(Space* space, 
                              Solution* sol_src, Solution* sol_dest, 
                              MatrixSolverType matrix_solver, ProjNormType proj_norm);
protected:

  // Underlying function for global orthogonal projection.
  // Not intended for the user. NOTE: the weak form here must be 
  // a special projection weak form, which is different from 
  // the weak form of the PDE. If you supply a weak form of the 
  // PDE, the PDE will just be solved. 
  static void project_internal(Hermes::vector<Space *> spaces, WeakForm *proj_wf, scalar* target_vec, MatrixSolverType matrix_solver = SOLVER_UMFPACK);

  template<typename Real, typename Scalar>
  static Scalar H1projection_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
  {
    Scalar result = 0;
    for (int i = 0; i < n; i++)
      result += wt[i] * (u->val[i] * v->val[i] + 
                         u->dx[i] * v->dx[i] + 
                         u->dy[i] * v->dy[i] +
                         u->dz[i] * v->dz[i]);
    return result;
  }

  template<typename Real, typename Scalar>
  static Scalar H1projection_liform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
  {
    Scalar result = 0;
    for (int i = 0; i < n; i++)
      result += wt[i] * (ext->fn[0]->val[i] * v->val[i] + 
                         ext->fn[0]->dx[i] * v->dx[i] + 
                         ext->fn[0]->dy[i] * v->dy[i] + 
                         ext->fn[0]->dz[i] * v->dz[i]);
    return result;
  }

  template<typename Real, typename Scalar>
  static Scalar H1_semi_projection_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
  {
    Scalar result = 0;
    for (int i = 0; i < n; i++)
      result += wt[i] * (u->dx[i] * v->dx[i] + 
                         u->dy[i] * v->dy[i] + 
                         u->dz[i] * v->dz[i]);
    return result;
  }

  template<typename Real, typename Scalar>
  static Scalar H1_semi_projection_liform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
  {
    Scalar result = 0;
    for (int i = 0; i < n; i++)
      result += wt[i] * (ext->fn[0]->dx[i] * v->dx[i] + 
                         ext->fn[0]->dy[i] * v->dy[i] + 
                         ext->fn[0]->dz[i] * v->dz[i]);
    return result;
  }

  template<typename Real, typename Scalar>
  static Scalar L2projection_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
  {
    Scalar result = 0;
    for (int i = 0; i < n; i++)
      result += wt[i] * (u->val[i] * v->val[i]);
    return result;
  }

  template<typename Real, typename Scalar>
  static Scalar L2projection_liform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
  {
    Scalar result = 0;
    for (int i = 0; i < n; i++)
      result += wt[i] * (ext->fn[0]->val[i] * v->val[i]);
    return result;
  }

  // Hcurl projections
  template<typename Real, typename Scalar>
  static Scalar Hcurlprojection_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                                Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
  {
    Scalar result = 0;
    for (int i = 0; i < n; i++) {
      result += wt[i] * (u->curl0[i] * conj(v->curl0[i]) + u->curl1[i] * conj(v->curl1[i]) + u->curl2[i] * conj(v->curl2[i]));
      result += wt[i] * (u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]) + u->val2[i] * conj(v->val2[i]));
    }
    return result;
  }

  template<typename Real, typename Scalar>
  static Scalar Hcurlprojection_liform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                                Geom<Real> *e, ExtData<Scalar> *ext)
  {
    Scalar result = 0;
    for (int i = 0; i < n; i++) {
      result += wt[i] * (ext->fn[0]->curl0[i] * conj(v->curl0[i]) + ext->fn[0]->curl1[i] * conj(v->curl1[i]) + ext->fn[0]->curl2[i] * conj(v->curl2[i]));
      result += wt[i] * (ext->fn[0]->val0[i] * conj(v->val0[i]) + ext->fn[0]->val1[i] * conj(v->val1[i]) + ext->fn[0]->val2[i] * conj(v->val2[i]));
    }
    return result;
  }
};

#endif
