// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.


#ifndef __H2D_OGPROJECTION_H
#define __H2D_OGPROJECTION_H

#include "discrete_problem.h"

class HERMES_API OGProjection
{
public:
  static void project_global(Hermes::Tuple<Space *> spaces, Hermes::Tuple<MeshFunction *> source_meshfns, 
                              scalar* target_vec, MatrixSolverType matrix_solver = SOLVER_UMFPACK, Hermes::Tuple<ProjNormType> proj_norms = Hermes::Tuple<ProjNormType>());

  static void project_global(Hermes::Tuple<Space *> spaces, 
                              Hermes::Tuple<Solution*> sols_src, Hermes::Tuple<Solution*> sols_dest, 
                              MatrixSolverType matrix_solver = SOLVER_UMFPACK, 
                              Hermes::Tuple<ProjNormType> proj_norms = Hermes::Tuple<ProjNormType>());

  static void project_global(Hermes::Tuple<Space *> spaces, Hermes::Tuple< std::pair<WeakForm::matrix_form_val_t, WeakForm::matrix_form_ord_t> > proj_biforms, 
                      Hermes::Tuple< std::pair<WeakForm::vector_form_val_t, WeakForm::vector_form_ord_t> > proj_liforms, Hermes::Tuple<MeshFunction*> source_meshfns, 
                      scalar* target_vec, MatrixSolverType matrix_solver = SOLVER_UMFPACK);

  static void project_global(Space *space, 
                      std::pair<WeakForm::matrix_form_val_t, WeakForm::matrix_form_ord_t> proj_biform,
                      std::pair<WeakForm::vector_form_val_t, WeakForm::vector_form_ord_t> proj_liform,
                      ExactFunction source_fn, scalar* target_vec, MatrixSolverType matrix_solver = SOLVER_UMFPACK);

  /// Global orthogonal projection of one vector-valued ExactFunction.
  static void project_global(Space *space, ExactFunction2 source_fn, scalar* target_vec, MatrixSolverType matrix_solver = SOLVER_UMFPACK);

  /// Global orthogonal projection of one scalar-valued ExactFunction.
  static void project_global(Space *space, ExactFunction source_fn, scalar* target_vec, MatrixSolverType matrix_solver = SOLVER_UMFPACK);

  /// Projection-based interpolation of an exact function. This is faster than the
  /// global projection since no global matrix problem is solved.
  static void project_local(Space *space, int proj_norm, ExactFunction source_fn, Mesh* mesh,
                   scalar* target_vec);

  // Underlying function for global orthogonal projection.
  // Not intended for the user. NOTE: the weak form here must be 
  // a special projection weak form, which is different from 
  // the weak form of the PDE. If you supply a weak form of the 
  // PDE, the PDE will just be solved. 
protected:
  static void project_internal(Hermes::Tuple<Space *> spaces, WeakForm *proj_wf, scalar* target_vec, MatrixSolverType matrix_solver = SOLVER_UMFPACK);

  // The projection functionality below is identical in H2D and H3D.
  template<typename Real, typename Scalar>
  static Scalar H1projection_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
  {
    _F_
    Scalar result = 0;
    for (int i = 0; i < n; i++)
      result += wt[i] * (u->val[i] * v->val[i] + u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
    return result;
  }

  template<typename Real, typename Scalar>
  static Scalar H1_semi_projection_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
  {
    _F_
    Scalar result = 0;
    for (int i = 0; i < n; i++)
      result += wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
    return result;
  }

  template<typename Real, typename Scalar>
  static Scalar H1projection_liform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
  {
    _F_
    Scalar result = 0;
    for (int i = 0; i < n; i++)
      result += wt[i] * (ext->fn[0]->val[i] * v->val[i] + ext->fn[0]->dx[i] * v->dx[i] + ext->fn[0]->dy[i] * v->dy[i]);
    return result;
  }

  template<typename Real, typename Scalar>
  static Scalar H1_semi_projection_liform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
  {
    _F_
    Scalar result = 0;
    for (int i = 0; i < n; i++)
      result += wt[i] * (ext->fn[0]->dx[i] * v->dx[i] + ext->fn[0]->dy[i] * v->dy[i]);
    return result;
  }

  template<typename Real, typename Scalar>
  static Scalar L2projection_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
  {
    _F_
    Scalar result = 0;
    for (int i = 0; i < n; i++)
      result += wt[i] * (u->val[i] * v->val[i]);
    return result;
  }

  template<typename Real, typename Scalar>
  static Scalar L2projection_liform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
  {
    _F_
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
    _F_
    Scalar result = 0;
    for (int i = 0; i < n; i++) {
      result += wt[i] * (u->curl[i] * conj(v->curl[i]));
      result += wt[i] * (u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
    }
    return result;
  }

  template<typename Real, typename Scalar>
  static Scalar Hcurlprojection_liform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                                Geom<Real> *e, ExtData<Scalar> *ext)
  {
    _F_
    Scalar result = 0;
    for (int i = 0; i < n; i++) {
      result += wt[i] * (ext->fn[0]->curl[i] * conj(v->curl[i]));
      result += wt[i] * (ext->fn[0]->val0[i] * conj(v->val0[i]) + ext->fn[0]->val1[i] * conj(v->val1[i]));
    }

    return result;
  }
 
};

#endif
