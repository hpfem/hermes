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

#include "function/solution.h"
#include "forms.h"
#include "weakform/weakform.h"

namespace Hermes
{
  namespace Hermes2D
  {
    const ProjNormType HERMES_DEFAULT_PROJ_NORM = HERMES_H1_NORM;

    template<typename Scalar>
    class HERMES_API OGProjection
    {
    public:
      OGProjection();

      static void project_global(Hermes::vector<Space<Scalar>*> spaces, Hermes::vector<MeshFunction<Scalar>*> source_meshfns,
        Scalar* target_vec, Hermes::MatrixSolverType matrix_solver = SOLVER_UMFPACK,
        Hermes::vector<ProjNormType> proj_norms = Hermes::vector<ProjNormType>());

      static void project_global(Hermes::vector<Space<Scalar>*> spaces, Hermes::vector<Solution<Scalar>*> source_sols,
        Scalar* target_vec, Hermes::MatrixSolverType matrix_solver = SOLVER_UMFPACK,
        Hermes::vector<ProjNormType> proj_norms = Hermes::vector<ProjNormType>());

      static void project_global(Space<Scalar>* space, MeshFunction<Scalar>* source_meshfn,
        Scalar* target_vec, Hermes::MatrixSolverType matrix_solver = SOLVER_UMFPACK,
        ProjNormType proj_norm = HERMES_H1_NORM);

      static void project_global(Hermes::vector<Space<Scalar>*> spaces,
        Hermes::vector<Solution<Scalar>*> sols_src, Hermes::vector<Solution<Scalar>*> sols_dest,
        Hermes::MatrixSolverType matrix_solver = SOLVER_UMFPACK,
        Hermes::vector<ProjNormType> proj_norms = Hermes::vector<ProjNormType>(), bool delete_old_mesh = false);

      static void project_global(Space<Scalar>* space,
        Solution<Scalar>* sol_src, Solution<Scalar>* sol_dest,
        Hermes::MatrixSolverType matrix_solver = SOLVER_UMFPACK,
        ProjNormType proj_norm = HERMES_UNSET_NORM);

      static void project_global(Hermes::vector<Space<Scalar>*> spaces,
        Hermes::vector<MatrixFormVol<Scalar> *> custom_projection_jacobian,
        Hermes::vector<VectorFormVol<Scalar> *> custom_projection_residual,
        Scalar* target_vec, Hermes::MatrixSolverType matrix_solver = SOLVER_UMFPACK);

      static void project_global(Hermes::vector<Space<Scalar> *> spaces,
                                 Hermes::vector<MatrixFormVol<Scalar> *> custom_projection_jacobian,
                                 Hermes::vector<VectorFormVol<Scalar> *> custom_projection_residual,
                                 Hermes::vector<Solution<Scalar> *> sols_dest,
                                 Hermes::MatrixSolverType matrix_solver = Hermes::SOLVER_UMFPACK);

      static void project_global(Space<Scalar>* space,
                                 MatrixFormVol<Scalar>* custom_projection_jacobian,
                                 VectorFormVol<Scalar>* custom_projection_residual,
                                 Solution<Scalar>* sol_dest,
                                 Hermes::MatrixSolverType matrix_solver = Hermes::SOLVER_UMFPACK);

      // Underlying function for global orthogonal projection.
      // Not intended for the user. NOTE: the weak form here must be
      // a special projection weak form, which is different from
      // the weak form of the PDE. If you supply a weak form of the
      // PDE, the PDE will just be solved.
    protected:
      static void project_internal(Hermes::vector<Space<Scalar>*> spaces, WeakForm<Scalar>*proj_wf, Scalar* target_vec,
        Hermes::MatrixSolverType matrix_solver = SOLVER_UMFPACK);

      // Jacobian matrix (same as stiffness matrix since projections are linear).
      class ProjectionMatrixFormVol : public MatrixFormVol<Scalar>
      {
      public:
        ProjectionMatrixFormVol(int i, int j, ProjNormType projNormType) : MatrixFormVol<Scalar>(i, j)
        {
          this->projNormType = projNormType;
        }

        Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
          Geom<double> *e, ExtData<Scalar> *ext) const
        {
          switch (projNormType)
          {
          case HERMES_L2_NORM:
            return l2_projection_biform<double, Scalar>(n, wt, u_ext, u, v, e, ext);
          case HERMES_H1_NORM:
            return h1_projection_biform<double, Scalar>(n, wt, u_ext, u, v, e, ext);
          case HERMES_H1_SEMINORM:
            return h1_semi_projection_biform<double, Scalar>(n, wt, u_ext, u, v, e, ext);
          case HERMES_HCURL_NORM:
            return hcurl_projection_biform<double, Scalar>(n, wt, u_ext, u, v, e, ext);
          case HERMES_HDIV_NORM:
            return hdiv_projection_biform<double, Scalar>(n, wt, u_ext, u, v, e, ext);
          default:
            error("Unknown projection type");
            return 0.0;
          }
        }

        Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u, Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const
        {
          switch (projNormType)
          {
          case HERMES_L2_NORM:
            return l2_projection_biform<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, u, v, e, ext);
          case HERMES_H1_NORM:
            return h1_projection_biform<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, u, v, e, ext);
          case HERMES_H1_SEMINORM:
            return h1_semi_projection_biform<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, u, v, e, ext);
          case HERMES_HCURL_NORM:
            return hcurl_projection_biform<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, u, v, e, ext);
          case HERMES_HDIV_NORM:
            return hdiv_projection_biform<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, u, v, e, ext);
          default:
            error("Unknown projection type");
            return Hermes::Ord();
          }
        }

      private:
        ProjNormType projNormType;

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        static SolFunctionDomain h1_projection_biform(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *u,
          Func<TestFunctionDomain> *v, Geom<TestFunctionDomain> *e, ExtData<SolFunctionDomain> *ext)
        {
          _F_
            SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++)
            result += wt[i] * (u->val[i] * v->val[i] + u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
          return result;
        }

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        static SolFunctionDomain h1_semi_projection_biform(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *u,
          Func<TestFunctionDomain> *v, Geom<TestFunctionDomain> *e, ExtData<SolFunctionDomain> *ext)
        {
          _F_
            SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++)
            result += wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
          return result;
        }

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        static SolFunctionDomain l2_projection_biform(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *u,
          Func<TestFunctionDomain> *v, Geom<TestFunctionDomain> *e, ExtData<SolFunctionDomain> *ext)
        {
          _F_
            SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++)
            result += wt[i] * (u->val[i] * v->val[i]);
          return result;
        }

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        static SolFunctionDomain hcurl_projection_biform(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *u,
          Func<TestFunctionDomain> *v, Geom<TestFunctionDomain> *e, ExtData<SolFunctionDomain> *ext)
        {
          _F_
            SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++) {
            result += wt[i] * (u->curl[i] * conj(v->curl[i]));
            result += wt[i] * (u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
          }
          return result;
        }

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        static SolFunctionDomain hdiv_projection_biform(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *u,
          Func<TestFunctionDomain> *v, Geom<TestFunctionDomain> *e, ExtData<SolFunctionDomain> *ext)
        {
          _F_
            SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++) {
            result += wt[i] * (u->div[i] * conj(v->div[i]));
            result += wt[i] * (u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
          }
          return result;
        }
      };

      // Residual.
      class ProjectionVectorFormVol : public VectorFormVol<Scalar>
      {
      public:
        ProjectionVectorFormVol(int i, MeshFunction<Scalar>* ext, ProjNormType projNormType) : VectorFormVol<Scalar>(i)
        {
          this->projNormType = projNormType;
          this->ext = Hermes::vector<MeshFunction<Scalar>*>();
          this->ext.push_back(ext);
        }

        Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
          Geom<double> *e, ExtData<Scalar> *ext) const
        {
          switch (projNormType)
          {
          case HERMES_L2_NORM:
            return l2_projection_residual<double, Scalar>(n, wt, u_ext, v, e, ext);
          case HERMES_H1_NORM:
            return h1_projection_residual<double, Scalar>(n, wt, u_ext, v, e, ext);
          case HERMES_H1_SEMINORM:
            return h1_semi_projection_residual<double, Scalar>(n, wt, u_ext, v, e, ext);
          case HERMES_HCURL_NORM:
            return hcurl_projection_residual<double, Scalar>(n, wt, u_ext, v, e, ext);
          case HERMES_HDIV_NORM:
            return hdiv_projection_residual<double, Scalar>(n, wt, u_ext, v, e, ext);
          default:
            error("Unknown projection type");
            return 0.0;
          }
        }

        Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const
        {
          switch (projNormType)
          {
          case HERMES_L2_NORM:
            return l2_projection_residual<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, v, e, ext);
          case HERMES_H1_NORM:
            return h1_projection_residual<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, v, e, ext);
          case HERMES_H1_SEMINORM:
            return h1_semi_projection_residual<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, v, e, ext);
          case HERMES_HCURL_NORM:
            return hcurl_projection_residual<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, v, e, ext);
          case HERMES_HDIV_NORM:
            return hdiv_projection_residual<Hermes::Ord, Hermes::Ord>(n, wt, u_ext, v, e, ext);
          default:
            error("Unknown projection type");
            return Hermes::Ord();
          }
        }

      private:
        ProjNormType projNormType;

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        SolFunctionDomain h1_projection_residual(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *v,
          Geom<TestFunctionDomain> *e, ExtData<SolFunctionDomain> *ext) const
        {
          _F_
            SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++)
            result += wt[i] * ((u_ext[this->i]->val[i] - ext->fn[0]->val[i]) * v->val[i]
          + (u_ext[this->i]->dx[i] - ext->fn[0]->dx[i]) * v->dx[i]
          + (u_ext[this->i]->dy[i] - ext->fn[0]->dy[i]) * v->dy[i]);
          return result;
        }

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        SolFunctionDomain h1_semi_projection_residual(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *v,
          Geom<TestFunctionDomain> *e, ExtData<SolFunctionDomain> *ext) const
        {
          _F_
            SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++)
            result += wt[i] * ((u_ext[this->i]->dx[i] - ext->fn[0]->dx[i]) * v->dx[i]
          + (u_ext[this->i]->dy[i] - ext->fn[0]->dy[i]) * v->dy[i]);
          return result;
        }

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        SolFunctionDomain l2_projection_residual(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *v,
          Geom<TestFunctionDomain> *e, ExtData<SolFunctionDomain> *ext) const
        {
          _F_
            SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++)
            result += wt[i] * (u_ext[this->i]->val[i] - ext->fn[0]->val[i]) * v->val[i];
          return result;
        }

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        SolFunctionDomain hcurl_projection_residual(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *v,
          Geom<TestFunctionDomain> *e, ExtData<SolFunctionDomain> *ext) const
        {
          _F_
            SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++) {
            result += wt[i] * (u_ext[this->i]->curl[i] - ext->fn[0]->curl[i]) * conj(v->curl[i]);
            result += wt[i] * ((u_ext[this->i]->val0[i] - ext->fn[0]->val0[i]) * conj(v->val0[i])
              + (u_ext[this->i]->val1[i] - ext->fn[0]->val1[i]) * conj(v->val1[i]));
          }

          return result;
        }

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        SolFunctionDomain hdiv_projection_residual(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *v,
          Geom<TestFunctionDomain> *e, ExtData<SolFunctionDomain> *ext) const
        {
          _F_
            SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++) {
            result += wt[i] * (u_ext[this->i]->div[i] - ext->fn[0]->div[i]) * conj(v->div[i]);
            result += wt[i] * ((u_ext[this->i]->val0[i] - ext->fn[0]->val0[i]) * conj(v->val0[i])
              + (u_ext[this->i]->val1[i] - ext->fn[0]->val1[i]) * conj(v->val1[i]));
          }

          return result;
        }
      };

      static int ndof;
    };
  }
}
#endif
