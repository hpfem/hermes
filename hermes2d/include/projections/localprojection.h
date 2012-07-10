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

#ifndef __H2D_LOCALPROJECTION_H
#define __H2D_LOCALPROJECTION_H

#include "../function/solution.h"
#include "../forms.h"
#include "../weakform/weakform.h"
#include "../views/scalar_view.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// @ingroup projections
    template<typename Scalar>
    class HERMES_API LocalProjection
    {
    public:
      LocalProjection();

      // Main functionality.
      static void project_local(const Space<Scalar>* space, MeshFunction<Scalar>* meshfn,
          Scalar* target_vec, ProjNormType proj_norm = HERMES_UNSET_NORM);

      // Wrapper that delivers a Solution instead of coefficient vector.
      static void project_local(const Space<Scalar>* space,
    Solution<Scalar>* source_sln, Solution<Scalar>* target_sln, ProjNormType proj_norm = HERMES_UNSET_NORM);

      // Wrapper that takes multiple MeshFunctions.
      static void project_local(Hermes::vector<const Space<Scalar>*> spaces, Hermes::vector<MeshFunction<Scalar>*> meshfns,
          Scalar* target_vec, Hermes::vector<ProjNormType> proj_norms = Hermes::vector<ProjNormType>());

      // Wrapper that takes multiple Solutions.
      static void project_local(Hermes::vector<const Space<Scalar>*> spaces, Hermes::vector<Solution<Scalar>*> slns,
          Scalar* target_vec, Hermes::vector<ProjNormType> proj_norms = Hermes::vector<ProjNormType>());

      // Wrapper that delivers Solutions instead of a coefficient vector.
      static void project_local(Hermes::vector<const Space<Scalar>*> spaces,
          Hermes::vector<Solution<Scalar>*> source_slns, Hermes::vector<Solution<Scalar>*> target_slns,
          Hermes::vector<ProjNormType> proj_norms = Hermes::vector<ProjNormType>(), bool delete_old_mesh = false);

    protected:

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
            throw Hermes::Exceptions::Exception("Unknown projection type");
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
            throw Hermes::Exceptions::Exception("Unknown projection type");
            return Hermes::Ord();
          }
        }

      private:
        ProjNormType projNormType;

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        static SolFunctionDomain h1_projection_biform(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *u,
          Func<TestFunctionDomain> *v, Geom<TestFunctionDomain> *e, ExtData<SolFunctionDomain> *ext)
        {
          SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++)
            result += wt[i] * (u->val[i] * v->val[i] + u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
          return result;
        }

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        static SolFunctionDomain h1_semi_projection_biform(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *u,
          Func<TestFunctionDomain> *v, Geom<TestFunctionDomain> *e, ExtData<SolFunctionDomain> *ext)
        {
          SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++)
            result += wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
          return result;
        }

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        static SolFunctionDomain l2_projection_biform(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *u,
          Func<TestFunctionDomain> *v, Geom<TestFunctionDomain> *e, ExtData<SolFunctionDomain> *ext)
        {
          SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++)
            result += wt[i] * (u->val[i] * v->val[i]);
          return result;
        }

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        static SolFunctionDomain hcurl_projection_biform(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *u,
          Func<TestFunctionDomain> *v, Geom<TestFunctionDomain> *e, ExtData<SolFunctionDomain> *ext)
        {
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
            throw Hermes::Exceptions::Exception("Unknown projection type");
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
            throw Hermes::Exceptions::Exception("Unknown projection type");
            return Hermes::Ord();
          }
        }

      private:
        ProjNormType projNormType;

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        SolFunctionDomain h1_projection_residual(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *v,
          Geom<TestFunctionDomain> *e, ExtData<SolFunctionDomain> *ext) const
        {
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
          SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++)
            result += wt[i] * (u_ext[this->i]->val[i] - ext->fn[0]->val[i]) * v->val[i];
          return result;
        }

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        SolFunctionDomain hcurl_projection_residual(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *v,
          Geom<TestFunctionDomain> *e, ExtData<SolFunctionDomain> *ext) const
        {
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