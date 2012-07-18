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

#include "../function/solution.h"
#include "../forms.h"
#include "../weakform/weakform.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /** @defgroup projections Projections
    * \brief Projection classes for various kinds of projecting a MeshFunction onto a Space.
    */

    /// @ingroup projections
    /// \brief Class for (global) orthogonal projecting. If the projection is not necessary (if a solution belongs to the space), then its solution vector is used.
    template<typename Scalar>
    class HERMES_API OGProjection : public Hermes::Mixins::Loggable
    {
    public:
      OGProjection();

      /// Main functionality is in the protected method project_internal().
      
      /// This method allows to specify your own OG-projection form.
      void project_global(const Space<Scalar>* space,
          MatrixFormVol<Scalar>* custom_projection_jacobian,
          VectorFormVol<Scalar>* custom_projection_residual,
          Scalar* target_vec);

      /// Wrapper that delivers a Solution instead of a coefficient vector.   
      void project_global(const Space<Scalar>* space,
          MatrixFormVol<Scalar>* custom_projection_jacobian,
          VectorFormVol<Scalar>* custom_projection_residual,
          Solution<Scalar>* target_sln, 
          double newton_tol = 1e-6, int newton_max_iter = 10);
      
      /// This method allows to specify your own multiple OG-projection forms.
      void project_global(const Hermes::vector<const Space<Scalar>*>& spaces,
          const Hermes::vector<MatrixFormVol<Scalar>*>& custom_projection_jacobian,
          const Hermes::vector<VectorFormVol<Scalar>*>& custom_projection_residual,
          Scalar* target_vec,
          double newton_tol = 1e-6, int newton_max_iter = 10);
          
      /// Wrapper that delivers a vector of Solutions instead of a coefficient vector.   
      void project_global(const Hermes::vector<const Space<Scalar>*>& spaces,
          const Hermes::vector<MatrixFormVol<Scalar>*>& custom_projection_jacobian,
          const Hermes::vector<VectorFormVol<Scalar>*>& custom_projection_residual,
          const Hermes::vector<Solution<Scalar>*>& target_slns,
          double newton_tol = 1e-6, int newton_max_iter = 10);
          
      /**
       \fn  static void OGProjection::project_global(Space<Scalar>* space,
        MeshFunction<Scalar>* source_meshfn, Scalar* target_vec,
        ProjNormType proj_norm = HERMES_UNSET_NORM, double newton_tol = 1e-6, int newton_max_iter = 10);

       \brief The method checks source_meshfn if it is an instance of Solution, if so, it checks its sln_vector, and space_seq
              if they can be used directly.

       \author  LK
       \date  10/29/2011

       \param[in]  space         If non-null, the space.
       \param[in]  source_meshfn If non-null, source meshfn.
       \param[out]  target_vec    If non-null, target vector.
       \param matrix_solver           (optional) the matrix solver.
       \param proj_norm               (optional) the project normalise.
       \param newton_tol              (optional) the newton tolerance.
       \param newton_max_iter         (optional) the newton maximum iterator.
       */
      void project_global(const Space<Scalar>* space, MeshFunction<Scalar>* source_meshfn,
          Scalar* target_vec, ProjNormType proj_norm = HERMES_UNSET_NORM);

      /// Wrapper that delivers a Solution instead of coefficient vector.
      void project_global(const Space<Scalar>* space,
          Solution<Scalar>* source_sln, Solution<Scalar>* target_sln,
          ProjNormType proj_norm = HERMES_UNSET_NORM);

      /// Wrapper for multiple source MeshFunctions that delivers coefficient vector.
      void project_global(Hermes::vector<const Space<Scalar>*> spaces, Hermes::vector<MeshFunction<Scalar>*> source_meshfns,
          Scalar* target_vec, Hermes::vector<ProjNormType> proj_norms = Hermes::vector<ProjNormType>());

      /// Wrapper for multiple source Solutions that delivers coefficient vector.
      void project_global(Hermes::vector<const Space<Scalar>*> spaces, Hermes::vector<Solution<Scalar>*> source_slns,
          Scalar* target_vec, Hermes::vector<ProjNormType> proj_norms = Hermes::vector<ProjNormType>());

      void project_global(Hermes::vector<const Space<Scalar>*> spaces,
          Hermes::vector<Solution<Scalar>*> source_slns, Hermes::vector<Solution<Scalar>*> target_slns,
          Hermes::vector<ProjNormType> proj_norms = Hermes::vector<ProjNormType>(), bool delete_old_mesh = false);

    protected:
      /// Underlying function for global orthogonal projection.
      /// Not intended for the user. NOTE: the weak form here must be
      /// a special projection weak form, which is different from
      /// the weak form of the PDE. If you supply a weak form of the
      /// PDE, the PDE will just be solved.
      void project_internal(const Space<Scalar>* space, WeakForm<Scalar>* proj_wf, Scalar* target_vec);

      /// Jacobian matrix (same as stiffness matrix since projections are linear).
      class ProjectionMatrixFormVol : public MatrixFormVol<Scalar>
      {
      public:
        ProjectionMatrixFormVol(int i, int j, ProjNormType projNormType) : MatrixFormVol<Scalar>(i, j)
        {
          this->projNormType = projNormType;
        }

        Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
          Geom<double> *e, Func<Scalar> **ext) const
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
          Geom<Hermes::Ord> *e, Func<Ord> **ext) const
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

        MatrixFormVol<Scalar>* clone() const
        {
          return new ProjectionMatrixFormVol(*this);
        }

      private:
        ProjNormType projNormType;

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        static SolFunctionDomain h1_projection_biform(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *u,
          Func<TestFunctionDomain> *v, Geom<TestFunctionDomain> *e, Func<SolFunctionDomain> **ext)
        {
          SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++)
            result += wt[i] * (u->val[i] * v->val[i] + u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
          return result;
        }

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        static SolFunctionDomain h1_semi_projection_biform(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *u,
          Func<TestFunctionDomain> *v, Geom<TestFunctionDomain> *e, Func<SolFunctionDomain> **ext)
        {
          SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++)
            result += wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
          return result;
        }

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        static SolFunctionDomain l2_projection_biform(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *u,
          Func<TestFunctionDomain> *v, Geom<TestFunctionDomain> *e, Func<SolFunctionDomain> **ext)
        {
          SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++)
            result += wt[i] * (u->val[i] * v->val[i]);
          return result;
        }

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        static SolFunctionDomain hcurl_projection_biform(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *u,
          Func<TestFunctionDomain> *v, Geom<TestFunctionDomain> *e, Func<SolFunctionDomain> **ext)
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
          Func<TestFunctionDomain> *v, Geom<TestFunctionDomain> *e, Func<SolFunctionDomain> **ext)
        {
          SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++) {
            result += wt[i] * (u->div[i] * conj(v->div[i]));
            result += wt[i] * (u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
          }
          return result;
        }
      };

      /// Residual.
      class ProjectionVectorFormVol : public VectorFormVol<Scalar>
      {
      public:
        ProjectionVectorFormVol(int i, ProjNormType projNormType) : VectorFormVol<Scalar>(i)
        {
          this->projNormType = projNormType;
        }

        Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
          Geom<double> *e, Func<Scalar> **ext) const
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
          Geom<Hermes::Ord> *e, Func<Ord> **ext) const
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

        VectorFormVol<Scalar>* clone() const
        {
          return new ProjectionVectorFormVol(*this);
        }

      private:
        ProjNormType projNormType;

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        SolFunctionDomain h1_projection_residual(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *v,
          Geom<TestFunctionDomain> *e, Func<SolFunctionDomain> **ext) const
        {
          SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++)
            result += wt[i] * ((ext[0]->val[i]) * v->val[i]
          + (ext[0]->dx[i]) * v->dx[i]
          + (ext[0]->dy[i]) * v->dy[i]);
          return result;
        }

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        SolFunctionDomain h1_semi_projection_residual(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *v,
          Geom<TestFunctionDomain> *e, Func<SolFunctionDomain> **ext) const
        {
          SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++)
            result += wt[i] * ((ext[0]->dx[i]) * v->dx[i]
          + (ext[0]->dy[i]) * v->dy[i]);
          return result;
        }

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        SolFunctionDomain l2_projection_residual(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *v,
          Geom<TestFunctionDomain> *e, Func<SolFunctionDomain> **ext) const
        {
          SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++)
            result += wt[i] * (ext[0]->val[i]) * v->val[i];
          return result;
        }

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        SolFunctionDomain hcurl_projection_residual(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *v,
          Geom<TestFunctionDomain> *e, Func<SolFunctionDomain> **ext) const
        {
          SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++) {
            result += wt[i] * (ext[0]->curl[i]) * conj(v->curl[i]);
            result += wt[i] * ((ext[0]->val0[i]) * conj(v->val0[i])
              + (ext[0]->val1[i]) * conj(v->val1[i]));
          }

          return result;
        }

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        SolFunctionDomain hdiv_projection_residual(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *v,
          Geom<TestFunctionDomain> *e, Func<SolFunctionDomain> **ext) const
        {
          SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++) {
            result += wt[i] * (ext[0]->div[i]) * conj(v->div[i]);
            result += wt[i] * ((ext[0]->val0[i]) * conj(v->val0[i])
              + (ext[0]->val1[i]) * conj(v->val1[i]));
          }

          return result;
        }
      };

      int ndof;
    };
  }
}
#endif
