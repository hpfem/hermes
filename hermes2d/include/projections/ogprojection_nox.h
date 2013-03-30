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
/*! \file ogprojection_nox.h
\brief Orthogonal projection via NOX (matrix-free).
*/
#ifndef __H2D_OGPROJECTION_NOX_H
#define __H2D_OGPROJECTION_NOX_H

#include "../function/solution.h"
#include "../forms.h"
#include "../weakform/weakform.h"

//#include "epetra.h"
#if(defined HAVE_NOX && defined HAVE_EPETRA && defined HAVE_TEUCHOS)
#include <NOX.H>
#ifdef _POSIX_C_SOURCE
# undef _POSIX_C_SOURCE  // pyconfig.h included by NOX_Epetra defines it
#endif
#ifdef _XOPEN_SOURCE
# undef _XOPEN_SOURCE  // pyconfig.h included by NOX_Epetra defines it
#endif
#include <NOX_Epetra.H>

namespace Hermes
{
  namespace Hermes2D
  {
    /// @ingroup projections
    template<typename Scalar>
    class HERMES_API OGProjectionNOX : public Hermes::Mixins::Loggable
    {
    public:
      OGProjectionNOX();

      /// Main functionality is in the protected method project_internal().
      /// This is a wrapper that allows the user to specify his own projection form.
      void project_global(SpaceSharedPtr<Scalar> space,
          MatrixFormVol<Scalar>* custom_projection_jacobian,
          VectorFormVol<Scalar>* custom_projection_residual,
          Scalar* target_vec, double newton_tol = 1e-6, int newton_max_iter = 10);

      /**
       \fn  static void OGProjection::project_global(SpaceSharedPtr<Scalar> space,
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
      void project_global(SpaceSharedPtr<Scalar> space, MeshFunction<Scalar>* source_meshfn,
          Scalar* target_vec, ProjNormType proj_norm = HERMES_UNSET_NORM,
          double newton_tol = 1e-6, int newton_max_iter = 10);

      /// Wrapper that accepts MeshFunctionSharedPtr instead of the ordinary MeshFunction pointer.
      void project_global(SpaceSharedPtr<Scalar> space, MeshFunctionSharedPtr<Scalar> source_meshfn,
          Scalar* target_vec, ProjNormType proj_norm = HERMES_UNSET_NORM,
          double newton_tol = 1e-6, int newton_max_iter = 10);
          
      /// Wrapper that delivers a MeshFunctionSharedPtr instead of coefficient vector.
      void project_global(SpaceSharedPtr<Scalar> space,
          MeshFunctionSharedPtr<Scalar> source_sln, MeshFunctionSharedPtr<Scalar> target_sln,
          ProjNormType proj_norm = HERMES_UNSET_NORM,
          double newton_tol = 1e-6, int newton_max_iter = 10);

      /// Wrapper for multiple source MeshFunction pointers that delivers coefficient vector.
      void project_global(Hermes::vector<SpaceSharedPtr<Scalar> > spaces, Hermes::vector<MeshFunction<Scalar>* > source_meshfns,
          Scalar* target_vec, Hermes::vector<ProjNormType> proj_norms = Hermes::vector<ProjNormType>(),
          double newton_tol = 1e-6, int newton_max_iter = 10);

      /// Wrapper for multiple source MeshFunctionSharedPtrs that delivers coefficient vector.
      void project_global(Hermes::vector<SpaceSharedPtr<Scalar> > spaces, Hermes::vector<MeshFunctionSharedPtr<Scalar> > source_slns,
          Scalar* target_vec, Hermes::vector<ProjNormType> proj_norms = Hermes::vector<ProjNormType>(),
          double newton_tol = 1e-6, int newton_max_iter = 10);

      void project_global(Hermes::vector<SpaceSharedPtr<Scalar> > spaces,
          Hermes::vector<MeshFunctionSharedPtr<Scalar> > source_slns, Hermes::vector<MeshFunctionSharedPtr<Scalar> > target_slns,
          Hermes::vector<ProjNormType> proj_norms = Hermes::vector<ProjNormType>(), bool delete_old_mesh = false,
          double newton_tol = 1e-6, int newton_max_iter = 10);

    protected:
      /// Underlying function for global orthogonal projection.
      /// Not intended for the user. NOTE: the weak form here must be
      /// a special projection weak form, which is different from
      /// the weak form of the PDE. If you supply a weak form of the
      /// PDE, the PDE will just be solved.
      void project_internal(SpaceSharedPtr<Scalar> space, WeakForm<Scalar>* proj_wf, Scalar* target_vec,
        double newton_tol = 1e-6, int newton_max_iter = 10);

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

        MatrixFormVol<Scalar>* clone()
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

        VectorFormVol<Scalar>* clone()
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
            result += wt[i] * ((u_ext[this->i]->val[i] - ext[0]->val[i]) * v->val[i]
          + (u_ext[this->i]->dx[i] - ext[0]->dx[i]) * v->dx[i]
          + (u_ext[this->i]->dy[i] - ext[0]->dy[i]) * v->dy[i]);
          return result;
        }

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        SolFunctionDomain h1_semi_projection_residual(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *v,
          Geom<TestFunctionDomain> *e, Func<SolFunctionDomain> **ext) const
        {
          SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++)
            result += wt[i] * ((u_ext[this->i]->dx[i] - ext[0]->dx[i]) * v->dx[i]
          + (u_ext[this->i]->dy[i] - ext[0]->dy[i]) * v->dy[i]);
          return result;
        }

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        SolFunctionDomain l2_projection_residual(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *v,
          Geom<TestFunctionDomain> *e, Func<SolFunctionDomain> **ext) const
        {
          SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++)
            result += wt[i] * (u_ext[this->i]->val[i] - ext[0]->val[i]) * v->val[i];
          return result;
        }

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        SolFunctionDomain hcurl_projection_residual(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *v,
          Geom<TestFunctionDomain> *e, Func<SolFunctionDomain> **ext) const
        {
          SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++) {
            result += wt[i] * (u_ext[this->i]->curl[i] - ext[0]->curl[i]) * conj(v->curl[i]);
            result += wt[i] * ((u_ext[this->i]->val0[i] - ext[0]->val0[i]) * conj(v->val0[i])
              + (u_ext[this->i]->val1[i] - ext[0]->val1[i]) * conj(v->val1[i]));
          }

          return result;
        }

        template<typename TestFunctionDomain, typename SolFunctionDomain>
        SolFunctionDomain hdiv_projection_residual(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<TestFunctionDomain> *v,
          Geom<TestFunctionDomain> *e, Func<SolFunctionDomain> **ext) const
        {
          SolFunctionDomain result = SolFunctionDomain(0);
          for (int i = 0; i < n; i++) {
            result += wt[i] * (u_ext[this->i]->div[i] - ext[0]->div[i]) * conj(v->div[i]);
            result += wt[i] * ((u_ext[this->i]->val0[i] - ext[0]->val0[i]) * conj(v->val0[i])
              + (u_ext[this->i]->val1[i] - ext[0]->val1[i]) * conj(v->val1[i]));
          }

          return result;
        }
      };

      int ndof;
    };
  }
}
#endif
#endif