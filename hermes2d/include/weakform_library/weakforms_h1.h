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

#ifndef __H2D_H1_WEAK_FORMS_H
#define __H2D_H1_WEAK_FORMS_H

#include "../integrals/h1.h"
#include "../weakform/weakform.h"
#include "../spline.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace WeakFormsH1
    {
      /* Default volumetric matrix form \int_{area} const_coeff * function_coeff(x, y) * u * v \bfx
      const_coeff... constant number
      function_coeff... (generally nonconstant) function of x, y
      */

      template<typename Scalar>
      class HERMES_API DefaultMatrixFormVol : public MatrixFormVol<Scalar>
      {
      public:
        DefaultMatrixFormVol(int i, int j, std::string area = HERMES_ANY,
          Hermes2DFunction<Scalar>* coeff = HERMES_ONE,
          SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

        DefaultMatrixFormVol(int i, int j, Hermes::vector<std::string> areas,
          Hermes2DFunction<Scalar>* coeff = HERMES_ONE,
          SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

        ~DefaultMatrixFormVol();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
          Geom<double> *e, ExtData<Scalar> *ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u,
          Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;

        virtual MatrixFormVol<Scalar>* clone();

      private:

        Hermes2DFunction<Scalar>* coeff;
        GeomType gt;
      };

      /* Default volumetric matrix form \int_{area} const_coeff * spline_coeff'(u_ext[0]) u \nabla u_ext[0] \cdot \nabla v
      + const_coeff * spline_coeff(u_ext[0]) * \nabla u \cdot \nabla v d\bfx
      const_coeff... constant number
      spline_coeff... nonconstant parameter given by cubic spline
      */

      template<typename Scalar>
      class HERMES_API DefaultJacobianDiffusion : public MatrixFormVol<Scalar>
      {
      public:
        DefaultJacobianDiffusion(int i, int j, std::string area = HERMES_ANY, Hermes1DFunction<Scalar>* coeff = HERMES_ONE,
          SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

        DefaultJacobianDiffusion(int i, int j, Hermes::vector<std::string> areas, Hermes1DFunction<Scalar>* coeff = HERMES_ONE,
          SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

        ~DefaultJacobianDiffusion();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
          Func<double> *v, Geom<double> *e, ExtData<Scalar> *ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u, Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;

        virtual MatrixFormVol<Scalar>* clone();

      private:
        int idx_j;

        Hermes1DFunction<Scalar>* coeff;
        GeomType gt;
      };

       template<typename Scalar>
      class HERMES_API DefaultMatrixFormDiffusion : public MatrixFormVol<Scalar>
      {
      public:
        DefaultMatrixFormDiffusion(int i, int j, std::string area = HERMES_ANY, Hermes1DFunction<Scalar>* coeff = HERMES_ONE,
          SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

        DefaultMatrixFormDiffusion(int i, int j, Hermes::vector<std::string> areas, Hermes1DFunction<Scalar>* coeff = HERMES_ONE,
          SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

        ~DefaultMatrixFormDiffusion();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
          Func<double> *v, Geom<double> *e, ExtData<Scalar> *ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u, Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;

        virtual MatrixFormVol<Scalar>* clone();

      private:
        int idx_j;

        Hermes1DFunction<Scalar>* coeff;
        GeomType gt;
      };

      /* Default volumetric matrix form
      \int_{area} spline_coeff1`(u_ext[0]) * u * u_ext[0]->dx * v
      + spline_coeff1(u_ext[0]) * u->dx * v
      + spline_coeff2`(u_ext[0]) * u * u_ext[0]->dy * v
      + spline_coeff2(u_ext[0]) * u->dy * v d\bfx.
      spline_coeff1, spline_coeff2... non-constant parameters given by cubic splines
      */

      template<typename Scalar>
      class HERMES_API DefaultJacobianAdvection : public MatrixFormVol<Scalar>
      {
      public:
        DefaultJacobianAdvection(int i, int j, std::string area = HERMES_ANY,
          Hermes1DFunction<Scalar>* coeff_1 = HERMES_ONE, Hermes1DFunction<Scalar>* coeff_2 = HERMES_ONE, GeomType gt = HERMES_PLANAR);

        DefaultJacobianAdvection(int i, int j, Hermes::vector<std::string> areas,
          Hermes1DFunction<Scalar>* coeff_1 = HERMES_ONE, Hermes1DFunction<Scalar>* coeff_2 = HERMES_ONE,
          GeomType gt = HERMES_PLANAR);

        ~DefaultJacobianAdvection();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
          Func<double> *v, Geom<double> *e, ExtData<Scalar> *ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u, Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;

        virtual MatrixFormVol<Scalar>* clone();

      private:
        int idx_j;
        Hermes1DFunction<Scalar>* coeff1, *coeff2;
        GeomType gt;
      };

      /* Default volumetric vector form \int_{area} const_coeff * function_coeff(x, y) * v d\bfx
      const_coeff... constant number
      function_coeff... (generally nonconstant) function of x, y
      */

      template<typename Scalar>
      class HERMES_API DefaultVectorFormVol : public VectorFormVol<Scalar>
      {
      public:
        DefaultVectorFormVol(int i, std::string area = HERMES_ANY, Hermes2DFunction<Scalar>* coeff = HERMES_ONE,
          GeomType gt = HERMES_PLANAR);

        DefaultVectorFormVol(int i, Hermes::vector<std::string> areas, Hermes2DFunction<Scalar>* coeff = HERMES_ONE,
          GeomType gt = HERMES_PLANAR);

        ~DefaultVectorFormVol();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
          Geom<double> *e, ExtData<Scalar> *ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;

        virtual VectorFormVol<Scalar>* clone();

      private:

        Hermes2DFunction<Scalar>* coeff;
        GeomType gt;
      };

      /* Default volumetric vector form \int_{area} const_coeff * function_coeff(x, y) * u_ext[0] * v d\bfx
      const_coeff... constant number
      function_coeff... (generally nonconstant) function of x, y
      */

      template<typename Scalar>
      class HERMES_API DefaultResidualVol : public VectorFormVol<Scalar>
      {
      public:
        DefaultResidualVol(int i, std::string area = HERMES_ANY, Hermes2DFunction<Scalar>* coeff = HERMES_ONE,
          GeomType gt = HERMES_PLANAR);
        DefaultResidualVol(int i, Hermes::vector<std::string> areas, Hermes2DFunction<Scalar>* coeff = HERMES_ONE,
          GeomType gt = HERMES_PLANAR);

        ~DefaultResidualVol();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
          Geom<double> *e, ExtData<Scalar> *ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;

        virtual VectorFormVol<Scalar>* clone();

      private:
        int idx_i;

        Hermes2DFunction<Scalar>* coeff;
        GeomType gt;
      };

      /* Default volumetric vector form \int_{area} const_coeff * spline_coeff(u_ext[0]) *
      \nabla u_ext[0] \cdot \nabla v d\bfx
      const_coeff... constant number
      spline_coeff... non-constant parameter given by a cubic spline
      */

      template<typename Scalar>
      class HERMES_API DefaultResidualDiffusion : public VectorFormVol<Scalar>
      {
      public:
        DefaultResidualDiffusion(int i, std::string area = HERMES_ANY, Hermes1DFunction<Scalar>* coeff = HERMES_ONE,
          GeomType gt = HERMES_PLANAR);

        DefaultResidualDiffusion(int i, Hermes::vector<std::string> areas, Hermes1DFunction<Scalar>* coeff = HERMES_ONE,
          GeomType gt = HERMES_PLANAR);

        ~DefaultResidualDiffusion();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
          Geom<double> *e, ExtData<Scalar> *ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;

        virtual VectorFormVol<Scalar>* clone();

      private:
        int idx_i;

        Hermes1DFunction<Scalar>* coeff;
        GeomType gt;
      };

      /* Default volumetric vector form \int_{area} spline_coeff1(u_ext[0]) * u->dx * v->val
      + spline_coeff2(u_ext[0]) * u->dy * v->val d\bfx
      spline_coeff1, spline_coeff2... non-constant parameters given by cubic splines
      */

      template<typename Scalar>
      class HERMES_API DefaultResidualAdvection : public VectorFormVol<Scalar>
      {
      public:
        DefaultResidualAdvection(int i, std::string area = HERMES_ANY,
          Hermes1DFunction<Scalar>* coeff_1 = HERMES_ONE, Hermes1DFunction<Scalar>* coeff_2 = HERMES_ONE,
          GeomType gt = HERMES_PLANAR);
        DefaultResidualAdvection(int i, Hermes::vector<std::string> areas,
          Hermes1DFunction<Scalar>* coeff_1 = HERMES_ONE, Hermes1DFunction<Scalar>* coeff_2 = HERMES_ONE, GeomType gt = HERMES_PLANAR);

        ~DefaultResidualAdvection();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
          Geom<double> *e, ExtData<Scalar> *ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;

        virtual VectorFormVol<Scalar>* clone();

      private:
        int idx_i;
        Hermes1DFunction<Scalar>* coeff1, *coeff2;
        GeomType gt;
      };

      /* Default surface matrix form \int_{area} const_coeff * function_coeff(x, y) * u * v dS
      const_coeff... constant number
      function_coeff... (generally nonconstant) function of x, y
      */

      template<typename Scalar>
      class HERMES_API DefaultMatrixFormSurf : public MatrixFormSurf<Scalar>
      {
      public:
        DefaultMatrixFormSurf(int i, int j, std::string area = HERMES_ANY,
          Hermes2DFunction<Scalar>* coeff = HERMES_ONE,
          GeomType gt = HERMES_PLANAR);

        DefaultMatrixFormSurf(int i, int j, Hermes::vector<std::string> areas,
          Hermes2DFunction<Scalar>* coeff = HERMES_ONE,
          GeomType gt = HERMES_PLANAR);

        ~DefaultMatrixFormSurf();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
          Geom<double> *e, ExtData<Scalar> *ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u,
          Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;

        virtual MatrixFormSurf<Scalar>* clone();

      private:

        Hermes2DFunction<Scalar>* coeff;
        GeomType gt;
      };

      /* Default surface matrix form \int_{area} const_coeff * spline_coeff'(u_ext[0]) * u_ext[0] * u * v
      + const_coeff * spline_coeff(u_ext[0]) * u * v dS
      spline_coeff... non-constant parameter given by a spline
      */

      template<typename Scalar>
      class HERMES_API DefaultJacobianFormSurf : public MatrixFormSurf<Scalar>
      {
      public:
        DefaultJacobianFormSurf(int i, int j, std::string area = HERMES_ANY, Hermes1DFunction<Scalar>* coeff = HERMES_ONE,
          GeomType gt = HERMES_PLANAR);
        DefaultJacobianFormSurf(int i, int j, Hermes::vector<std::string> areas, Hermes1DFunction<Scalar>* coeff = HERMES_ONE,
          GeomType gt = HERMES_PLANAR);

        ~DefaultJacobianFormSurf();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
          Geom<double> *e, ExtData<Scalar> *ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u,
          Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;

        virtual MatrixFormSurf<Scalar>* clone();

      private:
        int idx_j;

        Hermes1DFunction<Scalar>* coeff;
        GeomType gt;
      };

      /* Default surface vector form \int_{area} const_coeff * function_coeff(x, y) * v dS
      const_coeff... constant number
      function_coeff... (generally nonconstant) function of x, y
      */

      template<typename Scalar>
      class HERMES_API DefaultVectorFormSurf : public VectorFormSurf<Scalar>
      {
      public:
        DefaultVectorFormSurf(int i, std::string area = HERMES_ANY, Hermes2DFunction<Scalar>* coeff = HERMES_ONE,
          GeomType gt = HERMES_PLANAR);
        DefaultVectorFormSurf(int i, Hermes::vector<std::string> areas, Hermes2DFunction<Scalar>* coeff = HERMES_ONE,
          GeomType gt = HERMES_PLANAR);

        ~DefaultVectorFormSurf();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
          Geom<double> *e, ExtData<Scalar> *ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;

        virtual VectorFormSurf<Scalar>* clone();

      private:

        Hermes2DFunction<Scalar>* coeff;
        GeomType gt;
      };

      /* Default surface vector form \int_{area} const_coeff * function_coeff(x, y) * u_ext[0] v dS
      const_coeff... constant number
      function_coeff... (generally nonconstant) function of x, y
      */

      template<typename Scalar>
      class HERMES_API DefaultResidualSurf : public VectorFormSurf<Scalar>
      {
      public:
        DefaultResidualSurf(int i, std::string area = HERMES_ANY, Hermes2DFunction<Scalar>* coeff = HERMES_ONE,
          GeomType gt = HERMES_PLANAR);
        DefaultResidualSurf(int i, Hermes::vector<std::string> areas, Hermes2DFunction<Scalar>* coeff = HERMES_ONE,
          GeomType gt = HERMES_PLANAR);

        ~DefaultResidualSurf();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
          Geom<double> *e, ExtData<Scalar> *ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;

        virtual VectorFormSurf<Scalar>* clone();

      private:
        int idx_i;

        Hermes2DFunction<Scalar>* coeff;
        GeomType gt;
      };

      /* Default weak form for the Laplace equation -div(const_coeff spline_coeff(u) grad u) = 0. */

      template<typename Scalar>
      class HERMES_API DefaultWeakFormLaplace : public WeakForm<Scalar>
      {
      public:
        DefaultWeakFormLaplace(std::string area = HERMES_ANY, Hermes1DFunction<Scalar>* coeff = HERMES_ONE,
          GeomType gt = HERMES_PLANAR);
      };


      /* Default weak form for the Poisson equation -div(const_coeff spline_coeff(u) grad u) - rhs = 0. */

      template<typename Scalar>
      class HERMES_API DefaultWeakFormPoisson : public WeakForm<Scalar>
      {
      public:
        DefaultWeakFormPoisson();

        DefaultWeakFormPoisson(std::string area,
          Hermes1DFunction<Scalar>* coeff,
          Hermes2DFunction<Scalar>* f,
          GeomType gt = HERMES_PLANAR);
      };
    };
  }
}
#endif
