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

#ifndef __H2D_HCURL_WEAK_FORMS_H
#define __H2D_HCURL_WEAK_FORMS_H

#include "../integrals/hcurl.h"
#include "../weakform/weakform.h"
#include "../spline.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace WeakFormsHcurl
    {
      /* Default volumetric matrix form \int_{area} const_coeff * function_coeff(x, y) * E \cdot F d\bfx
      const_coeff... constant number
      function_coeff... (generally nonconstant) function of x, y
      */

      template<typename Scalar>
      class HERMES_API DefaultMatrixFormVol : public MatrixFormVol<Scalar>
      {
      public:
        DefaultMatrixFormVol<Scalar>(int i, int j, std::string area = HERMES_ANY,
          Scalar const_coeff = 1.0, Hermes2DFunction<Scalar>* f_coeff = HERMES_DEFAULT_FUNCTION,
          SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

        DefaultMatrixFormVol<Scalar>(int i, int j, Hermes::vector<std::string> areas,
          Scalar const_coeff = 1.0, Hermes2DFunction<Scalar>* f_coeff = HERMES_DEFAULT_FUNCTION,
          SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

        ~DefaultMatrixFormVol<Scalar>();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
          Func<double> *v, Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u, Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual MatrixFormVol<Scalar>* clone() const;

      private:
        Scalar const_coeff;
        Hermes2DFunction<Scalar>* function_coeff;
        GeomType gt;
      };

      /* FIXME
      Default volumetric matrix form \int_{area} const_coeff \curl E \curl F d\bfx
      coeff... constant number
      */

      template<typename Scalar>
      class HERMES_API DefaultJacobianCurlCurl : public MatrixFormVol<Scalar>
      {
      public:
        DefaultJacobianCurlCurl(int i, int j, std::string area = HERMES_ANY, Scalar const_coeff = 1.0,
          CubicSpline* c_spline = HERMES_DEFAULT_SPLINE,
          SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

        DefaultJacobianCurlCurl(int i, int j, Hermes::vector<std::string> areas, Scalar const_coeff = 1.0,
          CubicSpline* c_spline = HERMES_DEFAULT_SPLINE,
          SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR);

        ~DefaultJacobianCurlCurl();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
          Func<double> *v, Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u, Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual MatrixFormVol<Scalar>* clone() const;

      private:
        int idx_j;
        Scalar const_coeff;
        CubicSpline* spline_coeff;
        GeomType gt;
      };

      /* FIXME
      Default volumetric vector form \int_{area} (coeff0, coeff1) \cdot E d\bfx
      coeff0, coeff1... constant numbers
      */

      template<typename Scalar>
      class HERMES_API DefaultVectorFormVol : public VectorFormVol<Scalar>
      {
      public:
        DefaultVectorFormVol<Scalar>(int i, std::string area = HERMES_ANY,
          Scalar const_coeff0 = 1.0, Scalar const_coeff1 = 1.0,
          Hermes2DFunction<Scalar>* f_coeff0 = HERMES_DEFAULT_FUNCTION,
          Hermes2DFunction<Scalar>* f_coeff1 = HERMES_DEFAULT_FUNCTION,
          GeomType gt = HERMES_PLANAR);

        DefaultVectorFormVol<Scalar>(int i, Hermes::vector<std::string> areas,
          Scalar const_coeff0 = 1.0, Scalar const_coeff1 = 1.0,
          Hermes2DFunction<Scalar>* f_coeff0 = HERMES_DEFAULT_FUNCTION,
          Hermes2DFunction<Scalar>* f_coeff1 = HERMES_DEFAULT_FUNCTION,
          GeomType gt = HERMES_PLANAR);

        ~DefaultVectorFormVol<Scalar>();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
          Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual VectorFormVol<Scalar>* clone() const;

      private:
        Scalar const_coeff0, const_coeff1;
        Hermes2DFunction<Scalar>* function_coeff0, *function_coeff1;
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
        DefaultResidualVol(int i, std::string area = HERMES_ANY, Scalar const_coeff = 1.0,
          Hermes2DFunction<Scalar>* f_coeff = HERMES_DEFAULT_FUNCTION,
          GeomType gt = HERMES_PLANAR);

        DefaultResidualVol(int i, Hermes::vector<std::string> areas, Scalar const_coeff = 1.0,
          Hermes2DFunction<Scalar>* f_coeff = HERMES_DEFAULT_FUNCTION,
          GeomType gt = HERMES_PLANAR);

        ~DefaultResidualVol();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
          Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual VectorFormVol<Scalar>* clone() const;

      private:
        int idx_i;
        Scalar const_coeff;
        Hermes2DFunction<Scalar>* function_coeff;
        GeomType gt;
      };

      /* FIXME
      Default volumetric vector form \int_{area} const_coeff * spline_coeff(u_ext[0]) *
      \nabla u_ext[0] \cdot \nabla v d\bfx
      const_coeff... constant number
      spline_coeff... non-constant parameter given by a cubic spline
      */

      template<typename Scalar>
      class HERMES_API DefaultResidualCurlCurl : public VectorFormVol<Scalar>
      {
      public:
        DefaultResidualCurlCurl(int i, std::string area = HERMES_ANY, Scalar const_coeff = 1.0,
          CubicSpline* c_spline = HERMES_DEFAULT_SPLINE,
          GeomType gt = HERMES_PLANAR);

        DefaultResidualCurlCurl(int i, Hermes::vector<std::string> areas, Scalar const_coeff = 1.0,
          CubicSpline* c_spline = HERMES_DEFAULT_SPLINE,
          GeomType gt = HERMES_PLANAR);

        ~DefaultResidualCurlCurl();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
          Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual VectorFormVol<Scalar>* clone() const;

      private:
        int idx_i;
        Scalar const_coeff;
        CubicSpline* spline_coeff;
        GeomType gt;
      };

      /* FIXME
      Default surface matrix form \int_{area} coeff e tau f tau dS
      coeff... constant number
      */

      template<typename Scalar>
      class HERMES_API DefaultMatrixFormSurf : public MatrixFormSurf<Scalar>
      {
      public:
        DefaultMatrixFormSurf<Scalar>(int i, int j, std::string area = HERMES_ANY,
          Scalar const_coeff = 1.0, Hermes2DFunction<Scalar>* f_coeff = HERMES_DEFAULT_FUNCTION,
          GeomType gt = HERMES_PLANAR);

        DefaultMatrixFormSurf<Scalar>(int i, int j, Hermes::vector<std::string> areas,
          Scalar const_coeff = 1.0, Hermes2DFunction<Scalar>* f_coeff = HERMES_DEFAULT_FUNCTION,
          GeomType gt = HERMES_PLANAR);

        ~DefaultMatrixFormSurf<Scalar>();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
          Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u,
          Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual MatrixFormSurf<Scalar>* clone() const;

      private:
        Scalar const_coeff;
        Hermes2DFunction<Scalar>* function_coeff;
        GeomType gt;
      };

      /* FIXME
      Default surface vector form \int_{area} const_coeff * function_coeff(x, y) * v dS
      const_coeff... constant number
      function_coeff... (generally nonconstant) function of x, y
      */

      template<typename Scalar>
      class HERMES_API DefaultVectorFormSurf : public VectorFormSurf<Scalar>
      {
      public:
        DefaultVectorFormSurf<Scalar>(int i, std::string area = HERMES_ANY, Scalar const_coeff = 1.0,
          Hermes2DFunction<Scalar>* f_coeff = HERMES_DEFAULT_FUNCTION,
          GeomType gt = HERMES_PLANAR);

        DefaultVectorFormSurf<Scalar>(int i, Hermes::vector<std::string> areas, Scalar const_coeff = 1.0,
          Hermes2DFunction<Scalar>* f_coeff = HERMES_DEFAULT_FUNCTION,
          GeomType gt = HERMES_PLANAR);

        ~DefaultVectorFormSurf<Scalar>();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
          Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual VectorFormSurf<Scalar>* clone() const;

      private:
        Scalar const_coeff;
        Hermes2DFunction<Scalar>* function_coeff;
        GeomType gt;
      };

      /* FIXME
      Default surface residual form \int_{area} coeff u_ext[0] tau f tau dS
      coeff... constant number
      */

      template<typename Scalar>
      class HERMES_API DefaultResidualSurf : public VectorFormSurf<Scalar>
      {
      public:
        DefaultResidualSurf(int i, std::string area = HERMES_ANY,
          Scalar const_coeff = 1.0, Hermes2DFunction<Scalar>* f_coeff = HERMES_DEFAULT_FUNCTION,
          GeomType gt = HERMES_PLANAR);

        DefaultResidualSurf(int i, Hermes::vector<std::string> areas,
          Scalar const_coeff = 1.0, Hermes2DFunction<Scalar>* f_coeff = HERMES_DEFAULT_FUNCTION,
          GeomType gt = HERMES_PLANAR);

        ~DefaultResidualSurf();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
          Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[],
          Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual VectorFormSurf<Scalar>* clone() const;

      private:
        Scalar const_coeff;
        Hermes2DFunction<Scalar>* function_coeff;
        GeomType gt;
      };
    }
  }
}
#endif