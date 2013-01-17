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

#ifndef __H2D_CONST_H1_WEAK_FORMS_H
#define __H2D_CONST_H1_WEAK_FORMS_H

#include "../integrals/h1.h"
#include "../weakform/weakform.h"
#include "../spline.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace ConstantWeakFormsH1
    {
      template<typename Scalar>
      class HERMES_API ConstantMatrixFormVol : public MatrixFormVol<Scalar>
      {
      public:
        ConstantMatrixFormVol(int i, int j, std::string area = HERMES_ANY);

        ConstantMatrixFormVol(int i, int j, Hermes::vector<std::string> areas);

        void init_tables();

        ~ConstantMatrixFormVol();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
          Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u,
          Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual MatrixFormVol<Scalar>* clone() const;
      };

      template<typename Scalar>
      class HERMES_API ConstantMatrixFormDx : public MatrixFormVol<Scalar>
      {
      public:
        ConstantMatrixFormDx(int i, int j, std::string area = HERMES_ANY);

        ConstantMatrixFormDx(int i, int j, Hermes::vector<std::string> areas);

        void init_tables();

        ~ConstantMatrixFormDx();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
          Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u,
          Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual MatrixFormVol<Scalar>* clone() const;
      };

      template<typename Scalar>
      class HERMES_API ConstantMatrixFormDy : public MatrixFormVol<Scalar>
      {
      public:
        ConstantMatrixFormDy(int i, int j, std::string area = HERMES_ANY);

        ConstantMatrixFormDy(int i, int j, Hermes::vector<std::string> areas);

        void init_tables();

        ~ConstantMatrixFormDy();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
          Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u,
          Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual MatrixFormVol<Scalar>* clone() const;
      };

      template<typename Scalar>
      class HERMES_API ConstantMatrixFormDuDxValV : public MatrixFormVol<Scalar>
      {
      public:
        ConstantMatrixFormDuDxValV(int i, int j, std::string area = HERMES_ANY);

        ConstantMatrixFormDuDxValV(int i, int j, Hermes::vector<std::string> areas);

        void init_tables();

        ~ConstantMatrixFormDuDxValV();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
          Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u,
          Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual MatrixFormVol<Scalar>* clone() const;
      };

      template<typename Scalar>
      class HERMES_API ConstantMatrixFormDuDyValV : public MatrixFormVol<Scalar>
      {
      public:
        ConstantMatrixFormDuDyValV(int i, int j, std::string area = HERMES_ANY);

        ConstantMatrixFormDuDyValV(int i, int j, Hermes::vector<std::string> areas);

        void init_tables();

        ~ConstantMatrixFormDuDyValV();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
          Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u,
          Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual MatrixFormVol<Scalar>* clone() const;
      };

      template<typename Scalar>
      class HERMES_API ConstantMatrixFormValUDvDx : public MatrixFormVol<Scalar>
      {
      public:
        ConstantMatrixFormValUDvDx(int i, int j, std::string area = HERMES_ANY);

        ConstantMatrixFormValUDvDx(int i, int j, Hermes::vector<std::string> areas);

        void init_tables();

        ~ConstantMatrixFormValUDvDx();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
          Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u,
          Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual MatrixFormVol<Scalar>* clone() const;
      };

      template<typename Scalar>
      class HERMES_API ConstantMatrixFormValUDvDy : public MatrixFormVol<Scalar>
      {
      public:
        ConstantMatrixFormValUDvDy(int i, int j, std::string area = HERMES_ANY);

        ConstantMatrixFormValUDvDy(int i, int j, Hermes::vector<std::string> areas);

        void init_tables();

        ~ConstantMatrixFormValUDvDy();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
          Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u,
          Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual MatrixFormVol<Scalar>* clone() const;
      };

      /* Default volumetric vector form \int_{area} const_coeff * function_coeff(x, y) * v d\bfx
      const_coeff... constant number
      function_coeff... (generally nonconstant) function of x, y
      */

      template<typename Scalar>
      class HERMES_API ConstantVectorFormVol : public VectorFormVol<Scalar>
      {
      public:
        ConstantVectorFormVol(int i, std::string area = HERMES_ANY);

        ConstantVectorFormVol(int i, Hermes::vector<std::string> areas);

        void init_tables();

        ~ConstantVectorFormVol();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
          Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual VectorFormVol<Scalar>* clone() const;
      };

      template<typename Scalar>
      class HERMES_API ConstantVectorFormDx : public VectorFormVol<Scalar>
      {
      public:
        ConstantVectorFormDx(int i, std::string area = HERMES_ANY);

        ConstantVectorFormDx(int i, Hermes::vector<std::string> areas);

        void init_tables();

        ~ConstantVectorFormDx();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
          Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual VectorFormVol<Scalar>* clone() const;
      };

      template<typename Scalar>
      class HERMES_API ConstantVectorFormDy : public VectorFormVol<Scalar>
      {
      public:
        ConstantVectorFormDy(int i, std::string area = HERMES_ANY);

        ConstantVectorFormDy(int i, Hermes::vector<std::string> areas);

        void init_tables();

        ~ConstantVectorFormDy();

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
          Geom<double> *e, Func<Scalar> **ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, Func<Ord> **ext) const;

        virtual VectorFormVol<Scalar>* clone() const;
      };
    }
  }
}
#endif
  
