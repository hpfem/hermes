// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#include "weakforms_h1_const.h"
#include "api2d.h"
namespace Hermes
{
  namespace Hermes2D
  {
    namespace ConstantWeakFormsH1
    {
      template<typename Scalar>
      ConstantMatrixFormVol<Scalar>::ConstantMatrixFormVol
        (int i, int j, std::string area) : MatrixFormVol<Scalar>(i, j)
      {
        this->set_area(area);
        this->init_tables();
      }

      template<typename Scalar>
      ConstantMatrixFormVol<Scalar>::ConstantMatrixFormVol
        (int i, int j, Hermes::vector<std::string> areas) : MatrixFormVol<Scalar>(i, j)
      {
        this->set_areas(areas);
        this->init_tables();
      }

      template<typename Scalar>
      void ConstantMatrixFormVol<Scalar>::init_tables()
      {
        // Settings of precalculated values.
        this->set_h1_h1_const_tables(HERMES_MODE_TRIANGLE, "MatrixFormVolTriangle.h1h1");
        this->set_h1_h1_const_tables(HERMES_MODE_QUAD, "MatrixFormVolQuad.h1h1");
        this->set_l2_l2_const_tables(HERMES_MODE_TRIANGLE, "MatrixFormVolTriangle.l2l2");
        this->set_l2_l2_const_tables(HERMES_MODE_QUAD, "MatrixFormVolQuad.l2l2");

        /// \todo Cross-tables (h1 <-> l2)
        /// \todo Hcurl, Hdiv
      }

      template<typename Scalar>
      ConstantMatrixFormVol<Scalar>::~ConstantMatrixFormVol()
      {
      };

      template<typename Scalar>
      Scalar ConstantMatrixFormVol<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        for (int i = 0; i < n; i++)
          result += wt[i] * u->val[i] * v->val[i];
        return result;
      }

      template<typename Scalar>
      Ord ConstantMatrixFormVol<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
        Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        for (int i = 0; i < n; i++)
          result += wt[i] * u->val[i] * v->val[i];
        return result;
      }

      template<typename Scalar>
      MatrixFormVol<Scalar>* ConstantMatrixFormVol<Scalar>::clone() const
      {
        /// \todo Check that this copies the tables data.
        ConstantMatrixFormVol<Scalar>* toReturn = new ConstantMatrixFormVol<Scalar>(*this);
        toReturn->has_precalculated_tables = false;
        return toReturn;
      }

      template<typename Scalar>
      ConstantMatrixFormDx<Scalar>::ConstantMatrixFormDx
        (int i, int j, std::string area) : MatrixFormVol<Scalar>(i, j)
      {
        this->set_area(area);
        this->init_tables();
      }

      template<typename Scalar>
      ConstantMatrixFormDx<Scalar>::ConstantMatrixFormDx
        (int i, int j, Hermes::vector<std::string> areas) : MatrixFormVol<Scalar>(i, j)
      {
        this->set_areas(areas);
        this->init_tables();
      }

      template<typename Scalar>
      void ConstantMatrixFormDx<Scalar>::init_tables()
      {
        // Settings of precalculated values.
        this->set_h1_h1_const_tables(HERMES_MODE_TRIANGLE, "MatrixFormDxTriangle.h1h1");
        this->set_h1_h1_const_tables(HERMES_MODE_QUAD, "MatrixFormDxQuad.h1h1");
        this->set_l2_l2_const_tables(HERMES_MODE_TRIANGLE, "MatrixFormDxTriangle.l2l2");
        this->set_l2_l2_const_tables(HERMES_MODE_QUAD, "MatrixFormDxQuad.l2l2");

        /// \todo Cross-tables (h1 <-> l2)
        /// \todo Hcurl, Hdiv
      }

      template<typename Scalar>
      ConstantMatrixFormDx<Scalar>::~ConstantMatrixFormDx()
      {
      };

      template<typename Scalar>
      Scalar ConstantMatrixFormDx<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        for (int i = 0; i < n; i++)
          result += wt[i] * u->dx[i] * v->dx[i];
        return result;
      }

      template<typename Scalar>
      Ord ConstantMatrixFormDx<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
        Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        for (int i = 0; i < n; i++)
          result += wt[i] * u->dx[i] * v->dx[i];
        return result;
      }

      template<typename Scalar>
      MatrixFormVol<Scalar>* ConstantMatrixFormDx<Scalar>::clone() const
      {
        /// \todo Check that this copies the tables data.
        ConstantMatrixFormDx<Scalar>* toReturn = new ConstantMatrixFormDx<Scalar>(*this);
        toReturn->has_precalculated_tables = false;
        return toReturn;
      }

      template<typename Scalar>
      ConstantMatrixFormDy<Scalar>::ConstantMatrixFormDy
        (int i, int j, std::string area) : MatrixFormVol<Scalar>(i, j)
      {
        this->set_area(area);
        this->init_tables();
      }

      template<typename Scalar>
      ConstantMatrixFormDy<Scalar>::ConstantMatrixFormDy
        (int i, int j, Hermes::vector<std::string> areas) : MatrixFormVol<Scalar>(i, j)
      {
        this->set_areas(areas);
        this->init_tables();
      }

      template<typename Scalar>
      void ConstantMatrixFormDy<Scalar>::init_tables()
      {
        // Settings of precalculated values.
        this->set_h1_h1_const_tables(HERMES_MODE_TRIANGLE, "MatrixFormDyTriangle.h1h1");
        this->set_h1_h1_const_tables(HERMES_MODE_QUAD, "MatrixFormDyQuad.h1h1");
        this->set_l2_l2_const_tables(HERMES_MODE_TRIANGLE, "MatrixFormDyTriangle.l2l2");
        this->set_l2_l2_const_tables(HERMES_MODE_QUAD, "MatrixFormDyQuad.l2l2");

        /// \todo Cross-tables (h1 <-> l2)
        /// \todo Hcurl, Hdiv
      }

      template<typename Scalar>
      ConstantMatrixFormDy<Scalar>::~ConstantMatrixFormDy()
      {
      };

      template<typename Scalar>
      Scalar ConstantMatrixFormDy<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        for (int i = 0; i < n; i++)
          result += wt[i] * u->dy[i] * v->dy[i];
        return result;
      }

      template<typename Scalar>
      Ord ConstantMatrixFormDy<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
        Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        for (int i = 0; i < n; i++)
          result += wt[i] * u->dy[i] * v->dy[i];
        return result;
      }

      template<typename Scalar>
      MatrixFormVol<Scalar>* ConstantMatrixFormDy<Scalar>::clone() const
      {
        /// \todo Check that this copies the tables data.
        ConstantMatrixFormDy<Scalar>* toReturn = new ConstantMatrixFormDy<Scalar>(*this);
        toReturn->has_precalculated_tables = false;
        return toReturn;
      }

      template<typename Scalar>
      ConstantMatrixFormDuDxValV<Scalar>::ConstantMatrixFormDuDxValV
        (int i, int j, std::string area) : MatrixFormVol<Scalar>(i, j)
      {
        this->set_area(area);
        this->init_tables();
      }

      template<typename Scalar>
      ConstantMatrixFormDuDxValV<Scalar>::ConstantMatrixFormDuDxValV
        (int i, int j, Hermes::vector<std::string> areas) : MatrixFormVol<Scalar>(i, j)
      {
        this->set_areas(areas);
        this->init_tables();
      }

      template<typename Scalar>
      void ConstantMatrixFormDuDxValV<Scalar>::init_tables()
      {
        // Settings of precalculated values.
        this->set_h1_h1_const_tables(HERMES_MODE_TRIANGLE, "MatrixFormDuDxValVTriangle.h1h1");
        this->set_h1_h1_const_tables(HERMES_MODE_QUAD, "MatrixFormDuDxValVQuad.h1h1");
        this->set_l2_l2_const_tables(HERMES_MODE_TRIANGLE, "MatrixFormDuDxValVTriangle.l2l2");
        this->set_l2_l2_const_tables(HERMES_MODE_QUAD, "MatrixFormDuDxValVQuad.l2l2");

        /// \todo Cross-tables (h1 <-> l2)
        /// \todo Hcurl, Hdiv
      }

      template<typename Scalar>
      ConstantMatrixFormDuDxValV<Scalar>::~ConstantMatrixFormDuDxValV()
      {
      };

      template<typename Scalar>
      Scalar ConstantMatrixFormDuDxValV<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        for (int i = 0; i < n; i++)
          result += wt[i] * u->dx[i] * v->val[i];
        return result;
      }

      template<typename Scalar>
      Ord ConstantMatrixFormDuDxValV<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
        Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        for (int i = 0; i < n; i++)
          result += wt[i] * u->dx[i] * v->val[i];
        return result;
      }

      template<typename Scalar>
      MatrixFormVol<Scalar>* ConstantMatrixFormDuDxValV<Scalar>::clone() const
      {
        /// \todo Check that this copies the tables data.
        ConstantMatrixFormDuDxValV<Scalar>* toReturn = new ConstantMatrixFormDuDxValV<Scalar>(*this);
        toReturn->has_precalculated_tables = false;
        return toReturn;
      }

      template<typename Scalar>
      ConstantMatrixFormDuDyValV<Scalar>::ConstantMatrixFormDuDyValV
        (int i, int j, std::string area) : MatrixFormVol<Scalar>(i, j)
      {
        this->set_area(area);
        this->init_tables();
      }

      template<typename Scalar>
      ConstantMatrixFormDuDyValV<Scalar>::ConstantMatrixFormDuDyValV
        (int i, int j, Hermes::vector<std::string> areas) : MatrixFormVol<Scalar>(i, j)
      {
        this->set_areas(areas);
        this->init_tables();
      }

      template<typename Scalar>
      void ConstantMatrixFormDuDyValV<Scalar>::init_tables()
      {
        // Settings of precalculated values.
        this->set_h1_h1_const_tables(HERMES_MODE_TRIANGLE, "MatrixFormDuDyValVTriangle.h1h1");
        this->set_h1_h1_const_tables(HERMES_MODE_QUAD, "MatrixFormDuDyValVQuad.h1h1");
        this->set_l2_l2_const_tables(HERMES_MODE_TRIANGLE, "MatrixFormDuDyValVTriangle.l2l2");
        this->set_l2_l2_const_tables(HERMES_MODE_QUAD, "MatrixFormDuDyValVQuad.l2l2");

        /// \todo Cross-tables (h1 <-> l2)
        /// \todo Hcurl, Hdiv
      }

      template<typename Scalar>
      ConstantMatrixFormDuDyValV<Scalar>::~ConstantMatrixFormDuDyValV()
      {
      };

      template<typename Scalar>
      Scalar ConstantMatrixFormDuDyValV<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        for (int i = 0; i < n; i++)
          result += wt[i] * u->dy[i] * v->val[i];
        return result;
      }

      template<typename Scalar>
      Ord ConstantMatrixFormDuDyValV<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
        Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        for (int i = 0; i < n; i++)
          result += wt[i] * u->dy[i] * v->val[i];
        return result;
      }

      template<typename Scalar>
      MatrixFormVol<Scalar>* ConstantMatrixFormDuDyValV<Scalar>::clone() const
      {
        /// \todo Check that this copies the tables data.
        ConstantMatrixFormDuDyValV<Scalar>* toReturn = new ConstantMatrixFormDuDyValV<Scalar>(*this);
        toReturn->has_precalculated_tables = false;
        return toReturn;
      }

      template<typename Scalar>
      ConstantMatrixFormValUDvDx<Scalar>::ConstantMatrixFormValUDvDx
        (int i, int j, std::string area) : MatrixFormVol<Scalar>(i, j)
      {
        this->set_area(area);
        this->init_tables();
      }

      template<typename Scalar>
      ConstantMatrixFormValUDvDx<Scalar>::ConstantMatrixFormValUDvDx
        (int i, int j, Hermes::vector<std::string> areas) : MatrixFormVol<Scalar>(i, j)
      {
        this->set_areas(areas);
        this->init_tables();
      }

      template<typename Scalar>
      void ConstantMatrixFormValUDvDx<Scalar>::init_tables()
      {
        // Settings of precalculated values.
        this->set_h1_h1_const_tables(HERMES_MODE_TRIANGLE, "MatrixFormValUDvDxTriangle.h1h1");
        this->set_h1_h1_const_tables(HERMES_MODE_QUAD, "MatrixFormValUDvDxQuad.h1h1");
        this->set_l2_l2_const_tables(HERMES_MODE_TRIANGLE, "MatrixFormValUDvDxTriangle.l2l2");
        this->set_l2_l2_const_tables(HERMES_MODE_QUAD, "MatrixFormValUDvDxQuad.l2l2");

        /// \todo Cross-tables (h1 <-> l2)
        /// \todo Hcurl, Hdiv
      }

      template<typename Scalar>
      ConstantMatrixFormValUDvDx<Scalar>::~ConstantMatrixFormValUDvDx()
      {
      };

      template<typename Scalar>
      Scalar ConstantMatrixFormValUDvDx<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        for (int i = 0; i < n; i++)
          result += wt[i] * u->val[i] * v->dx[i];
        return result;
      }

      template<typename Scalar>
      Ord ConstantMatrixFormValUDvDx<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
        Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        for (int i = 0; i < n; i++)
          result += wt[i] * u->val[i] * v->dx[i];
        return result;
      }

      template<typename Scalar>
      MatrixFormVol<Scalar>* ConstantMatrixFormValUDvDx<Scalar>::clone() const
      {
        /// \todo Check that this copies the tables data.
        ConstantMatrixFormValUDvDx<Scalar>* toReturn = new ConstantMatrixFormValUDvDx<Scalar>(*this);
        toReturn->has_precalculated_tables = false;
        return toReturn;
      }

      template<typename Scalar>
      ConstantMatrixFormValUDvDy<Scalar>::ConstantMatrixFormValUDvDy
        (int i, int j, std::string area) : MatrixFormVol<Scalar>(i, j)
      {
        this->set_area(area);
        this->init_tables();
      }

      template<typename Scalar>
      ConstantMatrixFormValUDvDy<Scalar>::ConstantMatrixFormValUDvDy
        (int i, int j, Hermes::vector<std::string> areas) : MatrixFormVol<Scalar>(i, j)
      {
        this->set_areas(areas);
        this->init_tables();
      }

      template<typename Scalar>
      void ConstantMatrixFormValUDvDy<Scalar>::init_tables()
      {
        // Settings of precalculated values.
        this->set_h1_h1_const_tables(HERMES_MODE_TRIANGLE, "MatrixFormValUDvDyTriangle.h1h1");
        this->set_h1_h1_const_tables(HERMES_MODE_QUAD, "MatrixFormValUDvDyQuad.h1h1");
        this->set_l2_l2_const_tables(HERMES_MODE_TRIANGLE, "MatrixFormValUDvDyTriangle.l2l2");
        this->set_l2_l2_const_tables(HERMES_MODE_QUAD, "MatrixFormValUDvDyQuad.l2l2");

        /// \todo Cross-tables (h1 <-> l2)
        /// \todo Hcurl, Hdiv
      }

      template<typename Scalar>
      ConstantMatrixFormValUDvDy<Scalar>::~ConstantMatrixFormValUDvDy()
      {
      };

      template<typename Scalar>
      Scalar ConstantMatrixFormValUDvDy<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        for (int i = 0; i < n; i++)
          result += wt[i] * u->val[i] * v->dy[i];
        return result;
      }

      template<typename Scalar>
      Ord ConstantMatrixFormValUDvDy<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
        Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        for (int i = 0; i < n; i++)
          result += wt[i] * u->val[i] * v->dy[i];
        return result;
      }

      template<typename Scalar>
      MatrixFormVol<Scalar>* ConstantMatrixFormValUDvDy<Scalar>::clone() const
      {
        /// \todo Check that this copies the tables data.
        ConstantMatrixFormValUDvDy<Scalar>* toReturn = new ConstantMatrixFormValUDvDy<Scalar>(*this);
        toReturn->has_precalculated_tables = false;
        return toReturn;
      }

      template<typename Scalar>
      ConstantVectorFormVol<Scalar>::ConstantVectorFormVol
        (int i, std::string area) : VectorFormVol<Scalar>(i)
      {
        this->set_area(area);
        this->init_tables();
      }

      template<typename Scalar>
      ConstantVectorFormVol<Scalar>::ConstantVectorFormVol
        (int i, Hermes::vector<std::string> areas) : VectorFormVol<Scalar>(i)
      {
        this->set_areas(areas);
        this->init_tables();
      }
        
      template<typename Scalar>
      void ConstantVectorFormVol<Scalar>::init_tables()
      {
        // Settings of precalculated values.
        this->set_h1_const_tables(HERMES_MODE_TRIANGLE, "VectorFormVolTriangle.h1");
        this->set_h1_const_tables(HERMES_MODE_QUAD, "VectorFormVolQuad.h1");
        this->set_l2_const_tables(HERMES_MODE_TRIANGLE, "VectorFormVolTriangle.l2");
        this->set_l2_const_tables(HERMES_MODE_QUAD, "VectorFormVolQuad.l2");
        /// \todo Hcurl, Hdiv
      }

      template<typename Scalar>
      ConstantVectorFormVol<Scalar>::~ConstantVectorFormVol()
      {
      };

      template<typename Scalar>
      Scalar ConstantVectorFormVol<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        for (int i = 0; i < n; i++)
          result += wt[i] * v->val[i];
        return result;
      }

      template<typename Scalar>
      Ord ConstantVectorFormVol<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
        Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        for (int i = 0; i < n; i++)
          result += wt[i] * v->val[i];
        return result;
      }

      template<typename Scalar>
      VectorFormVol<Scalar>* ConstantVectorFormVol<Scalar>::clone() const
      {
        /// \todo Check that this copies the tables data.
        ConstantVectorFormVol<Scalar>* toReturn = new ConstantVectorFormVol<Scalar>(*this);
        toReturn->has_precalculated_tables = false;
        return toReturn;
      }

      template<typename Scalar>
      ConstantVectorFormDx<Scalar>::ConstantVectorFormDx
        (int i, std::string area) : VectorFormVol<Scalar>(i)
      {
        this->set_area(area);
        this->init_tables();
      }

      template<typename Scalar>
      ConstantVectorFormDx<Scalar>::ConstantVectorFormDx
        (int i, Hermes::vector<std::string> areas) : VectorFormVol<Scalar>(i)
      {
        this->set_areas(areas);
        this->init_tables();
      }

      template<typename Scalar>
      void ConstantVectorFormDx<Scalar>::init_tables()
      {
        // Settings of precalculated values.
        this->set_h1_const_tables(HERMES_MODE_TRIANGLE, "VectorFormDxTriangle.h1");
        this->set_h1_const_tables(HERMES_MODE_QUAD, "VectorFormDxQuad.h1");
        this->set_l2_const_tables(HERMES_MODE_TRIANGLE, "VectorFormDxTriangle.l2");
        this->set_l2_const_tables(HERMES_MODE_QUAD, "VectorFormDxQuad.l2");
        /// \todo Hcurl, Hdiv
      }

      template<typename Scalar>
      ConstantVectorFormDx<Scalar>::~ConstantVectorFormDx()
      {
      };

      template<typename Scalar>
      Scalar ConstantVectorFormDx<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        for (int i = 0; i < n; i++)
          result += wt[i] * v->dx[i];
        return result;
      }

      template<typename Scalar>
      Ord ConstantVectorFormDx<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
        Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        for (int i = 0; i < n; i++)
          result += wt[i] * v->dx[i];
        return result;
      }

      template<typename Scalar>
      VectorFormVol<Scalar>* ConstantVectorFormDx<Scalar>::clone() const
      {
        /// \todo Check that this copies the tables data.
        ConstantVectorFormDx<Scalar>* toReturn = new ConstantVectorFormDx<Scalar>(*this);
        toReturn->has_precalculated_tables = false;
        return toReturn;
      }

      template<typename Scalar>
      ConstantVectorFormDy<Scalar>::ConstantVectorFormDy
        (int i, std::string area) : VectorFormVol<Scalar>(i)
      {
        this->set_area(area);
        this->init_tables();
      }

      template<typename Scalar>
      ConstantVectorFormDy<Scalar>::ConstantVectorFormDy
        (int i, Hermes::vector<std::string> areas) : VectorFormVol<Scalar>(i)
      {
        this->set_areas(areas);
        this->init_tables();
      }

      template<typename Scalar>
      void ConstantVectorFormDy<Scalar>::init_tables()
      {
        // Settings of precalculated values.
        this->set_h1_const_tables(HERMES_MODE_TRIANGLE, "VectorFormDyTriangle.h1");
        this->set_h1_const_tables(HERMES_MODE_QUAD, "VectorFormDyQuad.h1");
        this->set_l2_const_tables(HERMES_MODE_TRIANGLE, "VectorFormDyTriangle.l2");
        this->set_l2_const_tables(HERMES_MODE_QUAD, "VectorFormDyQuad.l2");
        /// \todo Hcurl, Hdiv
      }

      template<typename Scalar>
      ConstantVectorFormDy<Scalar>::~ConstantVectorFormDy()
      {
      };

      template<typename Scalar>
      Scalar ConstantVectorFormDy<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        for (int i = 0; i < n; i++)
          result += wt[i] * v->dy[i];
        return result;
      }

      template<typename Scalar>
      Ord ConstantVectorFormDy<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
        Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        for (int i = 0; i < n; i++)
          result += wt[i] * v->dy[i];
        return result;
      }

      template<typename Scalar>
      VectorFormVol<Scalar>* ConstantVectorFormDy<Scalar>::clone() const
      {
        /// \todo Check that this copies the tables data.
        ConstantVectorFormDy<Scalar>* toReturn = new ConstantVectorFormDy<Scalar>(*this);
        toReturn->has_precalculated_tables = false;
        return toReturn;
      }

      template class HERMES_API ConstantMatrixFormVol<double>;
      template class HERMES_API ConstantMatrixFormVol<std::complex<double> >;
      template class HERMES_API ConstantMatrixFormDx<double>;
      template class HERMES_API ConstantMatrixFormDx<std::complex<double> >;
      template class HERMES_API ConstantMatrixFormDy<double>;
      template class HERMES_API ConstantMatrixFormDy<std::complex<double> >;
      template class HERMES_API ConstantVectorFormVol<double>;
      template class HERMES_API ConstantVectorFormVol<std::complex<double> >;
      template class HERMES_API ConstantVectorFormDx<double>;
      template class HERMES_API ConstantVectorFormDx<std::complex<double> >;
      template class HERMES_API ConstantVectorFormDy<double>;
      template class HERMES_API ConstantVectorFormDy<std::complex<double> >;
    };
  }
}
