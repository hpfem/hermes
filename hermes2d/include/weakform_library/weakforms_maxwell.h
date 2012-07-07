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

#ifndef __H2D_MAXWELL_WEAK_FORMS_H
#define __H2D_MAXWELL_WEAK_FORMS_H

#include "../integrals/h1.h"
#include "../weakform/weakform.h"
#include "../spline.h"
namespace Hermes
{
  namespace Hermes2D
  {
    namespace WeakFormsMaxwell {
      /* Default volumetric matrix form \int_{area} coeff_spline(u_ext[0]) \curl u \curl v d\bfx
      spline_coeff... nonconstant parameter given by cubic spline
      */

      template<typename Scalar>
      class HERMES_API DefaultJacobianMagnetostatics : public MatrixFormVol<Scalar>
      {
      public:
        DefaultJacobianMagnetostatics(int i, int j, std::string area = HERMES_ANY, Scalar const_coeff = 1.0,
          CubicSpline* c_spline = HERMES_DEFAULT_SPLINE, SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR,
          int order_increase = 3);
        DefaultJacobianMagnetostatics(int i, int j, Hermes::vector<std::string> areas,
          Scalar const_coeff, CubicSpline* c_spline = HERMES_DEFAULT_SPLINE,
          SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR, int order_increase = 3);

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
          Func<double> *v, Geom<double> *e, ExtData<Scalar> *ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u, Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;

        // This is to make the form usable in rk_time_step_newton().
        virtual MatrixFormVol<Scalar>* clone();

      private:
        int idx_j;
        Scalar const_coeff;
        CubicSpline* spline_coeff;
        GeomType gt;
        int order_increase;
      };

      template<typename Scalar>
      class HERMES_API DefaultResidualMagnetostatics : public VectorFormVol<Scalar>
      {
      public:
        DefaultResidualMagnetostatics(int i, std::string area = HERMES_ANY, Scalar const_coeff = 1.0,
          CubicSpline* c_spline = HERMES_DEFAULT_SPLINE,
          GeomType gt = HERMES_PLANAR,
          int order_increase = 3);
        DefaultResidualMagnetostatics(int i, Hermes::vector<std::string> areas, Scalar const_coeff = 1.0,
          CubicSpline* c_spline = HERMES_DEFAULT_SPLINE,
          GeomType gt = HERMES_PLANAR, int order_increase = 3);

        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
          Geom<double> *e, ExtData<Scalar> *ext) const;

        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
          Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;

        // This is to make the form usable in rk_time_step_newton().
        virtual VectorFormVol<Scalar>* clone();

      private:
        int idx_i;
        Scalar const_coeff;
        CubicSpline* spline_coeff;
        GeomType gt;
        int order_increase;
      };
    };
  }
}
#endif