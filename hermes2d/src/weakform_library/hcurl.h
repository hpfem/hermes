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

#include "../integrals/integrals_hcurl.h"

namespace WeakFormsHcurl {

  namespace VolumetricMatrixForms {

    /* Default volumetric matrix form \int_{area} coeff \curl E \curl F d\bfx
     coeff... constant number
    */

    class DefaultLinearCurlCurl : public WeakForm::MatrixFormVol
    {
    public:
      DefaultLinearCurlCurl(int i, int j, scalar coeff = 1.0, SymFlag sym = HERMES_SYM)
            : WeakForm::MatrixFormVol(i, j, sym), coeff(coeff) { }
      DefaultLinearCurlCurl(int i, int j, std::string area, scalar coeff = 1.0, SymFlag sym = HERMES_SYM)
            : WeakForm::MatrixFormVol(i, j, sym, area), coeff(coeff) { }

      template<typename Real, typename Scalar>
      Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                         Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
        return coeff * int_curl_e_curl_f<Real, Scalar>(n, wt, u, v);
      }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                           Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
        return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                      Geom<Ord> *e, ExtData<Ord> *ext) const {
        return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
      }

      virtual WeakForm::MatrixFormVol* clone() {
        return new DefaultLinearCurlCurl(*this);
      }

      private:
        scalar coeff;
    };

    /* Default volumetric matrix form \int_{area} coeff E \cdot F d\bfx
       coeff... constant number
    */

    class DefaultLinearMass : public WeakForm::MatrixFormVol
    {
    public:
      DefaultLinearMass(int i, int j, scalar coeff = 1.0, SymFlag sym = HERMES_SYM)
            : WeakForm::MatrixFormVol(i, j, sym), coeff(coeff) { }
      DefaultLinearMass(int i, int j, std::string area, scalar coeff = 1.0, SymFlag sym = HERMES_SYM)
            : WeakForm::MatrixFormVol(i, j, sym, area), coeff(coeff) { }

      template<typename Real, typename Scalar>
      Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                         Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
        return coeff * int_e_f<Real, Scalar>(n, wt, u, v);
      }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                           Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
        return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
              Geom<Ord> *e, ExtData<Ord> *ext) const {
        return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
      }

      virtual WeakForm::MatrixFormVol* clone() {
        return new DefaultLinearMass(*this);
      }

      private:
        scalar coeff;
    };
  }

  namespace SurfaceMatrixForms {

    /* Default surface matrix form \int_{area} coeff e tau f tau dS
       coeff... constant number
    */

    class DefaultMatrixFormSurf : public WeakForm::MatrixFormSurf
    {
    public:
      DefaultMatrixFormSurf(int i, int j, scalar coeff)
        : WeakForm::MatrixFormSurf(i, j), coeff(coeff) { }
      DefaultMatrixFormSurf(int i, int j, std::string area, scalar coeff)
        : WeakForm::MatrixFormSurf(i, j, area), coeff(coeff) { }

      template<typename Real, typename Scalar>
      Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                         Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
        return coeff * int_e_tau_f_tau<Real, Scalar>(n, wt, u, v, e);
      }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                           Geom<double> *e, ExtData<scalar> *ext) const {
        return matrix_form(n, wt, u_ext, u, v, e, ext);
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                      Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
        return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
      }

      virtual WeakForm::MatrixFormSurf* clone() {
        return new DefaultMatrixFormSurf(*this);
      }

      private:
        scalar coeff;
    };
  }

  namespace VolumetricVectorForms {

      /* Default volumetric vector form \int_{area} (coeff0, coeff1) \cdot E d\bfx
         coeff... constant number
      */

      class DefaultVectorFormConst : public WeakForm::VectorFormVol
      {
      public:
      DefaultVectorFormConst(int i, scalar coeff0, scalar coeff1)
       : WeakForm::VectorFormVol(i), coeff0(coeff0), coeff1(coeff1) { }
      DefaultVectorFormConst(int i, std::string area, scalar coeff0, scalar coeff1)
         : WeakForm::VectorFormVol(i, area), coeff0(coeff0), coeff1(coeff1) { }

        virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                             Geom<double> *e, ExtData<scalar> *ext) const {
          scalar int_v0 = 0, int_v1 = 0;
          for (int i = 0; i < n; i++) int_v0 += wt[i] * v->val0[i];
          for (int i = 0; i < n; i++) int_v1 += wt[i] * v->val1[i];
          return coeff0 * int_v0 + coeff1 * int_v1;
        }

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                Geom<Ord> *e, ExtData<Ord> *ext) const {
          Ord int_v0 = 0, int_v1 = 0;
          for (int i = 0; i < n; i++) int_v0 += wt[i] * v->val0[i];
          for (int i = 0; i < n; i++) int_v1 += wt[i] * v->val1[i];
          return coeff0 * int_v0 + coeff1 * int_v1;
        }

      virtual WeakForm::VectorFormVol* clone() {
        return new DefaultVectorFormConst(*this);
      }

      private:
        scalar coeff0, coeff1;
      };
  }
}
#endif
