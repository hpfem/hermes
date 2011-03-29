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

#include "../integrals/integrals_hcurl.h"

/* Default volumetric matrix form \int_{area} coeff \curl E \curl F d\bfx 
   coeff... constant number
*/

namespace WeakFormsHcurl {
  namespace VolumetricMatrixForms {
    class DefaultMatrixFormCurlCurl : public WeakForm::MatrixFormVol
    {
    public:
      DefaultMatrixFormCurlCurl(int i, int j, double coeff = 1.0, SymFlag sym = HERMES_SYM) 
            : WeakForm::MatrixFormVol(i, j, sym), coeff(coeff) { }
      DefaultMatrixFormCurlCurl(int i, int j, std::string area, double coeff = 1.0, SymFlag sym = HERMES_SYM) 
            : WeakForm::MatrixFormVol(i, j, sym, area), coeff(coeff) { }

      template<typename Real, typename Scalar>
      Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                         Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
        return coeff * int_curl_e_curl_f<Real, Scalar>(n, wt, u, v);
      }

      scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                   Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
        return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
      }

      Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
              Geom<Ord> *e, ExtData<Ord> *ext) {
        return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
      }

      private:
        double coeff;
    };

    /* Default volumetric matrix form \int_{area} coeff E \cdot F d\bfx 
       coeff... constant number
    */

    class DefaultMatrixFormMass : public WeakForm::MatrixFormVol
    {
    public:
      DefaultMatrixFormMass(int i, int j, double coeff = 1.0, SymFlag sym = HERMES_SYM) 
            : WeakForm::MatrixFormVol(i, j, sym), coeff(coeff) { }
      DefaultMatrixFormMass(int i, int j, std::string area, double coeff = 1.0, SymFlag sym = HERMES_SYM) 
            : WeakForm::MatrixFormVol(i, j, sym, area), coeff(coeff) { }

      template<typename Real, typename Scalar>
      Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                         Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
        return coeff * int_e_f<Real, Scalar>(n, wt, u, v);
      }

      scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                   Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
        return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
      }

      Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
              Geom<Ord> *e, ExtData<Ord> *ext) {
        return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
      }

      private:
        double coeff;
    };
  }

  namespace VolumetricVectorForms {
      /* Default volumetric vector form \int_{area} (coeff0, coeff1) \cdot E d\bfx 
         coeff... constant number
      */

      class DefaultVectorFormConst : public WeakForm::VectorFormVol
      {
      public:
      DefaultVectorFormConst(int i, double coeff0, coeff1) 
	     : WeakForm::VectorFormVol(i), coeff0(coeff0), coeff1(coeff1) { }
      DefaultVectorFormConst(int i, std::string area, double coeff0, double coeff1) 
	       : WeakForm::VectorFormVol(i, area), coeff0(coeff0), coeff1(coeff1) { }

        virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                             Geom<double> *e, ExtData<scalar> *ext) {
          return coeff0 * int_v0<scalar, scalar>(n, wt, v) + 
                 coeff1 * int_v1<scalar, scalar>(n, wt, v);
        }

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                Geom<Ord> *e, ExtData<Ord> *ext) {
          return coeff0 * int_v0<Ord, Ord>(n, wt, v) + 
                 coeff1 * int_v1<Ord, Ord>(n, wt, v);
        }

      private:
        double coeff0, coeff1;
      };
  }
}
#endif
