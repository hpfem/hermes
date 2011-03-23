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

#ifndef __H2D_SAMPLE_WEAK_FORMS_H
#define __H2D_SAMPLE_WEAK_FORMS_H

#include "../integrals/integrals_h1.h"

/* Default volumetric matrix form \int_{area} coeff \nabla u \cdot \nabla v d\bfx 
   coeff... constant number
*/

class DefaultMatrixFormVolConst : public WeakForm::MatrixFormVol
{
public:
  DefaultMatrixFormVolConst(int i, int j, double coeff = 1.0) 
        : WeakForm::MatrixFormVol(i, j, HERMES_SYM), coeff(coeff) { }
  DefaultMatrixFormVolConst(int i, int j, std::string area, double coeff = 1.0) 
        : WeakForm::MatrixFormVol(i, j, HERMES_SYM, area), coeff(coeff) { }

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                     Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
    return coeff * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
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

/* Default volumetric matrix form \int_{area} coeff u v d\bfx 
   coeff... constant number
*/

class DefaultMatrixFormVolMassConst : public WeakForm::MatrixFormVol
{
public:
  DefaultMatrixFormVolMassConst(int i, int j, double coeff = 1.0) 
        : WeakForm::MatrixFormVol(i, j, HERMES_SYM), coeff(coeff) { }
  DefaultMatrixFormVolMassConst(int i, int j, std::string area, double coeff = 1.0) 
        : WeakForm::MatrixFormVol(i, j, HERMES_SYM, area), coeff(coeff) { }

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                     Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
    return coeff * int_u_v<Real, Scalar>(n, wt, u, v);
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

/* Default volumetric vector form \int_{area} coeff v d\bfx 
   coeff... constant number
*/

class DefaultVectorFormVolConst : public WeakForm::VectorFormVol
{
public:
  DefaultVectorFormVolConst(int i, double coeff) 
               : WeakForm::VectorFormVol(i), coeff(coeff) { }
  DefaultVectorFormVolConst(int i, std::string area, double coeff) 
               : WeakForm::VectorFormVol(i, area), coeff(coeff) { }

  virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                       Geom<double> *e, ExtData<scalar> *ext) {
    return coeff * int_v<scalar, scalar>(n, wt, v);
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
          Geom<Ord> *e, ExtData<Ord> *ext) {
    return int_v<Ord, Ord>(n, wt, v);
  }

private:
  double coeff;
};

/* Default surface matrix form \int_{area} coeff u v dS
   coeff... constant number
*/

class DefaultMatrixFormSurfConst : public WeakForm::MatrixFormSurf
{
public:
  DefaultMatrixFormSurfConst(int i, int j, double coeff) 
        : WeakForm::MatrixFormSurf(i, j), coeff(coeff) { }
  DefaultMatrixFormSurfConst(int i, int j, std::string area, double coeff) 
        : WeakForm::MatrixFormSurf(i, j, area), coeff(coeff) { }

  template<typename Real, typename Scalar>
  Scalar matrix_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                          Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
    return coeff * int_u_v<Real, Scalar>(n, wt, u, v);
  }

  scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
               Geom<double> *e, ExtData<scalar> *ext) {
    return matrix_form_surf<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
    return matrix_form_surf<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
  }

  private:
    double coeff;
};


/* Default surface vector form \int_{area} coeff v dS
   coeff... constant number
*/

class DefaultVectorFormSurfConst : public WeakForm::VectorFormSurf
{
public:
  DefaultVectorFormSurfConst(int i, double coeff) 
         : WeakForm::VectorFormSurf(i), coeff(coeff) { }
  DefaultVectorFormSurfConst(int i, std::string area, double coeff) 
         : WeakForm::VectorFormSurf(i, area), coeff(coeff) { }

  template<typename Real, typename Scalar>
  Scalar vector_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
    return coeff * int_v<Real, Scalar>(n, wt, v);
  }

  scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
    return vector_form_surf<scalar, scalar>(n, wt, u_ext, v, e, ext);
  }

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
    return vector_form_surf<Ord, Ord>(n, wt, u_ext, v, e, ext);
  }

private:
  double coeff;
};


#endif
