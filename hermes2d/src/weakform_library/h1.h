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

#include "../integrals/integrals_h1.h"

namespace WeakFormsH1 {

  namespace VolumetricMatrixForms {

    /* Default volumetric matrix form \int_{area} coeff \nabla u \cdot \nabla v d\bfx
       coeff... constant number
    */

    class DefaultLinearDiffusion : public WeakForm::MatrixFormVol
    {
    public:
      DefaultLinearDiffusion(int i, int j, scalar coeff = 1.0,
                             SymFlag sym = HERMES_SYM, GeomType gt = HERMES_PLANAR)
            : WeakForm::MatrixFormVol(i, j, sym), coeff(coeff), gt(gt) { }
      DefaultLinearDiffusion(int i, int j, std::string area, scalar coeff = 1.0,
                             SymFlag sym = HERMES_SYM, GeomType gt = HERMES_PLANAR)
      : WeakForm::MatrixFormVol(i, j, sym, area), coeff(coeff), gt(gt) { }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                   Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
        scalar result = 0;
        if (gt == HERMES_PLANAR) result = int_grad_u_grad_v<double, scalar>(n, wt, u, v);
        else {
          if (gt == HERMES_AXISYM_X) result = int_y_grad_u_grad_v<double, scalar>(n, wt, u, v, e);
          else result = int_x_grad_u_grad_v<double, scalar>(n, wt, u, v, e);
        }
        return coeff * result;
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
              Geom<Ord> *e, ExtData<Ord> *ext) const {
        Ord result;
        if (gt == HERMES_PLANAR) result = int_grad_u_grad_v<Ord, Ord>(n, wt, u, v);
        else {
          if (gt == HERMES_AXISYM_X) result = int_y_grad_u_grad_v<Ord, Ord>(n, wt, u, v, e);
          else result = int_x_grad_u_grad_v<Ord, Ord>(n, wt, u, v, e);
        }
        return result;
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::MatrixFormVol* clone() {
        return new DefaultLinearDiffusion(*this);
      }

      private:
        scalar coeff;
        GeomType gt;
    };

    /* Default volumetric matrix form \int_{area} coeff_spline'(u_ext[0]) u \nabla u_ext[0] \cdot \nabla v
       + spline_coeff(u_ext[0]) * \nabla u \cdot \nabla v d\bfx
       coeff_spline... nonconstant parameter given by cubic spline
    */

    class DefaultJacobianNonlinearDiffusion : public WeakForm::MatrixFormVol
    {
    public:
      DefaultJacobianNonlinearDiffusion(int i, int j, CubicSpline* spline_coeff, scalar const_coeff = 1.0,
                                        SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR)
  : WeakForm::MatrixFormVol(i, j, sym), spline_coeff(spline_coeff), const_coeff(const_coeff), gt(gt) { }
      DefaultJacobianNonlinearDiffusion(int i, int j, std::string area, CubicSpline* spline_coeff, scalar const_coeff = 1.0,
                                        SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR)
  : WeakForm::MatrixFormVol(i, j, sym, area), spline_coeff(spline_coeff), const_coeff(const_coeff), gt(gt) { }

      template<typename Real, typename Scalar>
      Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                         Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
        Scalar result = 0;
        for (int i = 0; i < n; i++) {
          result += wt[i] * (const_coeff*spline_coeff->get_derivative(u_ext[0]->val[i]) * u->val[i] *
                 (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i])
                 + const_coeff*spline_coeff->get_value(u_ext[0]->val[i])
                                   * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
        }
        return result;
      }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                   Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
        return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
              Geom<Ord> *e, ExtData<Ord> *ext) const {
        return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::MatrixFormVol* clone() {
        return new DefaultJacobianNonlinearDiffusion(*this);
      }

      private:
        CubicSpline* spline_coeff;
        scalar const_coeff;
        GeomType gt;
    };

    /* Default volumetric matrix form \int_{area} coeff u v d\bfx
       coeff... constant number
    */

    class DefaultLinearMass : public WeakForm::MatrixFormVol
    {
    public:
      DefaultLinearMass(int i, int j, scalar coeff = 1.0,
                        SymFlag sym = HERMES_SYM, GeomType gt = HERMES_PLANAR)
  : WeakForm::MatrixFormVol(i, j, sym), coeff(coeff), gt(gt) { }
      DefaultLinearMass(int i, int j, std::string area, scalar coeff = 1.0,
                        SymFlag sym = HERMES_SYM, GeomType gt = HERMES_PLANAR)
  : WeakForm::MatrixFormVol(i, j, sym, area), coeff(coeff), gt(gt) { }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                   Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
        scalar result = 0;
        if (gt == HERMES_PLANAR) result = int_u_v<double, scalar>(n, wt, u, v);
        else if (gt == HERMES_AXISYM_X) result = int_y_u_v<double, scalar>(n, wt, u, v, e);
        else result = int_x_u_v<double, scalar>(n, wt, u, v, e);

        return coeff * result;
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
              Geom<Ord> *e, ExtData<Ord> *ext) const {
        Ord result = 0;
        if (gt == HERMES_PLANAR) result = int_u_v<Ord, Ord>(n, wt, u, v);
        else if (gt == HERMES_AXISYM_X) result = int_y_u_v<Ord, Ord>(n, wt, u, v, e);
        else result = int_x_u_v<Ord, Ord>(n, wt, u, v, e);

        return result;
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::MatrixFormVol* clone() {
        return new DefaultLinearMass(*this);
      }

      private:
        scalar coeff;
        GeomType gt;
    };

    /* Default volumetric matrix form \int_{area} coeff_spline(u_ext[0]) u v d\bfx
       spline_coeff... non-constant parameter given by a cubic spline
    */

    class DefaultJacobianNonlinearMass : public WeakForm::MatrixFormVol
    {
    public:
      DefaultJacobianNonlinearMass(int i, int j, CubicSpline* spline_coeff, scalar const_coeff = 1.0,
                                   SymFlag sym = HERMES_SYM, GeomType gt = HERMES_PLANAR)
  : WeakForm::MatrixFormVol(i, j, sym), spline_coeff(spline_coeff), const_coeff(const_coeff), gt(gt) { }
      DefaultJacobianNonlinearMass(int i, int j, std::string area,
                                   CubicSpline* spline_coeff, scalar const_coeff = 1.0, SymFlag sym = HERMES_SYM,
                                   GeomType gt = HERMES_PLANAR)
  : WeakForm::MatrixFormVol(i, j, sym, area), spline_coeff(spline_coeff), const_coeff(const_coeff), gt(gt) { }

      template<typename Real, typename Scalar>
      Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                         Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
        Scalar result = 0;
        for (int i = 0; i < n; i++) {
          result += wt[i] * (const_coeff*spline_coeff->get_value(u_ext[0]->val[i]) * (u->val[i] * v->val[i]));
        }
        return result;
      }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                   Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
        return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
              Geom<Ord> *e, ExtData<Ord> *ext) const {
        return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::MatrixFormVol* clone() {
        return new DefaultJacobianNonlinearMass(*this);
      }

      private:
        CubicSpline* spline_coeff;
        scalar const_coeff;
        GeomType gt;
    };

    /* Default volumetric matrix form \int_{area} (coeff1, coeff2) \cdot \nabla u vd\bfx
       coeff1, coeff2... constant numbers
    */

    class DefaultLinearAdvection : public WeakForm::MatrixFormVol
    {
    public:
      DefaultLinearAdvection(int i, int j, scalar coeff1, scalar coeff2, GeomType gt = HERMES_PLANAR)
       : WeakForm::MatrixFormVol(i, j, HERMES_NONSYM), coeff1(coeff1), coeff2(coeff2), gt(gt) { }
      DefaultLinearAdvection(int i, int j, std::string area, scalar coeff1, scalar coeff2, GeomType gt = HERMES_PLANAR)
       : WeakForm::MatrixFormVol(i, j, HERMES_NONSYM, area), coeff1(coeff1), coeff2(coeff2), gt(gt) { }

      template<typename Real, typename Scalar>
      Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                         Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
        return   coeff1 * int_dudx_v<Real, Scalar>(n, wt, u, v)
               + coeff2 * int_dudy_v<Real, Scalar>(n, wt, u, v);
      }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                   Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
        return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
              Geom<Ord> *e, ExtData<Ord> *ext) const {
        return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::MatrixFormVol* clone() {
        return new DefaultLinearAdvection(*this);
      }

      private:
      scalar coeff1, coeff2;
      GeomType gt;
    };

    /* Default volumetric matrix form
       \int_{area} spline_coeff1`(u_ext[0]) * u * u_ext[0]->dx * v
       + spline_coeff1(u_ext[0]) * u->dx * v
       + spline_coeff2`(u_ext[0]) * u * u_ext[0]->dy * v
       + spline_coeff2(u_ext[0]) * u->dy * v d\bfx.
       spline_coeff1, spline_coeff2... non-constant parameters given by cubic splines
    */

    class DefaultJacobianNonlinearAdvection : public WeakForm::MatrixFormVol
    {
    public:
     DefaultJacobianNonlinearAdvection(int i, int j, CubicSpline* spline_coeff1,
                                       CubicSpline* spline_coeff2, scalar const_coeff = 1.0, GeomType gt = HERMES_PLANAR)
       : WeakForm::MatrixFormVol(i, j, HERMES_NONSYM), spline_coeff1(spline_coeff1),
  spline_coeff2(spline_coeff2), const_coeff(const_coeff), gt(gt) { }
     DefaultJacobianNonlinearAdvection(int i, int j, std::string area,
                                       CubicSpline* spline_coeff1, CubicSpline* spline_coeff2,
                                       scalar const_coeff = 1.0, GeomType gt = HERMES_PLANAR)
       : WeakForm::MatrixFormVol(i, j, HERMES_NONSYM, area), spline_coeff1(spline_coeff1),
                           spline_coeff2(spline_coeff2), const_coeff(const_coeff), gt(gt) { }

      template<typename Real, typename Scalar>
      Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                         Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
        Scalar result = 0;
        for (int i = 0; i < n; i++) {
          result += wt[i] * (  const_coeff*spline_coeff1->get_derivative(u_ext[0]->val[i]) * u->val[i] * u_ext[0]->dx[i] * v->val[i]
                             + const_coeff*spline_coeff1->get_value(u_ext[0]->val[i]) * u->dx[i] * v->val[i]
                             + const_coeff*spline_coeff2->get_derivative(u_ext[0]->val[i]) * u->val[i] * u_ext[0]->dy[i] * v->val[i]
           + const_coeff*spline_coeff2->get_value(u_ext[0]->val[i]) * u->dy[i] * v->val[i]);
        }
        return result;
      }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                   Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
        return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
              Geom<Ord> *e, ExtData<Ord> *ext) const {
        return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::MatrixFormVol* clone() {
        return new DefaultJacobianNonlinearAdvection(*this);
      }

      private:
      CubicSpline* spline_coeff1, *spline_coeff2;
      scalar const_coeff;
      GeomType gt;
    };
  }

  namespace RightHandSides {
    // Generic class for non-constant right-hand side.
    class DefaultNonConstRightHandSide
    {
    public:
      DefaultNonConstRightHandSide() { };

      virtual scalar value(double x, double y) const = 0;
      virtual Ord ord(Ord x, Ord y) const = 0;
    };
  }

  namespace VolumetricVectorForms {

    /* Default volumetric vector form \int_{area} coeff v d\bfx
       coeff... constant number
    */

    class DefaultVectorFormConst : public WeakForm::VectorFormVol
    {
    public:
      DefaultVectorFormConst(int i, scalar coeff = 1.0, GeomType gt = HERMES_PLANAR)
             : WeakForm::VectorFormVol(i), coeff(coeff), gt(gt) { }
      DefaultVectorFormConst(int i, std::string area, scalar coeff = 1.0, 
                             GeomType gt = HERMES_PLANAR)
             : WeakForm::VectorFormVol(i, area), coeff(coeff), gt(gt) { }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                           Geom<double> *e, ExtData<scalar> *ext) const {
        if (gt == HERMES_PLANAR) return coeff * int_v<double>(n, wt, v);
        else {
          if (gt == HERMES_AXISYM_X) return coeff * int_y_v<double>(n, wt, v, e);
          else return coeff * int_x_v<double>(n, wt, v, e);
        }
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                      Geom<Ord> *e, ExtData<Ord> *ext) const {
        if (gt == HERMES_PLANAR) return int_v<Ord>(n, wt, v);
        else {
          if (gt == HERMES_AXISYM_X) return int_y_v<Ord>(n, wt, v, e);
          else return int_x_v<Ord>(n, wt, v, e);
        }
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::VectorFormVol* clone() {
        return new DefaultVectorFormConst(*this);
      }

      private:
        scalar coeff;
        GeomType gt;
    };

    /* Default volumetric vector form \int_{area} coeff
       u_ext[0] v d\bfx
       coeff... constant parameter
    */

    class DefaultResidualLinearMass : public WeakForm::VectorFormVol
    {
    public:
      DefaultResidualLinearMass(int i, scalar coeff, GeomType gt = HERMES_PLANAR)
  : WeakForm::VectorFormVol(i), coeff(coeff), gt(gt) { }
      DefaultResidualLinearMass(int i, std::string area, scalar coeff, GeomType gt = HERMES_PLANAR)
  : WeakForm::VectorFormVol(i, area), coeff(coeff), gt(gt) { }

      template<typename Real, typename Scalar>
      Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                         Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
        Scalar result = 0;
        Func<Scalar>* u_prev = u_ext[0];
        for (int i = 0; i < n; i++) {
          result += wt[i] * coeff * u_prev->val[i] * v->val[i];
        }
        return result;
      }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                           Geom<double> *e, ExtData<scalar> *ext) const {
        return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
              Geom<Ord> *e, ExtData<Ord> *ext) const {
        return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::VectorFormVol* clone() {
        return new DefaultResidualLinearMass(*this);
      }

      private:
        scalar coeff;
        GeomType gt;
    };

    /* Default volumetric vector form \int_{area} coeff
       \nabla u_ext[0] \cdot \nabla v d\bfx
       coeff... constant parameter
    */

    class DefaultResidualLinearDiffusion : public WeakForm::VectorFormVol
    {
    public:
      DefaultResidualLinearDiffusion(int i, scalar coeff = 1.0, GeomType gt = HERMES_PLANAR)
  : WeakForm::VectorFormVol(i), coeff(coeff), gt(gt) { }
      DefaultResidualLinearDiffusion(int i, std::string area, scalar coeff = 1.0, GeomType gt = HERMES_PLANAR)
  : WeakForm::VectorFormVol(i, area), coeff(coeff), gt(gt) { }

      template<typename Real, typename Scalar>
      Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                         Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
        Scalar result = 0;
        Func<Scalar>* u_prev = u_ext[0];
        for (int i = 0; i < n; i++) {
          result += wt[i] * coeff * (u_prev->dx[i] * v->dx[i] + u_prev->dy[i] * v->dy[i]);
        }
        return result;
      }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                           Geom<double> *e, ExtData<scalar> *ext) const {
        return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
              Geom<Ord> *e, ExtData<Ord> *ext) const {
        return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::VectorFormVol* clone() {
        return new DefaultResidualLinearDiffusion(*this);
      }

      private:
        scalar coeff;
        GeomType gt;
    };

    /* Default volumetric vector form \int_{area} spline_coeff(u_ext[0])
       \nabla u_ext[0] \cdot \nabla v d\bfx
       spline_coeff... non-constant parameter given by a cubic spline
    */

    class DefaultResidualNonlinearDiffusion : public WeakForm::VectorFormVol
    {
    public:
      DefaultResidualNonlinearDiffusion(int i, CubicSpline* spline_coeff, scalar const_coeff = 1.0, GeomType gt = HERMES_PLANAR)
  : WeakForm::VectorFormVol(i), spline_coeff(spline_coeff), const_coeff(const_coeff), gt(gt) { }
      DefaultResidualNonlinearDiffusion(int i, std::string area, CubicSpline* spline_coeff, scalar const_coeff = 1.0, GeomType gt = HERMES_PLANAR)
  : WeakForm::VectorFormVol(i, area), spline_coeff(spline_coeff), const_coeff(const_coeff), gt(gt) { }

      template<typename Real, typename Scalar>
      Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                         Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
        Scalar result = 0;
        Func<Scalar>* u_prev = u_ext[0];
        for (int i = 0; i < n; i++) {
          result += wt[i] * (const_coeff*spline_coeff->get_value(u_prev->val[i]) *
                             (u_prev->dx[i] * v->dx[i] + u_prev->dy[i] * v->dy[i]));
        }
        return result;
      }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                           Geom<double> *e, ExtData<scalar> *ext) const {
        return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
              Geom<Ord> *e, ExtData<Ord> *ext) const {
        return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::VectorFormVol* clone() {
        return new DefaultResidualNonlinearDiffusion(*this);
      }

      private:
        CubicSpline* spline_coeff;
        scalar const_coeff;
        GeomType gt;
    };

    /* Default volumetric vector form \int_{area} coeff1 * u->dx * v->val
       + coeff2 * u->dy * v->val d\bfx
       coeff1, coeff2... constant parameters
    */

    class DefaultResidualLinearAdvection : public WeakForm::VectorFormVol
    {
    public:
    DefaultResidualLinearAdvection(int i, scalar coeff1, scalar coeff2, GeomType gt = HERMES_PLANAR)
      : WeakForm::VectorFormVol(i), coeff1(coeff1), coeff2(coeff2), gt(gt) { }
      DefaultResidualLinearAdvection(int i, std::string area, scalar coeff1, scalar coeff2, GeomType gt = HERMES_PLANAR)
  : WeakForm::VectorFormVol(i, area), coeff1(coeff1), coeff2(coeff2), gt(gt) { }

      template<typename Real, typename Scalar>
      Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                         Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
        Scalar result = 0;
        Func<Scalar>* u_prev = u_ext[0];
        for (int i = 0; i < n; i++) {
          result += wt[i] * (coeff1 * (u_prev->dx[i] * v->val[i])
                             + coeff2 * (u_prev->dy[i] * v->val[i]));
        }
        return result;
      }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                           Geom<double> *e, ExtData<scalar> *ext) const {
        return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
              Geom<Ord> *e, ExtData<Ord> *ext) const {
        return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::VectorFormVol* clone() {
        return new DefaultResidualLinearAdvection(*this);
      }

      private:
        scalar coeff1, coeff2;
        GeomType gt;
    };

    /* Default volumetric vector form \int_{area} spline_coeff1(u_ext[0]) * u->dx * v->val
       + spline_coeff2(u_ext[0]) * u->dy * v->val d\bfx
       spline_coeff1, spline_coeff2... non-constant parameters given by cubic splines
    */

    class DefaultResidualNonlinearAdvection : public WeakForm::VectorFormVol
    {
    public:
    DefaultResidualNonlinearAdvection(int i, CubicSpline* spline_coeff1,
                                      CubicSpline* spline_coeff2, scalar const_coeff = 1.0, GeomType gt = HERMES_PLANAR)
      : WeakForm::VectorFormVol(i), spline_coeff1(spline_coeff1), spline_coeff2(spline_coeff2), const_coeff(const_coeff), gt(gt) { }
      DefaultResidualNonlinearAdvection(int i, std::string area, CubicSpline* spline_coeff1,
                                        CubicSpline* spline_coeff2, scalar const_coeff = 1.0, GeomType gt = HERMES_PLANAR)
  : WeakForm::VectorFormVol(i, area), spline_coeff1(spline_coeff1), spline_coeff2(spline_coeff2), const_coeff(const_coeff),
  gt(gt) { }

      template<typename Real, typename Scalar>
      Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                         Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
        Scalar result = 0;
        Func<Scalar>* u_prev = u_ext[0];
        for (int i = 0; i < n; i++) {
          result += wt[i] * (const_coeff*spline_coeff1->get_value(u_prev->val[i]) * (u_prev->dx[i] * v->val[i])
                             + const_coeff*spline_coeff2->get_value(u_prev->val[i]) * (u_prev->dy[i] * v->val[i]));
        }
        return result;
      }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                           Geom<double> *e, ExtData<scalar> *ext) const {
        return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
              Geom<Ord> *e, ExtData<Ord> *ext) const {
        return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::VectorFormVol* clone() {
        return new DefaultResidualNonlinearAdvection(*this);
      }

      private:
        CubicSpline* spline_coeff1, *spline_coeff2;
        scalar const_coeff;
        GeomType gt;
    };

    /* Default volumetric vector form \int_{area} rhs(x, y) v d\bfx
       rhs(x, y)... non-constant right-hand side
    */
    class DefaultVectorFormNonConst : public WeakForm::VectorFormVol
    {
    public:
      DefaultVectorFormNonConst(int i, RightHandSides::DefaultNonConstRightHandSide* rhs,
                                GeomType gt = HERMES_PLANAR)
  : WeakForm::VectorFormVol(i), rhs(rhs), gt(gt) { }
      DefaultVectorFormNonConst(int i, std::string area, RightHandSides::DefaultNonConstRightHandSide* rhs,
                                GeomType gt = HERMES_PLANAR)
  : WeakForm::VectorFormVol(i, area), rhs(rhs), gt(gt) { }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                   Geom<double> *e, ExtData<scalar> *ext) const {
        scalar result = 0;
        for (int i = 0; i < n; i++)
          result += wt[i] * (rhs->value(e->x[i], e->y[i]) * v->val[i]);
        return result;
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
              Geom<Ord> *e, ExtData<Ord> *ext) const {
        Ord result = 0;
        for (int i = 0; i < n; i++)
          result += wt[i] * (rhs->ord(e->x[i], e->y[i]) * v->val[i]);
        return result;
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::VectorFormVol* clone() {
        return new DefaultVectorFormNonConst(*this);
      }

      private:
        RightHandSides::DefaultNonConstRightHandSide* rhs;
        GeomType gt;
    };
  }

  namespace SurfaceMatrixForms {

    /* Default surface matrix form \int_{area} coeff u v dS
       coeff... constant number
    */

    class DefaultMatrixFormSurf : public WeakForm::MatrixFormSurf
    {
    public:
      DefaultMatrixFormSurf(int i, int j, scalar coeff, GeomType gt = HERMES_PLANAR)
  : WeakForm::MatrixFormSurf(i, j), coeff(coeff), gt(gt) { }
      DefaultMatrixFormSurf(int i, int j, std::string area, scalar coeff, GeomType gt = HERMES_PLANAR)
  : WeakForm::MatrixFormSurf(i, j, area), coeff(coeff), gt(gt) { }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                           Geom<double> *e, ExtData<scalar> *ext) const {
        scalar result = 0;
        if (gt == HERMES_PLANAR) result = int_u_v<double, scalar>(n, wt, u, v);
        else if (gt == HERMES_AXISYM_X) result = int_y_u_v<double, scalar>(n, wt, u, v, e);
        else result = int_x_u_v<double, scalar>(n, wt, u, v, e);
        return coeff * result;
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                      Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
        Ord result = 0;
        if (gt == HERMES_PLANAR) result = int_u_v<Ord, Ord>(n, wt, u, v);
        else if (gt == HERMES_AXISYM_X) result = int_y_u_v<Ord, Ord>(n, wt, u, v, e);
        else result = int_x_u_v<Ord, Ord>(n, wt, u, v, e);
        return result;
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::MatrixFormSurf* clone() {
        return new DefaultMatrixFormSurf(*this);
      }

      private:
        scalar coeff;
        GeomType gt;
    };

    /* Default surface matrix form \int_{area} spline_coeff'(u_ext[0]) u_ext[0] u v + spline_coeff(u_ext[0]) u v dS
       spline_coeff... non-constant parameter given by a spline
    */

    class DefaultJacobianFormSurf : public WeakForm::MatrixFormSurf
    {
    public:
      DefaultJacobianFormSurf(int i, int j, CubicSpline* spline_coeff, scalar const_coeff = 1.0,
                              GeomType gt = HERMES_PLANAR)
  : WeakForm::MatrixFormSurf(i, j), spline_coeff(spline_coeff), const_coeff(const_coeff), gt(gt) { }
      DefaultJacobianFormSurf(int i, int j, std::string area, CubicSpline* spline_coeff,
                              scalar const_coeff = 1.0, GeomType gt = HERMES_PLANAR)
  : WeakForm::MatrixFormSurf(i, j, area), spline_coeff(spline_coeff), const_coeff(const_coeff), gt(gt) { }

      template<typename Real, typename Scalar>
      Scalar matrix_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
        Scalar result = 0;
        for (int i = 0; i < n; i++) {
          result += wt[i] * (const_coeff*spline_coeff->get_derivative(u_ext[0]->val[i]) * u_ext[0]->val[i]
                             + const_coeff*spline_coeff->get_value(u_ext[0]->val[i]))
                    * u->val[i] * v->val[i];
        }
        return result;
      }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                   Geom<double> *e, ExtData<scalar> *ext) const {
        return matrix_form_surf<double, scalar>(n, wt, u_ext, u, v, e, ext);
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
              Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
        return matrix_form_surf<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::MatrixFormSurf* clone() {
        return new DefaultJacobianFormSurf(*this);
      }

      private:
        CubicSpline* spline_coeff;
        scalar const_coeff;
        GeomType gt;
    };
  }

  namespace SurfaceVectorForms {

    /* Default surface vector form \int_{area} coeff v dS
       coeff... constant number
    */

    class DefaultVectorFormSurf : public WeakForm::VectorFormSurf
    {
    public:
      DefaultVectorFormSurf(int i, scalar coeff, GeomType gt = HERMES_PLANAR)
      : WeakForm::VectorFormSurf(i), coeff(coeff), gt(gt) { }
      DefaultVectorFormSurf(int i, std::string area, scalar coeff, GeomType gt = HERMES_PLANAR)
      : WeakForm::VectorFormSurf(i, area), coeff(coeff), gt(gt) { }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                           Geom<double> *e, ExtData<scalar> *ext) const {
        scalar result = 0;
        if (gt == HERMES_PLANAR) result = int_v<double>(n, wt, v);
        else if (gt == HERMES_AXISYM_X) result = int_y_v<double>(n, wt, v, e);
        else result = int_x_v<double>(n, wt, v, e);

        return coeff * result;
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e,
                      ExtData<Ord> *ext) const {
        Ord result = 0;
        if (gt == HERMES_PLANAR) result = int_v<Ord>(n, wt, v);
        else if (gt == HERMES_AXISYM_X) result = int_y_v<Ord>(n, wt, v, e);
        else result = int_x_v<Ord>(n, wt, v, e);

        return result;
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::VectorFormSurf* clone() {
        return new DefaultVectorFormSurf(*this);
      }

      private:
        scalar coeff;
        GeomType gt;
    };

    class MultiComponentDefaultVectorFormSurf : public WeakForm::MultiComponentVectorFormSurf
    {
    public:
      MultiComponentDefaultVectorFormSurf(Hermes::vector<unsigned int> coordinates,
                                          Hermes::vector<scalar> coeffs, GeomType gt = HERMES_PLANAR)
      : WeakForm::MultiComponentVectorFormSurf(coordinates), coeffs(coeffs), gt(gt) { }
      MultiComponentDefaultVectorFormSurf(Hermes::vector<unsigned int> coordinates, std::string area,
                                          Hermes::vector<scalar> coeffs, GeomType gt = HERMES_PLANAR)
      : WeakForm::MultiComponentVectorFormSurf(coordinates, area), coeffs(coeffs), gt(gt) { }

      virtual void value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                          Geom<double> *e, ExtData<scalar> *ext, Hermes::vector<scalar>& result) const {
        scalar result_base = 0;
        if (gt == HERMES_PLANAR)
          result_base = int_v<double>(n, wt, v);
        else
          if (gt == HERMES_AXISYM_X)
            result_base = int_y_v<double>(n, wt, v, e);
          else
            result_base = int_x_v<double>(n, wt, v, e);

        for(unsigned int result_i = 0; result_i < this->coordinates.size(); result_i++)
          result.push_back(result_base * coeffs[result_i]);
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                          Geom<Ord> *e, ExtData<Ord> *ext) const {
        if (gt == HERMES_PLANAR)
          return int_v<Ord>(n, wt, v);
        else
          if (gt == HERMES_AXISYM_X)
            return int_y_v<Ord>(n, wt, v, e);
          else
            return int_x_v<Ord>(n, wt, v, e);
      }

      /*
      // This is to make the form usable in rk_time_step().
      virtual WeakForm::VectorFormSurf* clone() {
        return new MultiComponentDefaultVectorFormSurf(*this);
      }
      */

      private:
        Hermes::vector<scalar> coeffs;
        GeomType gt;
    };

    /* Default surface vector form \int_{area} coeff * u_ext[0] v dS
       coeff... constant parameter
    */

    class DefaultResidualSurfConst : public WeakForm::VectorFormSurf
    {
    public:
      DefaultResidualSurfConst(int i, scalar coeff = 1.0,
                               GeomType gt = HERMES_PLANAR)
  : WeakForm::VectorFormSurf(i), coeff(coeff), gt(gt) { }
      DefaultResidualSurfConst(int i, std::string area, scalar coeff = 1.0,
                               GeomType gt = HERMES_PLANAR)
  : WeakForm::VectorFormSurf(i, area), coeff(coeff), gt(gt) { }

      template<typename Real, typename Scalar>
      Scalar vector_form_surf(int n, double *wt, Func<Scalar> *u_ext[],
                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
        Scalar result = 0;
        for (int i = 0; i < n; i++) {
          result += wt[i] * u_ext[0]->val[i] * v->val[i];
        }
        return coeff * result;
      }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                   Geom<double> *e, ExtData<scalar> *ext) const {
        return vector_form_surf<double, scalar>(n, wt, u_ext, v, e, ext);
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                      Geom<Ord> *e, ExtData<Ord> *ext) const {
        return vector_form_surf<Ord, Ord>(n, wt, u_ext, v, e, ext);
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::VectorFormSurf* clone() {
        return new DefaultResidualSurfConst(*this);
      }

      private:
        scalar coeff;
        GeomType gt;
    };

    /* Default surface vector form \int_{area} spline_coeff(u_ext[0]) v dS
       spline_coeff... non-constant parameter given by cubic spline
    */

    class DefaultResidualSurfSpline : public WeakForm::VectorFormSurf
    {
    public:
      DefaultResidualSurfSpline(int i, CubicSpline* spline_coeff, scalar const_coeff = 1.0,
                              GeomType gt = HERMES_PLANAR)
  : WeakForm::VectorFormSurf(i), spline_coeff(spline_coeff), const_coeff(const_coeff), gt(gt) { }
      DefaultResidualSurfSpline(int i, std::string area, CubicSpline* spline_coeff, scalar const_coeff = 1.0,
                              GeomType gt = HERMES_PLANAR)
  : WeakForm::VectorFormSurf(i, area), spline_coeff(spline_coeff), const_coeff(const_coeff), gt(gt) { }

      template<typename Real, typename Scalar>
      Scalar vector_form_surf(int n, double *wt, Func<Scalar> *u_ext[],
                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
        Scalar result = 0;
        for (int i = 0; i < n; i++) {
          result += wt[i] * const_coeff*spline_coeff->get_value(u_ext[0]->val[i]) * v->val[i];
        }
        return result;
      }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                   Geom<double> *e, ExtData<scalar> *ext) const {
        return vector_form_surf<double, scalar>(n, wt, u_ext, v, e, ext);
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                      Geom<Ord> *e, ExtData<Ord> *ext) const {
        return vector_form_surf<Ord, Ord>(n, wt, u_ext, v, e, ext);
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::VectorFormSurf* clone() {
        return new DefaultResidualSurfSpline(*this);
      }

      private:
        CubicSpline* spline_coeff;
        scalar const_coeff;
        GeomType gt;
    };
  }

  namespace WeakForms {
    /* Default weak form for the Laplace equation -Laplace u = 0
    */

    class DefaultWeakFormLaplace : public WeakForm
    {
    public:
      DefaultWeakFormLaplace() : WeakForm(1)
      {
        add_matrix_form(new VolumetricMatrixForms::DefaultLinearDiffusion(0, 0));
      };
    };
  }
}

#endif
