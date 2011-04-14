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

#include "../integrals/integrals_h1.h"

namespace WeakFormsMaxwell {

  namespace VolumetricMatrixForms {

    /* Default volumetric matrix form \int_{area} coeff \curl u \cdot \curl v d\bfx
       coeff... constant number
    */

    class DefaultLinearMagnetostatics : public WeakForm::MatrixFormVol
    {
    public:
      // The optional order_increase takes into account the axisymmetric part.
      DefaultLinearMagnetostatics(int i, int j, scalar coeff = 1.0,
                                  SymFlag sym = HERMES_SYM, GeomType gt = HERMES_PLANAR, 
                                  int order_increase = 3)
             : WeakForm::MatrixFormVol(i, j, sym), coeff(coeff), gt(gt), order_increase(order_increase) { }
      DefaultLinearMagnetostatics(int i, int j, std::string area, scalar coeff = 1.0,
                                  SymFlag sym = HERMES_SYM, GeomType gt = HERMES_PLANAR, int order_increase = 3)
             : WeakForm::MatrixFormVol(i, j, sym, area), coeff(coeff), gt(gt), order_increase(order_increase) { }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                           Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
        scalar planar_part = int_grad_u_grad_v<double, scalar>(n, wt, u, v);
        scalar axisym_part = 0;
        if (gt == HERMES_AXISYM_X)
          axisym_part = int_u_dvdy_over_y<double, scalar>(n, wt, u, v, e);
        else if (gt == HERMES_AXISYM_Y)
          axisym_part = int_u_dvdx_over_x<double, scalar>(n, wt, u, v, e);

        return coeff * (planar_part + axisym_part);
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
              Geom<Ord> *e, ExtData<Ord> *ext) const {
        Ord planar_part = int_grad_u_grad_v<Ord, Ord>(n, wt, u, v);

        // This increase is for the axisymmetric part. We are not letting the
        // Ord class do it since it would automatically choose the highest order
        // due to the nonpolynomial 1/r term.
        return planar_part * Ord(order_increase);
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::MatrixFormVol* clone() {
        return new DefaultLinearMagnetostatics(*this);
      }

      private:
        scalar coeff;
        GeomType gt;
        int order_increase;
    };

    /* Default volumetric matrix form \int_{area} coeff_spline(u_ext[0]) \curl u \curl v d\bfx
       spline_coeff... nonconstant parameter given by cubic spline
    */

    class DefaultJacobianNonlinearMagnetostatics : public WeakForm::MatrixFormVol
    {
    public:
      DefaultJacobianNonlinearMagnetostatics(int i, int j, CubicSpline* spline_coeff, 
                                             scalar const_coeff = 1.0,
                                             SymFlag sym = HERMES_NONSYM, 
                                             GeomType gt = HERMES_PLANAR,
                                             int order_increase = 3)
             : WeakForm::MatrixFormVol(i, j, sym), spline_coeff(spline_coeff), 
                                       const_coeff(const_coeff), gt(gt),
                                       order_increase(order_increase) { }
      DefaultJacobianNonlinearMagnetostatics(int i, int j, std::string area,
                                             CubicSpline* spline_coeff, 
                                             scalar const_coeff = 1.0, 
                                             SymFlag sym = HERMES_NONSYM,
                                             GeomType gt = HERMES_PLANAR, 
                                             int order_increase = 3)
             : WeakForm::MatrixFormVol(i, j, sym, area), spline_coeff(spline_coeff), 
                                       const_coeff(const_coeff), gt(gt),
                                       order_increase(order_increase) { }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                           Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
        scalar planar_part = 0;
        scalar axisym_part = 0;
        for (int i = 0; i < n; i++) {
          scalar B_i = sqrt(sqr(u_ext[0]->dx[i]) + sqr(u_ext[0]->dy[i]));
          if (std::abs(B_i) > 1e-12) {
            planar_part += wt[i] * const_coeff*spline_coeff->get_derivative(B_i) / B_i
                                 * (u_ext[0]->dx[i] * u->dx[i] + u_ext[0]->dy[i] * u->dy[i])
                                 * (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i]);
            if (gt == HERMES_AXISYM_X) {
              axisym_part += wt[i] * const_coeff*spline_coeff->get_derivative(B_i) / B_i / e->y[i]
       		                   * (u_ext[0]->dx[i] * u->dx[i] + u_ext[0]->dy[i] * u->dy[i])
                                   * u_ext[0]->val[i] * v->dy[i];
            }
            else if (gt == HERMES_AXISYM_Y) {
              axisym_part += wt[i] * const_coeff*spline_coeff->get_derivative(B_i) / B_i / e->x[i]
 		                   * (u_ext[0]->dx[i] * u->dx[i] + u_ext[0]->dy[i] * u->dy[i])
                                   * u_ext[0]->val[i] * v->dx[i];
            }
          }
          planar_part += wt[i] * const_coeff*spline_coeff->get_value(B_i)
                               * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
          if (gt == HERMES_AXISYM_X) {
            axisym_part += wt[i] * const_coeff*spline_coeff->get_value(B_i) / e->y[i]
                                 * u->val[i] * v->dy[i];
          }
          else if (gt == HERMES_AXISYM_Y) {
            axisym_part += wt[i] * const_coeff*spline_coeff->get_value(B_i) / e->x[i]
                                 * u->val[i] * v->dx[i];
          }
        }

        return planar_part + axisym_part;
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
              Geom<Ord> *e, ExtData<Ord> *ext) const {
        Ord planar_part = 0;
        for (int i = 0; i < n; i++) {
          Ord B_i = sqrt(sqr(u_ext[0]->dx[i]) + sqr(u_ext[0]->dy[i]));
          planar_part += wt[i] * const_coeff*spline_coeff->get_derivative(B_i) / B_i
                               * (u_ext[0]->dx[i] * u->dx[i] + u_ext[0]->dy[i] * u->dy[i])
                               * (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i]);
          planar_part += wt[i] * const_coeff*spline_coeff->get_value(B_i)
                               * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
        }

        // This increase is for the axisymmetric part. We are not letting the
        // Ord class do it since it would automatically choose the highest order
        // due to the nonpolynomial 1/r term.
        return planar_part * Ord(order_increase);
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::MatrixFormVol* clone() {
        return new DefaultJacobianNonlinearMagnetostatics(*this);
      }

      private:
        CubicSpline* spline_coeff;
        scalar const_coeff;
        GeomType gt;
        int order_increase;
    };
  }

  namespace VolumetricVectorForms {

    /* Default volumetric vector form \int_{area} coeff
       \nabla u_ext[0] \cdot \nabla v d\bfx
       coeff... constant parameter
    */

    class DefaultResidualLinearMagnetostatics : public WeakForm::VectorFormVol
    {
    public:
      DefaultResidualLinearMagnetostatics(int i, scalar coeff, GeomType gt = HERMES_PLANAR,
                                          int order_increase = 3)
             : WeakForm::VectorFormVol(i), coeff(coeff), gt(gt), order_increase(order_increase) { }
      DefaultResidualLinearMagnetostatics(int i, std::string area, scalar coeff,
                                          GeomType gt = HERMES_PLANAR, int order_increase = 3)
             : WeakForm::VectorFormVol(i, area), coeff(coeff), gt(gt), order_increase(order_increase) { }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                           Geom<double> *e, ExtData<scalar> *ext) const {
        scalar planar_part = int_grad_u_grad_v<double, scalar>(n, wt, u_ext[0], v);
        scalar axisym_part = 0;
        if (gt == HERMES_AXISYM_X)
          axisym_part = int_u_dvdy_over_y<double, scalar>(n, wt, u_ext[0], v, e);
        else if (gt == HERMES_AXISYM_Y)
          axisym_part = int_u_dvdx_over_x<double, scalar>(n, wt, u_ext[0], v, e);

        return coeff * (planar_part + axisym_part);
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
              Geom<Ord> *e, ExtData<Ord> *ext) const {
        Ord planar_part = int_grad_u_grad_v<Ord, Ord>(n, wt, u_ext[0], v);
        return planar_part * Ord(order_increase);
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::VectorFormVol* clone() {
        return new DefaultResidualLinearMagnetostatics(*this);
      }

      private:
        scalar coeff;
        GeomType gt;
        int order_increase;
    };

    /* Default volumetric vector form \int_{area} spline_coeff(u_ext[0])
       \nabla u_ext[0] \cdot \nabla v d\bfx
       spline_coeff... non-constant parameter given by a cubic spline
    */

    class DefaultResidualNonlinearMagnetostatics : public WeakForm::VectorFormVol
    {
    public:
      DefaultResidualNonlinearMagnetostatics(int i, CubicSpline* spline_coeff, 
                                             scalar const_coeff = 1.0, 
                                             GeomType gt = HERMES_PLANAR,
                                             int order_increase = 3)
             : WeakForm::VectorFormVol(i), spline_coeff(spline_coeff), const_coeff(const_coeff), 
                                           gt(gt), order_increase(order_increase) { }
      DefaultResidualNonlinearMagnetostatics(int i, std::string area, CubicSpline* spline_coeff, 
                                             scalar const_coeff = 1.0,
                                             GeomType gt = HERMES_PLANAR, int order_increase = 3)
             : WeakForm::VectorFormVol(i, area), spline_coeff(spline_coeff), const_coeff(const_coeff), gt(gt), 
                                       order_increase(order_increase) { }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                           Geom<double> *e, ExtData<scalar> *ext) const {
        scalar planar_part = 0;
        scalar axisym_part = 0;
        for (int i = 0; i < n; i++) {
          scalar B_i = sqrt(sqr(u_ext[0]->dx[i]) + sqr(u_ext[0]->dy[i]));
          planar_part += wt[i] * const_coeff*spline_coeff->get_value(B_i) *
                                 (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i]);
          if (gt == HERMES_AXISYM_X) axisym_part += wt[i] * const_coeff*spline_coeff->get_value(B_i) / e->y[i]
                                     * u_ext[0]->val[i] * v->dy[i];
          else if (gt == HERMES_AXISYM_Y) axisym_part += wt[i] * const_coeff*spline_coeff->get_value(B_i) / e->x[i]
                                     * u_ext[0]->val[i] * v->dx[i];
        }
        return planar_part + axisym_part;
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                      Geom<Ord> *e, ExtData<Ord> *ext) const {
        Ord planar_part = 0;
        for (int i = 0; i < n; i++) {
          Ord B_i = sqrt(sqr(u_ext[0]->dx[i]) + sqr(u_ext[0]->dy[i]));
          planar_part += wt[i] * const_coeff*spline_coeff->get_value(B_i) *
                                 (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i]);
        }
        return planar_part * Ord(order_increase);

      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::VectorFormVol* clone() {
        return new DefaultResidualNonlinearMagnetostatics(*this);
      }

      private:
        CubicSpline* spline_coeff;
        scalar const_coeff;
        GeomType gt;
        int order_increase;
    };
  }
}

#endif
