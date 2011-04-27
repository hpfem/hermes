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

  // Generic class for functions of x, y in weak forms.
  class DefaultFunction
  {
  public:
    DefaultFunction()
    {
      this->is_constant = false;
      this->const_value = -9999;
    }
    DefaultFunction(scalar value)
    {
      this->is_constant = true;
      this->const_value = value;
    };

    virtual scalar value(double x, double y) const
    {
      return const_value;
    };
    virtual Ord ord(Ord x, Ord y) const
    {
    return Ord(0);
    };

  protected:
    bool is_constant;
    scalar const_value;
  };

namespace WeakFormsH1 {

  namespace VolumetricMatrixForms {

    /* Default volumetric matrix form \int_{area} const_coeff * function_coeff(x, y) * u * v \bfx
       const_coeff... constant number
       function_coeff... (generally nonconstant) function of x, y
    */

    class DefaultMatrixFormVol : public WeakForm::MatrixFormVol
    {
    public:
      DefaultMatrixFormVol(int i, int j, std::string area = HERMES_ANY,
                           scalar const_coeff = 1.0, DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                           SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR)
        : WeakForm::MatrixFormVol(i, j, area, sym), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
      {
        // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
        if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
      }
      DefaultMatrixFormVol(int i, int j, Hermes::vector<std::string> areas,
                            scalar const_coeff, DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                            SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR)
        : WeakForm::MatrixFormVol(i, j, areas, sym), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
      {
        // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
        if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
      }

      ~DefaultMatrixFormVol() {
        if (function_coeff != HERMES_DEFAULT_FUNCTION) delete function_coeff;
      };

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                           Geom<double> *e, ExtData<scalar> *ext) const {
        scalar result = 0;
        if (gt == HERMES_PLANAR) {
    for (int i = 0; i < n; i++) {
      result += wt[i] * function_coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
    }
        }
        else {
          if (gt == HERMES_AXISYM_X) {
      for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * function_coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
      }
          }
          else {
      for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * function_coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
      }
          }
        }

  return const_coeff * result;
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                      Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
        Ord result = 0;
        if (gt == HERMES_PLANAR) {
    for (int i = 0; i < n; i++) {
      result += wt[i] * function_coeff->ord(e->x[i], e->y[i]) * u->val[i] * v->val[i];
    }
        }
        else {
          if (gt == HERMES_AXISYM_X) {
      for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * function_coeff->ord(e->x[i], e->y[i]) * u->val[i] * v->val[i];
      }
          }
          else {
      for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * function_coeff->ord(e->x[i], e->y[i]) * u->val[i] * v->val[i];
      }
          }
        }

  return result;
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::MatrixFormVol* clone() {
        return new DefaultMatrixFormVol(*this);
      }

      private:
        scalar const_coeff;
        DefaultFunction* function_coeff;
        GeomType gt;
    };

    /* Default volumetric matrix form \int_{area} const_coeff * spline_coeff'(u_ext[0]) u \nabla u_ext[0] \cdot \nabla v
       + const_coeff * spline_coeff(u_ext[0]) * \nabla u \cdot \nabla v d\bfx
       const_coeff... constant number
       spline_coeff... nonconstant parameter given by cubic spline
    */

    class DefaultJacobianDiffusion : public WeakForm::MatrixFormVol
    {
    public:
      DefaultJacobianDiffusion(int i, int j, std::string area = HERMES_ANY, scalar const_coeff = 1.0,
                               CubicSpline* c_spline = HERMES_DEFAULT_SPLINE,
                               SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR)
        : WeakForm::MatrixFormVol(i, j, area, sym), const_coeff(const_coeff), spline_coeff(c_spline), gt(gt)
      {
        // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
        if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
      };

    DefaultJacobianDiffusion(int i, int j, Hermes::vector<std::string> areas, scalar const_coeff = 1.0,
                             CubicSpline* c_spline = HERMES_DEFAULT_SPLINE,
                             SymFlag sym = HERMES_NONSYM, GeomType gt = HERMES_PLANAR)
        : WeakForm::MatrixFormVol(i, j, areas, sym), const_coeff(const_coeff), spline_coeff(c_spline), gt(gt)
      {
        // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
        if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
      }

      ~DefaultJacobianDiffusion() {
        if (spline_coeff != HERMES_DEFAULT_SPLINE) delete spline_coeff;
      };

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                   Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const
      {
        scalar result = 0;
        if (gt == HERMES_PLANAR) {
          for (int i = 0; i < n; i++) {
            result += wt[i] * (const_coeff*spline_coeff->get_derivative(u_ext[0]->val[i]) * u->val[i] *
                   (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i])
                   + const_coeff*spline_coeff->get_value(u_ext[0]->val[i])
                                     * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
          }
        }
        else {
          if (gt == HERMES_AXISYM_X) {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->y[i] * (const_coeff*spline_coeff->get_derivative(u_ext[0]->val[i]) * u->val[i] *
                        (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i])
                        + const_coeff*spline_coeff->get_value(u_ext[0]->val[i])
         * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
      }
          }
          else {
      for (int i = 0; i < n; i++) {
              result += wt[i] * e->x[i] * (const_coeff*spline_coeff->get_derivative(u_ext[0]->val[i]) * u->val[i] *
                        (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i])
                        + const_coeff*spline_coeff->get_value(u_ext[0]->val[i])
               * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
      }
          }
        }

        return result;
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
              Geom<Ord> *e, ExtData<Ord> *ext) const
      {
        Ord result = 0;
        if (gt == HERMES_PLANAR) {
          for (int i = 0; i < n; i++) {
            result += wt[i] * (const_coeff*spline_coeff->get_derivative(u_ext[0]->val[i]) * u->val[i] *
                      (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i])
                      + const_coeff*spline_coeff->get_value(u_ext[0]->val[i])
                                     * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
          }
        }
        else {
          if (gt == HERMES_AXISYM_X) {
      for (int i = 0; i < n; i++) {
              result += wt[i] * e->y[i] * (const_coeff*spline_coeff->get_derivative(u_ext[0]->val[i]) * u->val[i] *
                        (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i])
                        + const_coeff*spline_coeff->get_value(u_ext[0]->val[i])
             * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
      }
          }
          else {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->x[i] * (const_coeff*spline_coeff->get_derivative(u_ext[0]->val[i]) * u->val[i] *
                        (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i])
                        + const_coeff*spline_coeff->get_value(u_ext[0]->val[i])
                         * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
            }
          }
        }

        return result;
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::MatrixFormVol* clone() {
        return new DefaultJacobianDiffusion(*this);
      }

      private:
        scalar const_coeff;
        CubicSpline* spline_coeff;
        GeomType gt;
    };

    /* Default volumetric matrix form
       \int_{area} spline_coeff1`(u_ext[0]) * u * u_ext[0]->dx * v
       + spline_coeff1(u_ext[0]) * u->dx * v
       + spline_coeff2`(u_ext[0]) * u * u_ext[0]->dy * v
       + spline_coeff2(u_ext[0]) * u->dy * v d\bfx.
       spline_coeff1, spline_coeff2... non-constant parameters given by cubic splines
    */

    class DefaultJacobianAdvection : public WeakForm::MatrixFormVol
    {
    public:
      DefaultJacobianAdvection(int i, int j, std::string area = HERMES_ANY, scalar const_coeff = 1.0,
                               CubicSpline* c_spline1 = HERMES_DEFAULT_SPLINE,
                               CubicSpline* c_spline2 = HERMES_DEFAULT_SPLINE, GeomType gt = HERMES_PLANAR)
        : WeakForm::MatrixFormVol(i, j, area, HERMES_NONSYM), const_coeff(const_coeff), spline_coeff1(c_spline1),
                                  spline_coeff2(c_spline2), gt(gt)
      {
        if (gt != HERMES_PLANAR) error("Axisymmetric advection forms not implemented yet.");

        // If spline1 == HERMES_DEFAULT_SPLINE or spline2 == HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
        if (c_spline1 == HERMES_DEFAULT_SPLINE) this->spline_coeff1 = new CubicSpline(1.0);
        if (c_spline2 == HERMES_DEFAULT_SPLINE) this->spline_coeff2 = new CubicSpline(1.0);
      }

     DefaultJacobianAdvection(int i, int j, Hermes::vector<std::string> areas, scalar const_coeff = 1.0,
                              CubicSpline* c_spline1 = HERMES_DEFAULT_SPLINE,
                              CubicSpline* c_spline2 = HERMES_DEFAULT_SPLINE,
                              GeomType gt = HERMES_PLANAR)
       : WeakForm::MatrixFormVol(i, j, areas, HERMES_NONSYM), const_coeff(const_coeff), spline_coeff1(spline_coeff1),
                           spline_coeff2(spline_coeff2), gt(gt)
     {
       if (gt != HERMES_PLANAR) error("Axisymmetric advection forms not implemented yet.");

       // If spline1 == HERMES_DEFAULT_SPLINE or spline2 == HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
       if (c_spline1 == HERMES_DEFAULT_SPLINE) this->spline_coeff1 = new CubicSpline(1.0);
       if (c_spline2 == HERMES_DEFAULT_SPLINE) this->spline_coeff2 = new CubicSpline(1.0);
     }

     ~DefaultJacobianAdvection() {
       if (spline_coeff1 != HERMES_DEFAULT_SPLINE) delete spline_coeff1;
       if (spline_coeff2 != HERMES_DEFAULT_SPLINE) delete spline_coeff2;
     };

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
        return new DefaultJacobianAdvection(*this);
      }

      private:
      scalar const_coeff;
      CubicSpline* spline_coeff1, *spline_coeff2;
      GeomType gt;
    };
  }

  namespace VolumetricVectorForms {

    /* Default volumetric vector form \int_{area} const_coeff * function_coeff(x, y) * v d\bfx
       const_coeff... constant number
       function_coeff... (generally nonconstant) function of x, y
    */

    class DefaultVectorFormVol : public WeakForm::VectorFormVol
    {
    public:
      DefaultVectorFormVol(int i, std::string area = HERMES_ANY, scalar const_coeff = 1.0,
                           DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                           GeomType gt = HERMES_PLANAR)
             : WeakForm::VectorFormVol(i, area), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
      {
        // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
        if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
      }
      DefaultVectorFormVol(int i, Hermes::vector<std::string> areas, scalar const_coeff = 1.0,
                           DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                           GeomType gt = HERMES_PLANAR)
             : WeakForm::VectorFormVol(i, areas), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
      {
        // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
        if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
      }

      ~DefaultVectorFormVol() {
        if (function_coeff != HERMES_DEFAULT_FUNCTION) delete function_coeff;
      };

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                           Geom<double> *e, ExtData<scalar> *ext) const {
        scalar result = 0;
        if (gt == HERMES_PLANAR) {
          for (int i = 0; i < n; i++) {
            result += wt[i] * function_coeff->value(e->x[i], e->y[i]) * v->val[i];
          }
        }
        else {
          if (gt == HERMES_AXISYM_X) {
            for (int i = 0; i < n; i++) {
                result += wt[i] * e->y[i] * function_coeff->value(e->x[i], e->y[i]) * v->val[i];
            }
          }
          else {
            for (int i = 0; i < n; i++) {
                result += wt[i] * e->x[i] * function_coeff->value(e->x[i], e->y[i]) * v->val[i];
            }
          }
        }
        return const_coeff * result;
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                      Geom<Ord> *e, ExtData<Ord> *ext) const {
        Ord result = 0;
        if (gt == HERMES_PLANAR) {
    for (int i = 0; i < n; i++) {
      result += wt[i] * function_coeff->ord(e->x[i], e->y[i]) * v->val[i];
    }
        }
        else {
          if (gt == HERMES_AXISYM_X) {
      for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * function_coeff->ord(e->x[i], e->y[i]) * v->val[i];
      }
          }
          else {
      for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * function_coeff->ord(e->x[i], e->y[i]) * v->val[i];
      }
          }
        }

  return result;
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::VectorFormVol* clone() {
        return new DefaultVectorFormVol(*this);
      }

      private:
        scalar const_coeff;
        DefaultFunction* function_coeff;
        GeomType gt;
    };

    /* Default volumetric vector form \int_{area} const_coeff * function_coeff(x, y) * u_ext[0] * v d\bfx
       const_coeff... constant number
       function_coeff... (generally nonconstant) function of x, y
    */

    class DefaultResidualVol : public WeakForm::VectorFormVol
    {
    public:
      DefaultResidualVol(int i, std::string area = HERMES_ANY, scalar const_coeff = 1.0,
                           DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                           GeomType gt = HERMES_PLANAR)
             : WeakForm::VectorFormVol(i, area), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
      {
        // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
        if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
      }
      DefaultResidualVol(int i, Hermes::vector<std::string> areas, scalar const_coeff = 1.0,
                           DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                           GeomType gt = HERMES_PLANAR)
             : WeakForm::VectorFormVol(i, areas), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
      {
        // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
        if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
      }

      ~DefaultResidualVol() {
        if (function_coeff != HERMES_DEFAULT_FUNCTION) delete function_coeff;
      };

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                           Geom<double> *e, ExtData<scalar> *ext) const {
        scalar result = 0;
        if (gt == HERMES_PLANAR) {
          for (int i = 0; i < n; i++) {
            result += wt[i] * function_coeff->value(e->x[i], e->y[i]) * u_ext[0]->val[i] * v->val[i];
          }
        }
        else {
          if (gt == HERMES_AXISYM_X) {
            for (int i = 0; i < n; i++) {
                result += wt[i] * e->y[i] * function_coeff->value(e->x[i], e->y[i]) * u_ext[0]->val[i] * v->val[i];
            }
          }
          else {
            for (int i = 0; i < n; i++) {
                result += wt[i] * e->x[i] * function_coeff->value(e->x[i], e->y[i]) * u_ext[0]->val[i] * v->val[i];
            }
          }
        }
        return const_coeff * result;
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                      Geom<Ord> *e, ExtData<Ord> *ext) const {
        Ord result = 0;
        if (gt == HERMES_PLANAR) {
    for (int i = 0; i < n; i++) {
      result += wt[i] * function_coeff->ord(e->x[i], e->y[i]) * u_ext[0]->val[i] * v->val[i];
    }
        }
        else {
          if (gt == HERMES_AXISYM_X) {
      for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * function_coeff->ord(e->x[i], e->y[i]) * u_ext[0]->val[i] * v->val[i];
      }
          }
          else {
      for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * function_coeff->ord(e->x[i], e->y[i]) * u_ext[0]->val[i] * v->val[i];
      }
          }
        }

  return result;
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::VectorFormVol* clone() {
        return new DefaultResidualVol(*this);
      }

      private:
        scalar const_coeff;
        DefaultFunction* function_coeff;
        GeomType gt;
    };

    /* Default volumetric vector form \int_{area} const_coeff * spline_coeff(u_ext[0]) *
       \nabla u_ext[0] \cdot \nabla v d\bfx
       const_coeff... constant number
       spline_coeff... non-constant parameter given by a cubic spline
    */

    class DefaultResidualDiffusion : public WeakForm::VectorFormVol
    {
    public:
      DefaultResidualDiffusion(int i, std::string area = HERMES_ANY, scalar const_coeff = 1.0,
                               CubicSpline* c_spline = HERMES_DEFAULT_SPLINE,
                               GeomType gt = HERMES_PLANAR)
        : WeakForm::VectorFormVol(i, area), const_coeff(const_coeff), spline_coeff(c_spline), gt(gt)
      {
        // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
        if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
      };

      DefaultResidualDiffusion(int i, Hermes::vector<std::string> areas, scalar const_coeff,
                               CubicSpline* c_spline = HERMES_DEFAULT_SPLINE, GeomType gt = HERMES_PLANAR)
        : WeakForm::VectorFormVol(i, areas), const_coeff(const_coeff), spline_coeff(c_spline), gt(gt)
      {
        // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
        if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
      }

      ~DefaultResidualDiffusion() {
        if (spline_coeff != HERMES_DEFAULT_SPLINE) delete spline_coeff;
      };

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
                           Geom<double> *e, ExtData<scalar> *ext) const
      {
        scalar result = 0;
        if (gt == HERMES_PLANAR) {
          for (int i = 0; i < n; i++) {
            result += wt[i] * const_coeff * spline_coeff->get_value(u_ext[0]->val[i])
                      * (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i]);
          }
        }
        else {
          if (gt == HERMES_AXISYM_X) {
            for (int i = 0; i < n; i++) {
                    result += wt[i] * e->y[i] * const_coeff * spline_coeff->get_value(u_ext[0]->val[i])
                              * (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i]);
            }
          }
          else {
            for (int i = 0; i < n; i++) {
                    result += wt[i] * e->x[i] * const_coeff * spline_coeff->get_value(u_ext[0]->val[i])
                              * (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i]);
            }
          }
        }

        return result;
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                      Geom<Ord> *e, ExtData<Ord> *ext) const {
        Ord result = 0;
        // Planar base.
        for (int i = 0; i < n; i++) {
          result += wt[i] * const_coeff * spline_coeff->get_value(u_ext[0]->val[i])
                   * (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i]);
        }
        if (gt != HERMES_PLANAR) result = result * Ord(1);

        return result;
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::VectorFormVol* clone() {
        return new DefaultResidualDiffusion(*this);
      }

      private:
        scalar const_coeff;
        CubicSpline* spline_coeff;
        GeomType gt;
    };

    /* Default volumetric vector form \int_{area} spline_coeff1(u_ext[0]) * u->dx * v->val
       + spline_coeff2(u_ext[0]) * u->dy * v->val d\bfx
       spline_coeff1, spline_coeff2... non-constant parameters given by cubic splines
    */

    class DefaultResidualAdvection : public WeakForm::VectorFormVol
    {
    public:
      DefaultResidualAdvection(int i, std::string area = HERMES_ANY, scalar const_coeff = 1.0,
                               CubicSpline* c_spline1 = HERMES_DEFAULT_SPLINE,
                               CubicSpline* c_spline2 = HERMES_DEFAULT_SPLINE, 
                               GeomType gt = HERMES_PLANAR)
        : WeakForm::VectorFormVol(i, area), const_coeff(const_coeff), spline_coeff1(c_spline1), 
          spline_coeff2(c_spline2), gt(gt)
      {
        if (gt != HERMES_PLANAR) error("Axisymmetric advection forms not implemented yet.");

        // If spline1 == HERMES_DEFAULT_SPLINE or spline2 == HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
        if (c_spline1 == HERMES_DEFAULT_SPLINE) this->spline_coeff1 = new CubicSpline(1.0);
        if (c_spline2 == HERMES_DEFAULT_SPLINE) this->spline_coeff2 = new CubicSpline(1.0);
      }
      DefaultResidualAdvection(int i, Hermes::vector<std::string> areas, scalar const_coeff = 1.0,
                               CubicSpline* c_spline1 = HERMES_DEFAULT_SPLINE,
                               CubicSpline* c_spline2 = HERMES_DEFAULT_SPLINE, GeomType gt = HERMES_PLANAR)
        : WeakForm::VectorFormVol(i, areas), const_coeff(const_coeff), spline_coeff1(spline_coeff1),
                                             spline_coeff2(spline_coeff2), gt(gt)
      {
        if (gt != HERMES_PLANAR) error("Axisymmetric advection forms not implemented yet.");

        // If spline1 == HERMES_DEFAULT_SPLINE or spline2 == HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
        if (c_spline1 == HERMES_DEFAULT_SPLINE) this->spline_coeff1 = new CubicSpline(1.0);
        if (c_spline2 == HERMES_DEFAULT_SPLINE) this->spline_coeff2 = new CubicSpline(1.0);
      }

      ~DefaultResidualAdvection() {
        if (spline_coeff1 != HERMES_DEFAULT_SPLINE) delete spline_coeff1;
        if (spline_coeff2 != HERMES_DEFAULT_SPLINE) delete spline_coeff2;
      };

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
        return new DefaultResidualAdvection(*this);
      }

      private:
        scalar const_coeff;
        CubicSpline* spline_coeff1, *spline_coeff2;
        GeomType gt;
    };
  }

  namespace SurfaceMatrixForms {

    /* Default surface matrix form \int_{area} const_coeff * function_coeff(x, y) * u * v dS
       const_coeff... constant number
       function_coeff... (generally nonconstant) function of x, y
    */

    class DefaultMatrixFormSurf : public WeakForm::MatrixFormSurf
    {
    public:
      DefaultMatrixFormSurf(int i, int j, std::string area = HERMES_ANY,
                            scalar const_coeff = 1.0, DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                            GeomType gt = HERMES_PLANAR)
        : WeakForm::MatrixFormSurf(i, j, area), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
      {
        // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
        if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
      }
      DefaultMatrixFormSurf(int i, int j, Hermes::vector<std::string> areas,
                            scalar const_coeff, DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                            GeomType gt = HERMES_PLANAR)
        : WeakForm::MatrixFormSurf(i, j, areas), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
      {
        // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
        if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
      }

      ~DefaultMatrixFormSurf() {
        if (function_coeff != HERMES_DEFAULT_FUNCTION) delete function_coeff;
      };

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                           Geom<double> *e, ExtData<scalar> *ext) const {
        scalar result = 0;
        if (gt == HERMES_PLANAR) {
    for (int i = 0; i < n; i++) {
      result += wt[i] * function_coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
    }
        }
        else {
          if (gt == HERMES_AXISYM_X) {
      for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * function_coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
      }
          }
          else {
      for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * function_coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
      }
          }
        }

  return const_coeff * result;
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                      Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
        Ord result = 0;
        if (gt == HERMES_PLANAR) {
    for (int i = 0; i < n; i++) {
      result += wt[i] * function_coeff->ord(e->x[i], e->y[i]) * u->val[i] * v->val[i];
    }
        }
        else {
          if (gt == HERMES_AXISYM_X) {
      for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * function_coeff->ord(e->x[i], e->y[i]) * u->val[i] * v->val[i];
      }
          }
          else {
      for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * function_coeff->ord(e->x[i], e->y[i]) * u->val[i] * v->val[i];
      }
          }
        }

  return result;
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::MatrixFormSurf* clone() {
        return new DefaultMatrixFormSurf(*this);
      }

      private:
        scalar const_coeff;
        DefaultFunction* function_coeff;
        GeomType gt;
    };

    /* Default surface matrix form \int_{area} const_coeff * spline_coeff'(u_ext[0]) * u_ext[0] * u * v
       + const_coeff * spline_coeff(u_ext[0]) * u * v dS
       spline_coeff... non-constant parameter given by a spline
    */

    class DefaultJacobianFormSurf : public WeakForm::MatrixFormSurf
    {
    public:
      DefaultJacobianFormSurf(int i, int j, std::string area = HERMES_ANY, scalar const_coeff = 1.0,
                              CubicSpline* c_spline = HERMES_DEFAULT_SPLINE,
                              GeomType gt = HERMES_PLANAR)
        : WeakForm::MatrixFormSurf(i, j, area), const_coeff(const_coeff), spline_coeff(c_spline), gt(gt)
      {
        // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
        if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
      }
      DefaultJacobianFormSurf(int i, int j, Hermes::vector<std::string> areas, scalar const_coeff = 1.0,
                              CubicSpline* c_spline = HERMES_DEFAULT_SPLINE,
                              GeomType gt = HERMES_PLANAR)
        : WeakForm::MatrixFormSurf(i, j, areas), const_coeff(const_coeff), spline_coeff(spline_coeff), gt(gt)
      {
        // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
        if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
      }

      ~DefaultJacobianFormSurf() {
        if (spline_coeff != HERMES_DEFAULT_SPLINE) delete spline_coeff;
      };

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
        scalar const_coeff;
        CubicSpline* spline_coeff;
        GeomType gt;
    };
  }

  namespace SurfaceVectorForms {

    /* Default surface vector form \int_{area} const_coeff * function_coeff(x, y) * v dS
       const_coeff... constant number
       function_coeff... (generally nonconstant) function of x, y
    */

    class DefaultVectorFormSurf : public WeakForm::VectorFormSurf
    {
    public:
      DefaultVectorFormSurf(int i, std::string area = HERMES_ANY, scalar const_coeff = 1.0,
                            DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                            GeomType gt = HERMES_PLANAR)
        : WeakForm::VectorFormSurf(i, area), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
      {
        // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
        if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
      }
      DefaultVectorFormSurf(int i, Hermes::vector<std::string> areas, scalar const_coeff = 1.0,
                            DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                            GeomType gt = HERMES_PLANAR)
        : WeakForm::VectorFormSurf(i, areas), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
      {
        // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
        if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
      }

      ~DefaultVectorFormSurf() {
        if (function_coeff != HERMES_DEFAULT_FUNCTION) delete function_coeff;
      };

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                           Geom<double> *e, ExtData<scalar> *ext) const {
        scalar result = 0;
        if (gt == HERMES_PLANAR) {
    for (int i = 0; i < n; i++) {
      result += wt[i] * function_coeff->value(e->x[i], e->y[i]) * v->val[i];
    }
        }
        else {
          if (gt == HERMES_AXISYM_X) {
      for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * function_coeff->value(e->x[i], e->y[i]) * v->val[i];
      }
          }
          else {
      for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * function_coeff->value(e->x[i], e->y[i]) * v->val[i];
      }
          }
        }

  return const_coeff * result;
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                      Geom<Ord> *e, ExtData<Ord> *ext) const {
        Ord result = 0;
        if (gt == HERMES_PLANAR) {
    for (int i = 0; i < n; i++) {
      result += wt[i] * function_coeff->ord(e->x[i], e->y[i]) * v->val[i];
    }
        }
        else {
          if (gt == HERMES_AXISYM_X) {
      for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * function_coeff->ord(e->x[i], e->y[i]) * v->val[i];
      }
          }
          else {
      for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * function_coeff->ord(e->x[i], e->y[i]) * v->val[i];
      }
          }
        }

  return result;
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::VectorFormSurf* clone() {
        return new DefaultVectorFormSurf(*this);
      }

      private:
        scalar const_coeff;
        DefaultFunction* function_coeff;
        GeomType gt;
    };

    class DefaultMultiComponentVectorFormSurf : public WeakForm::MultiComponentVectorFormSurf
    {
    public:
      DefaultMultiComponentVectorFormSurf(Hermes::vector<unsigned int> coordinates, 
                                          std::string area = HERMES_ANY,
                                          Hermes::vector<scalar> coeffs = Hermes::vector<scalar>(1.0), 
                                          GeomType gt = HERMES_PLANAR)
      : WeakForm::MultiComponentVectorFormSurf(coordinates, area), coeffs(coeffs), gt(gt) { }
      DefaultMultiComponentVectorFormSurf(Hermes::vector<unsigned int> coordinates, 
                                          Hermes::vector<std::string> areas,
                                          Hermes::vector<scalar> coeffs, GeomType gt = HERMES_PLANAR)
      : WeakForm::MultiComponentVectorFormSurf(coordinates, areas), coeffs(coeffs), gt(gt) { }

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
        return new DefaultMultiComponentVectorFormSurf(*this);
      }
      */

      private:
        Hermes::vector<scalar> coeffs;
        GeomType gt;
    };

    /* Default surface vector form \int_{area} const_coeff * function_coeff(x, y) * u_ext[0] v dS
       const_coeff... constant number
       function_coeff... (generally nonconstant) function of x, y
    */

    class DefaultResidualSurf : public WeakForm::VectorFormSurf
    {
    public:
      DefaultResidualSurf(int i, std::string area = HERMES_ANY, scalar const_coeff = 1.0,
                          DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                          GeomType gt = HERMES_PLANAR)
        : WeakForm::VectorFormSurf(i, area), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
      {
        // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
        if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
      }
      DefaultResidualSurf(int i, Hermes::vector<std::string> areas, scalar const_coeff = 1.0,
                          DefaultFunction* f_coeff = HERMES_DEFAULT_FUNCTION,
                          GeomType gt = HERMES_PLANAR)
        : WeakForm::VectorFormSurf(i, areas), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
      {
        // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
        if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
      }

      ~DefaultResidualSurf() {
        if (function_coeff != HERMES_DEFAULT_FUNCTION) delete function_coeff;
      };

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                           Geom<double> *e, ExtData<scalar> *ext) const {
        scalar result = 0;
        if (gt == HERMES_PLANAR) {
    for (int i = 0; i < n; i++) {
      result += wt[i] * function_coeff->value(e->x[i], e->y[i]) * u_ext[0]->val[i] * v->val[i];
    }
        }
        else {
          if (gt == HERMES_AXISYM_X) {
      for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * function_coeff->value(e->x[i], e->y[i]) * u_ext[0]->val[i] * v->val[i];
      }
          }
          else {
      for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * function_coeff->value(e->x[i], e->y[i]) * u_ext[0]->val[i] * v->val[i];
      }
          }
        }

  return const_coeff * result;
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                      Geom<Ord> *e, ExtData<Ord> *ext) const {
        Ord result = 0;
        if (gt == HERMES_PLANAR) {
    for (int i = 0; i < n; i++) {
      result += wt[i] * function_coeff->ord(e->x[i], e->y[i]) * u_ext[0]->val[i] * v->val[i];
    }
        }
        else {
          if (gt == HERMES_AXISYM_X) {
      for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * function_coeff->ord(e->x[i], e->y[i]) * u_ext[0]->val[i] * v->val[i];
      }
          }
          else {
      for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * function_coeff->ord(e->x[i], e->y[i]) * u_ext[0]->val[i] * v->val[i];
      }
          }
        }

  return result;
      }

      // This is to make the form usable in rk_time_step().
      virtual WeakForm::VectorFormSurf* clone() {
        return new DefaultResidualSurf(*this);
      }

      private:
        scalar const_coeff;
        DefaultFunction* function_coeff;
        GeomType gt;
    };
  }

  /* Default weak form for the Laplace equation -div(S(u) grad u) = 0. */

  class DefaultWeakFormLaplace : public WeakForm
  {
  public:
    DefaultWeakFormLaplace(std::string area = HERMES_ANY, scalar const_coeff = 1.0,
                           CubicSpline* c_spline = HERMES_DEFAULT_SPLINE,
                           GeomType gt = HERMES_PLANAR) : WeakForm()
    {
      // Jacobian.
      add_matrix_form(new VolumetricMatrixForms::DefaultJacobianDiffusion(0, 0, area, const_coeff,
                                                                          HERMES_DEFAULT_SPLINE, HERMES_SYM, gt));
      // Residual.
      add_vector_form(new VolumetricVectorForms::DefaultResidualDiffusion(0, area, const_coeff,
                                                                          HERMES_DEFAULT_SPLINE, gt));
    };
  };


  /* Default weak form for the Poisson equation -div(S(u) grad u) - rhs = 0. */

  class DefaultWeakFormPoisson : public WeakForm
  {
  public:
  DefaultWeakFormPoisson(DefaultFunction* rhs = HERMES_DEFAULT_FUNCTION,
                         std::string area = HERMES_ANY, scalar const_coeff = 1.0, 
                         CubicSpline* c_spline = HERMES_DEFAULT_SPLINE,
                         GeomType gt = HERMES_PLANAR) : WeakForm()
    {
      // Jacobian.
      add_matrix_form(new VolumetricMatrixForms::DefaultJacobianDiffusion(0, 0, area, const_coeff,
                                                                          c_spline, HERMES_SYM, gt));
      // Residual.
      add_vector_form(new VolumetricVectorForms::DefaultResidualDiffusion(0, area, const_coeff,
                                                                          c_spline, gt));
      add_vector_form(new VolumetricVectorForms::DefaultVectorFormVol(0, HERMES_ANY, -1.0, rhs, gt));
    };
  };

}

#endif
