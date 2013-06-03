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

#include "weakform_library/weakforms_maxwell.h"
#include "weakform_library/integrals_h1.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace WeakFormsMaxwell
    {
      template<typename Scalar>
      DefaultJacobianMagnetostatics<Scalar>::DefaultJacobianMagnetostatics(int i, int j, std::string area, Scalar const_coeff,
        CubicSpline* c_spline,
        SymFlag sym,
        GeomType gt,
        int order_increase)
        : MatrixFormVol<Scalar>(i, j), idx_j(j), const_coeff(const_coeff),
        spline_coeff(c_spline), gt(gt),
        order_increase(order_increase)
      {
        this->set_area(area);
        this->setSymFlag(sym);

        // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
        if(c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
      }
      template<typename Scalar>
      DefaultJacobianMagnetostatics<Scalar>::DefaultJacobianMagnetostatics(int i, int j, Hermes::vector<std::string> areas,
        Scalar const_coeff,
        CubicSpline* c_spline,
        SymFlag sym,
        GeomType gt,
        int order_increase)
        : MatrixFormVol<Scalar>(i, j), idx_j(j), const_coeff(const_coeff),
        spline_coeff(c_spline), gt(gt),
        order_increase(order_increase)
      {
        this->set_areas(areas);
        this->setSymFlag(sym);

        // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
        if(c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
      }

      template<typename Scalar>
      Scalar DefaultJacobianMagnetostatics<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
        Func<double> *v, Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar planar_part = 0;
        Scalar axisym_part = 0;
        for (int i = 0; i < n; i++)
        {
          Scalar B_i = sqrt(sqr(u_ext[idx_j]->dx[i]) + sqr(u_ext[idx_j]->dy[i]));
          if(std::abs(B_i) > Hermes::epsilon)
          {
            planar_part += wt[i] * const_coeff*spline_coeff->derivative(B_i) / B_i
              * (u_ext[idx_j]->dx[i] * u->dx[i] + u_ext[idx_j]->dy[i] * u->dy[i])
              * (u_ext[idx_j]->dx[i] * v->dx[i] + u_ext[idx_j]->dy[i] * v->dy[i]);
            if(gt == HERMES_AXISYM_X)
            {
              axisym_part += wt[i] * const_coeff*spline_coeff->derivative(B_i) / B_i / e->y[i]
              * (u_ext[idx_j]->dx[i] * u->dx[i] + u_ext[idx_j]->dy[i] * u->dy[i])
                * u_ext[idx_j]->val[i] * v->dy[i];
            }
            else if(gt == HERMES_AXISYM_Y)
            {
              axisym_part += wt[i] * const_coeff*spline_coeff->derivative(B_i) / B_i / e->x[i]
              * (u_ext[idx_j]->dx[i] * u->dx[i] + u_ext[idx_j]->dy[i] * u->dy[i])
                * u_ext[idx_j]->val[i] * v->dx[i];
            }
          }
          planar_part += wt[i] * const_coeff*spline_coeff->value(B_i)
            * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
          if(gt == HERMES_AXISYM_X)
          {
            axisym_part += wt[i] * const_coeff*spline_coeff->value(B_i) / e->y[i]
            * u->val[i] * v->dy[i];
          }
          else if(gt == HERMES_AXISYM_Y)
          {
            axisym_part += wt[i] * const_coeff*spline_coeff->value(B_i) / e->x[i]
            * u->val[i] * v->dx[i];
          }
        }

        return planar_part + axisym_part;
      }

      template<typename Scalar>
      Ord DefaultJacobianMagnetostatics<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
        Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord planar_part(0);
        for (int i = 0; i < n; i++)
        {
          Ord B_i = sqrt(sqr(u_ext[idx_j]->dx[i]) + sqr(u_ext[idx_j]->dy[i]));
          planar_part += wt[i] * const_coeff*spline_coeff->derivative(B_i) / B_i
            * (u_ext[idx_j]->dx[i] * u->dx[i] + u_ext[idx_j]->dy[i] * u->dy[i])
            * (u_ext[idx_j]->dx[i] * v->dx[i] + u_ext[idx_j]->dy[i] * v->dy[i]);
          planar_part += wt[i] * const_coeff*spline_coeff->value(B_i)
            * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
        }

        // This increase is for the axisymmetric part. We are not letting the
        // Ord class do it since it would automatically choose the highest order
        // due to the nonpolynomial 1/r term.
        return planar_part * Ord(order_increase);
      }

      template<typename Scalar>
      MatrixFormVol<Scalar>* DefaultJacobianMagnetostatics<Scalar>::clone() const
      {
        return new DefaultJacobianMagnetostatics(*this);
      }

      template<typename Scalar>
      DefaultResidualMagnetostatics<Scalar>::DefaultResidualMagnetostatics(int i, std::string area, Scalar const_coeff,
        CubicSpline* c_spline,
        GeomType gt,
        int order_increase)
        : VectorFormVol<Scalar>(i), idx_i(i), const_coeff(const_coeff), spline_coeff(c_spline),
        gt(gt), order_increase(order_increase)
      {
        this->set_area(area);
        // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
        if(c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
      }

      template<typename Scalar>
      DefaultResidualMagnetostatics<Scalar>::DefaultResidualMagnetostatics(int i, Hermes::vector<std::string> areas, Scalar const_coeff,
        CubicSpline* c_spline,
        GeomType gt, int order_increase)
        : VectorFormVol<Scalar>(i), idx_i(i), const_coeff(const_coeff), spline_coeff(c_spline), gt(gt),
        order_increase(order_increase)
      {
        this->set_areas(areas);

        // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
        if(c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
      }

      template<typename Scalar>
      Scalar DefaultResidualMagnetostatics<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar planar_part = 0;
        Scalar axisym_part = 0;
        for (int i = 0; i < n; i++)
        {
          Scalar B_i = sqrt(sqr(u_ext[idx_i]->dx[i]) + sqr(u_ext[idx_i]->dy[i]));
          planar_part += wt[i] * const_coeff*spline_coeff->value(B_i) *
            (u_ext[idx_i]->dx[i] * v->dx[i] + u_ext[idx_i]->dy[i] * v->dy[i]);
          if(gt == HERMES_AXISYM_X) axisym_part += wt[i] * const_coeff*spline_coeff->value(B_i) / e->y[i]
          * u_ext[idx_i]->val[i] * v->dy[i];
          else if(gt == HERMES_AXISYM_Y) axisym_part += wt[i] * const_coeff*spline_coeff->value(B_i) / e->x[i]
          * u_ext[idx_i]->val[i] * v->dx[i];
        }
        return planar_part + axisym_part;
      }

      template<typename Scalar>
      Ord DefaultResidualMagnetostatics<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
        Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord planar_part(0);
        for (int i = 0; i < n; i++)
        {
          Ord B_i = sqrt(sqr(u_ext[idx_i]->dx[i]) + sqr(u_ext[idx_i]->dy[i]));
          planar_part += wt[i] * const_coeff*spline_coeff->value(B_i) *
            (u_ext[idx_i]->dx[i] * v->dx[i] + u_ext[idx_i]->dy[i] * v->dy[i]);
        }
        return planar_part * Ord(order_increase);
      }

      // This is to make the form usable in rk_time_step_newton().
      template<typename Scalar>
      VectorFormVol<Scalar>* DefaultResidualMagnetostatics<Scalar>::clone() const
      {
        return new DefaultResidualMagnetostatics(*this);
      }

      template class HERMES_API DefaultJacobianMagnetostatics<double>;
      template class HERMES_API DefaultResidualMagnetostatics<double>;
    }
  }
}