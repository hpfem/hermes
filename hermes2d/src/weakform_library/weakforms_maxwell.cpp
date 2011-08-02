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

#include "weakforms_maxwell.h"
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
        : MatrixFormVol<Scalar>(i, j, area, sym), idx_j(j), const_coeff(const_coeff), 
        spline_coeff(c_spline), gt(gt),
        order_increase(order_increase) 
      {
        // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
        if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
      }
      template<typename Scalar>
      DefaultJacobianMagnetostatics<Scalar>::DefaultJacobianMagnetostatics(int i, int j, Hermes::vector<std::string> areas,
        Scalar const_coeff,
        CubicSpline* c_spline,
        SymFlag sym,
        GeomType gt,
        int order_increase)
        : MatrixFormVol<Scalar>(i, j, areas, sym), idx_j(j), const_coeff(const_coeff), 
        spline_coeff(c_spline), gt(gt),
        order_increase(order_increase) 
      {
        // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
        if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
      }

      template<typename Scalar>
      Scalar DefaultJacobianMagnetostatics<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
        Func<double> *v, Geom<double> *e, ExtData<Scalar> *ext) const 
      {
        Scalar planar_part = 0;
        Scalar axisym_part = 0;
        for (int i = 0; i < n; i++) 
        {
          Scalar B_i = sqrt(sqr(u_ext[idx_j]->dx[i]) + sqr(u_ext[idx_j]->dy[i]));
          if (std::abs(B_i) > 1e-12) 
          {
            planar_part += wt[i] * const_coeff*spline_coeff->derivative(B_i) / B_i
              * (u_ext[idx_j]->dx[i] * u->dx[i] + u_ext[idx_j]->dy[i] * u->dy[i])
              * (u_ext[idx_j]->dx[i] * v->dx[i] + u_ext[idx_j]->dy[i] * v->dy[i]);
            if (gt == HERMES_AXISYM_X) 
            {
              axisym_part += wt[i] * const_coeff*spline_coeff->derivative(B_i) / B_i / e->y[i]
              * (u_ext[idx_j]->dx[i] * u->dx[i] + u_ext[idx_j]->dy[i] * u->dy[i])
                * u_ext[idx_j]->val[i] * v->dy[i];
            }
            else if (gt == HERMES_AXISYM_Y) 
            {
              axisym_part += wt[i] * const_coeff*spline_coeff->derivative(B_i) / B_i / e->x[i]
              * (u_ext[idx_j]->dx[i] * u->dx[i] + u_ext[idx_j]->dy[i] * u->dy[i])
                * u_ext[idx_j]->val[i] * v->dx[i];
            }
          }
          planar_part += wt[i] * const_coeff*spline_coeff->value(B_i)
            * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
          if (gt == HERMES_AXISYM_X) 
          {
            axisym_part += wt[i] * const_coeff*spline_coeff->value(B_i) / e->y[i]
            * u->val[i] * v->dy[i];
          }
          else if (gt == HERMES_AXISYM_Y) 
          {
            axisym_part += wt[i] * const_coeff*spline_coeff->value(B_i) / e->x[i]
            * u->val[i] * v->dx[i];
          }
        }

        return planar_part + axisym_part;
      }

      template<typename Scalar>
      Ord DefaultJacobianMagnetostatics<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
        Geom<Ord> *e, ExtData<Ord> *ext) const 
      {
        Ord planar_part = 0;
        for (int i = 0; i < n; i++) 
        {
          Ord B_i = sqrt(sqr(u_ext[idx_j]->dx[i]) + sqr(u_ext[idx_j]->dy[i]));
          planar_part += wt[i] * const_coeff*spline_coeff->derivative_ord(B_i) / B_i
            * (u_ext[idx_j]->dx[i] * u->dx[i] + u_ext[idx_j]->dy[i] * u->dy[i])
            * (u_ext[idx_j]->dx[i] * v->dx[i] + u_ext[idx_j]->dy[i] * v->dy[i]);
          planar_part += wt[i] * const_coeff*spline_coeff->value_ord(B_i)
            * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
        }

        // This increase is for the axisymmetric part. We are not letting the
        // Ord class do it since it would automatically choose the highest order
        // due to the nonpolynomial 1/r term.
        return planar_part * Ord(order_increase);
      }

      template<typename Scalar>
      MatrixFormVol<Scalar>* DefaultJacobianMagnetostatics<Scalar>::clone() 
      {
        return new DefaultJacobianMagnetostatics(*this);
      }


      template<typename Scalar>
      DefaultResidualMagnetostatics<Scalar>::DefaultResidualMagnetostatics(int i, std::string area, Scalar const_coeff,
        CubicSpline* c_spline,
        GeomType gt,
        int order_increase)
        : VectorFormVol<Scalar>(i, area), idx_i(i), const_coeff(const_coeff), spline_coeff(c_spline), 
        gt(gt), order_increase(order_increase) 
      {
        // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
        if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
      }

      template<typename Scalar>
      DefaultResidualMagnetostatics<Scalar>::DefaultResidualMagnetostatics(int i, Hermes::vector<std::string> areas, Scalar const_coeff, 
        CubicSpline* c_spline,
        GeomType gt, int order_increase)
        : VectorFormVol<Scalar>(i, areas), idx_i(i), const_coeff(const_coeff), spline_coeff(c_spline), gt(gt),
        order_increase(order_increase) 
      {
        // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
        if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
      }

      template<typename Scalar>
      Scalar DefaultResidualMagnetostatics<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, ExtData<Scalar> *ext) const 
      {
        Scalar planar_part = 0;
        Scalar axisym_part = 0;
        for (int i = 0; i < n; i++) 
        {
          Scalar B_i = sqrt(sqr(u_ext[idx_i]->dx[i]) + sqr(u_ext[idx_i]->dy[i]));
          planar_part += wt[i] * const_coeff*spline_coeff->value(B_i) *
            (u_ext[idx_i]->dx[i] * v->dx[i] + u_ext[idx_i]->dy[i] * v->dy[i]);
          if (gt == HERMES_AXISYM_X) axisym_part += wt[i] * const_coeff*spline_coeff->value(B_i) / e->y[i]
          * u_ext[idx_i]->val[i] * v->dy[i];
          else if (gt == HERMES_AXISYM_Y) axisym_part += wt[i] * const_coeff*spline_coeff->value(B_i) / e->x[i]
          * u_ext[idx_i]->val[i] * v->dx[i];
        }
        return planar_part + axisym_part;
      }

      template<typename Scalar>
      Ord DefaultResidualMagnetostatics<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
        Geom<Ord> *e, ExtData<Ord> *ext) const 
      {
        Ord planar_part = 0;
        for (int i = 0; i < n; i++) 
        {
          Ord B_i = sqrt(sqr(u_ext[idx_i]->dx[i]) + sqr(u_ext[idx_i]->dy[i]));
          planar_part += wt[i] * const_coeff*spline_coeff->value_ord(B_i) *
            (u_ext[idx_i]->dx[i] * v->dx[i] + u_ext[idx_i]->dy[i] * v->dy[i]);
        }
        return planar_part * Ord(order_increase);
      }

      // This is to make the form usable in rk_time_step().
      template<typename Scalar>
      VectorFormVol<Scalar>* DefaultResidualMagnetostatics<Scalar>::clone() 
      {
        return new DefaultResidualMagnetostatics(*this);
      }
    };


    /* Default volumetric matrix form \int_{area} u->val[i] * ((vel_x - e->y[i] * vel_ang) *
    v->dx[i] + (vel_y + e->x[i] * vel_ang) * v->dy[i])
    vel_x, vel_y... velocity components
    vel_ang... velocity angle
    */

    /*
    template<typename Scalar>
    class DefaultLinearMagnetostaticsVelocity : public MatrixFormVol<Scalar>
    {
    public:
    DefaultLinearMagnetostaticsVelocity(int i, int j, double gamma, double vel_x, double vel_y, double vel_ang = 0.0)
    : MatrixFormVol<Scalar>(i, j), gamma(gamma), vel_x(vel_x), vel_y(vel_y), vel_ang(vel_ang) { }

    DefaultLinearMagnetostaticsVelocity(int i, int j, std::string area, double gamma, double vel_x, double vel_y, double vel_ang = 0.0)
    : MatrixFormVol<Scalar>(i, j, area, HERMES_NONSYM), gamma(gamma), vel_x(vel_x), vel_y(vel_y), vel_ang(vel_ang) { }

    Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
    Geom<double> *e, ExtData<Scalar> *ext) const 
    {
    Scalar result = 0;
    for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * ((vel_x - e->y[i] * vel_ang) * v->dx[i] +
    (vel_y + e->x[i] * vel_ang) * v->dy[i]);

    return -gamma * result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
    Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
    Ord result = 0;
    for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * (v->dx[i] + v->dy[i]);

    return result;
    }

    // This is to make the form usable in rk_time_step().
    MatrixFormVol<Scalar>* clone() 
    {
    return new DefaultLinearMagnetostaticsVelocity(*this);
    }

    private:
    double gamma, vel_x, vel_y, vel_ang;
    };
    */


    /* Default volumetric vector form \int_{area} coeff
    \nabla u_ext[idx_i] \cdot \nabla v d\bfx
    coeff... constant parameter
    */

    /*
    template<typename Scalar>
    class DefaultResidualLinearMagnetostatics : public VectorFormVol<Scalar>
    {
    public:
    DefaultResidualLinearMagnetostatics(int i, Scalar coeff, GeomType gt,
    int order_increase)
    : VectorFormVol<Scalar>(i), idx_i(i), coeff(coeff), gt(gt), order_increase(order_increase) { }
    DefaultResidualLinearMagnetostatics(int i, std::string area, Scalar coeff,
    GeomType gt, int order_increase)
    : VectorFormVol<Scalar>(i, area), idx_i(i), coeff(coeff), gt(gt), order_increase(order_increase) { }

    Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
    Geom<double> *e, ExtData<Scalar> *ext) const 
    {
    Scalar planar_part = int_grad_u_ext_grad_v<double, Scalar>(n, wt, u_ext[idx_i], v);
    Scalar axisym_part = 0;
    if (gt == HERMES_AXISYM_X)
    axisym_part = int_u_dvdy_over_y<double, Scalar>(n, wt, u_ext[idx_i], v, e);
    else if (gt == HERMES_AXISYM_Y)
    axisym_part = int_u_dvdx_over_x<double, Scalar>(n, wt, u_ext[idx_i], v, e);

    return coeff * (planar_part + axisym_part);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
    Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
    Ord planar_part = int_grad_u_ext_grad_v<Ord, Ord>(n, wt, u_ext[idx_i], v);
    return planar_part * Ord(order_increase);
    }

    // This is to make the form usable in rk_time_step().
    VectorFormVol<Scalar>* clone() 
    {
    return new DefaultResidualLinearMagnetostatics(*this);
    }

    private:
    int idx_i;
    Scalar coeff;
    GeomType gt;
    int order_increase;
    };
    */

    /* Default volumetric vector form \int_{area} spline_coeff(u_ext[0])
    \nabla u_ext[idx_i] \cdot \nabla v d\bfx
    spline_coeff... non-constant parameter given by a cubic spline
    */



    /* Default volumetric vector form \int_{area} rem/perm * (- sin(rem_ang / 180.0 * M_PI) *
    v->dx[i + cos(rem_ang / 180.0 * M_PI) * v->dy[i])
    rem... remanent induction
    rem_ang... remanent induction angle
    per... permeability
    */

    /*
    template<typename Scalar>
    class DefaultLinearMagnetostaticsRemanence : public VectorFormVol<Scalar>
    {
    public:
    DefaultLinearMagnetostaticsRemanence(int i, double perm, double rem, double rem_ang, GeomType gt)
    : VectorFormVol<Scalar>(i), perm(perm), rem(rem), rem_ang(rem_ang), gt(gt) { }

    DefaultLinearMagnetostaticsRemanence(int i, std::string area, double perm, double rem, double rem_ang, GeomType gt)
    : VectorFormVol<Scalar>(i, area), perm(perm), rem(rem), rem_ang(rem_ang), gt(gt) { }

    Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
    Geom<double> *e, ExtData<Scalar> *ext) const 
    {
    Scalar result = 0;
    for (int i = 0; i < n; i++)
    result += wt[i] * rem/perm * (- sin(rem_ang / 180.0 * M_PI) * v->dx[i]
    + cos(rem_ang / 180.0 * M_PI) * v->dy[i]);

    return (gt == HERMES_PLANAR ? -1.0 : 1.0) * result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
    Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
    Ord result = 0;
    for (int i = 0; i < n; i++)
    result += wt[i] * (v->dx[i] + v->dy[i]);

    return result;
    }

    // This is to make the form usable in rk_time_step().
    VectorFormVol<Scalar>* clone() 
    {
    return new DefaultLinearMagnetostaticsRemanence(*this);
    }

    private:
    double perm, rem, rem_ang;
    GeomType gt;
    };
    */
  }
}