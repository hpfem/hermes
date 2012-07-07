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

#ifndef __H2D_INTEGRALS_H1_H
#define __H2D_INTEGRALS_H1_H

#include "../quadrature/limit_order.h"
#include "../forms.h"
#include "../function/function.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Real>
    Real int_v(int n, double *wt, Func<Real> *v)
    {
      Real result = Real(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * v->val[i];
      return result;
    }

    template<typename Real>
    Real int_x_v(int n, double *wt, Func<Real> *v, Geom<Real> *e)
    {
      Real result = Real(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * e->x[i] * v->val[i];
      return result;
    }

    template<typename Real>
    Real int_y_v(int n, double *wt, Func<Real> *v, Geom<Real> *e)
    {
      Real result = Real(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * e->y[i] * v->val[i];
      return result;
    }

    template<typename Real, typename Scalar>
    Scalar int_u_v(int n, double *wt, Func<Real> *u, Func<Real> *v)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * u->val[i] * v->val[i];
      return result;
    }

    // For residual forms.
    template<typename Real, typename Scalar>
    Scalar int_u_ext_v(int n, double *wt, Func<Scalar> *u_ext, Func<Real> *v)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * u_ext->val[i] * v->val[i];
      return result;
    }

    template<typename Real, typename Scalar>
    Scalar int_x_u_v(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * e->x[i] * u->val[i] * v->val[i];
      return result;
    }

    // For residual forms.
    template<typename Real, typename Scalar>
    Scalar int_x_u_ext_v(int n, double *wt, Func<Scalar> *u_ext, Func<Real> *v, Geom<Real> *e)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * e->x[i] * u_ext->val[i] * v->val[i];
      return result;
    }

    template<typename Real, typename Scalar>
    Scalar int_y_u_v(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * e->y[i] * u->val[i] * v->val[i];
      return result;
    }

    // For residual forms.
    template<typename Real, typename Scalar>
    Scalar int_y_u_ext_v(int n, double *wt, Func<Scalar> *u_ext, Func<Real> *v, Geom<Real> *e)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * e->y[i] * u_ext->val[i] * v->val[i];
      return result;
    }

    template<typename Real, typename Scalar>
    Scalar int_F_v(int n, double *wt, Real (*F)(Real x, Real y), Func<Real> *v, Geom<Real> *e)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * ((*F)(e->x[i], e->y[i]) * v->val[i]);
      return result;
    }

    template<typename Real, typename Scalar>
    Scalar int_grad_u_grad_v(int n, double *wt, Func<Real> *u, Func<Real> *v)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
      return result;
    }

    // For residual forms.
    template<typename Real, typename Scalar>
    Scalar int_grad_u_ext_grad_v(int n, double *wt, Func<Scalar> *u_ext, Func<Real> *v)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (u_ext->dx[i] * v->dx[i] + u_ext->dy[i] * v->dy[i]);
      return result;
    }

    template<typename Real, typename Scalar>
    Scalar int_x_grad_u_grad_v(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * e->x[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
      return result;
    }

    // For residual forms.
    template<typename Real, typename Scalar>
    Scalar int_x_grad_u_ext_grad_v(int n, double *wt, Func<Scalar> *u_ext, Func<Real> *v, Geom<Real> *e)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * e->x[i] * (u_ext->dx[i] * v->dx[i] + u_ext->dy[i] * v->dy[i]);
      return result;
    }

    template<typename Real, typename Scalar>
    Scalar int_y_grad_u_grad_v(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * e->y[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
      return result;
    }

    // For residual forms.
    template<typename Real, typename Scalar>
    Scalar int_y_grad_u_ext_grad_v(int n, double *wt, Func<Scalar> *u_ext, Func<Real> *v, Geom<Real> *e)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * e->y[i] * (u_ext->dx[i] * v->dx[i] + u_ext->dy[i] * v->dy[i]);
      return result;
    }

    template<typename Real, typename Scalar>
    Scalar int_dudx_v(int n, double *wt, Func<Real> *u, Func<Real> *v)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->dx[i] * v->val[i]);
      return result;
    }

    template<typename Real, typename Scalar>
    Scalar int_dudy_v(int n, double *wt, Func<Real> *u, Func<Real> *v)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->dy[i] * v->val[i]);
      return result;
    }

    template<typename Real, typename Scalar>
    Scalar int_u_dvdx(int n, double *wt, Func<Real> *u, Func<Real> *v)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * u->val[i] * v->dx[i];
      return result;
    }

    template<typename Real, typename Scalar>
    Scalar int_u_dvdx_over_x(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * ((e->x[i] > 0) ? u->val[i] * v->dx[i] / e->x[i] : 0.0);
      return result;
    }

    // For residual forms.
    template<typename Real, typename Scalar>
    Scalar int_u_ext_dvdx_over_x(int n, double *wt, Func<Scalar> *u_ext, Func<Real> *v, Geom<Real> *e)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * u_ext->val[i] * v->dx[i] / e->x[i];
      return result;
    }

    template<typename Real, typename Scalar>
    Scalar int_u_dvdy(int n, double *wt, Func<Real> *u, Func<Real> *v)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * u->val[i] * v->dy[i];
      return result;
    }

    template<typename Real, typename Scalar>
    Scalar int_u_dvdy_over_y(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * ((e->y[i] > 0) ? u->val[i] * v->dy[i] / e->y[i] : 0.0);
      return result;
    }

    // For residual forms.
    template<typename Real, typename Scalar>
    Scalar int_u_ext_dvdy_over_y(int n, double *wt, Func<Scalar> *u_ext, Func<Real> *v, Geom<Real> *e)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * u_ext->val[i] * v->dy[i] / e->y[i];
      return result;
    }

    template<typename Real, typename Scalar>
    Scalar int_dudx_dvdx(int n, double *wt, Func<Real> *u, Func<Real> *v)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->dx[i] * v->dx[i]);
      return result;
    }

    template<typename Real, typename Scalar>
    Scalar int_dudy_dvdy(int n, double *wt, Func<Real> *u, Func<Real> *v)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->dy[i] * v->dy[i]);
      return result;
    }

    template<typename Real, typename Scalar>
    Scalar int_dudx_dvdy(int n, double *wt, Func<Real> *u, Func<Real> *v)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->dx[i] * v->dy[i]);
      return result;
    }

    template<typename Real, typename Scalar>
    Scalar int_dudy_dvdx(int n, double *wt, Func<Real> *u, Func<Real> *v)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (v->dx[i] * u->dy[i]);
      return result;
    }

    template<typename Real, typename Scalar>
    Scalar int_w_nabla_u_v(int n, double *wt, Func<Real> *w1, Func<Real> *w2,
      Func<Real> *u, Func<Real> *v)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (w1->val[i] * u->dx[i] + w2->val[i] * u->dy[i]) * v->val[i];
      return result;
    }

    //// error & norm integrals  ////////////////////////////////////////////////////////////////////////

    // the inner integration loops for both constant and non-constant jacobian elements
    // for expression without partial derivatives - the variables e, quad, o must be already
    // defined and initialized
#define h1_integrate_expression(exp) \
    {double3* pt = quad->get_points(o, ru->get_active_element()->get_mode()); \
    int np = quad->get_num_points(o, ru->get_active_element()->get_mode()); \
    if(ru->is_jacobian_const()){ \
    for (int i = 0; i < np; i++) \
    result += pt[i][2] * (exp); \
    result *= ru->get_const_jacobian(); \
    } \
  else { \
  double* jac = ru->get_jacobian(o); \
  for (int i = 0; i < np; i++) \
  result += pt[i][2] * jac[i] * (exp); \
    }}

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    template<typename Scalar>
    inline double int_h1_error(Function<Scalar>* fu, Function<Scalar>* fv, RefMap* ru, RefMap* rv)
    {
      Quad2D* quad = fu->get_quad_2d();
      assert(quad == fv->get_quad_2d());

      int o = std::max(2*fu->get_fn_order(), 2*fv->get_fn_order()) + ru->get_inv_ref_order();
      limit_order(o, ru->get_active_element()->get_mode());
      fu->set_quad_order(o);
      fv->set_quad_order(o);

      Scalar* fnu = fu->get_fn_values();
      Scalar* fnv = fv->get_fn_values();

      Scalar *dudx, *dudy, *dvdx, *dvdy;
      fu->get_dx_dy_values(dudx, dudy);
      fv->get_dx_dy_values(dvdx, dvdy);

      double result = 0.0;
      h1_integrate_expression(sqr(fnu[i] - fnv[i]) + sqr(dudx[i] - dvdx[i]) + sqr(dudy[i] - dvdy[i]));
      return result;
    }

    template<typename Scalar>
    inline double int_h1_semi_error(Function<Scalar>* fu, Function<Scalar>* fv, RefMap* ru, RefMap* rv)
    {
      Quad2D* quad = fu->get_quad_2d();
      assert(quad == fv->get_quad_2d());

      int o = std::max(2*fu->get_fn_order(), 2*fv->get_fn_order()) + ru->get_inv_ref_order();
      limit_order(o, ru->get_active_element()->get_mode());
      fu->set_quad_order(o);
      fv->set_quad_order(o);

      Scalar* fnu = fu->get_fn_values();
      Scalar* fnv = fv->get_fn_values();

      Scalar *dudx, *dudy, *dvdx, *dvdy;
      fu->get_dx_dy_values(dudx, dudy);
      fv->get_dx_dy_values(dvdx, dvdy);

      double result = 0.0;
      h1_integrate_expression(sqr(dudx[i] - dvdx[i]) + sqr(dudy[i] - dvdy[i]));
      return result;
    }

    template<typename Scalar>
    inline double int_l2_error(Function<Scalar>* fu, Function<Scalar>* fv, RefMap* ru, RefMap* rv)
    {
      Quad2D* quad = fu->get_quad_2d();
      assert(quad == fv->get_quad_2d());

      int o = std::max(2*fu->get_fn_order(), 2*fv->get_fn_order()) + ru->get_inv_ref_order();
      limit_order(o, ru->get_active_element()->get_mode());
      fu->set_quad_order(o, H2D_FN_VAL);
      fv->set_quad_order(o, H2D_FN_VAL);

      Scalar* fnu = fu->get_fn_values();
      Scalar* fnv = fv->get_fn_values();

      double result = 0.0;
      h1_integrate_expression(sqr(fnu[i] - fnv[i]));
      return result;
    }

    template<typename Scalar>
    inline double int_dx_error(Function<Scalar>* fu, Function<Scalar>* fv, RefMap* ru, RefMap* rv)
    {
      Quad2D* quad = fu->get_quad_2d();
      assert(quad == fv->get_quad_2d());

      int o = std::max(2*fu->get_fn_order(), 2*fv->get_fn_order()) + ru->get_inv_ref_order();
      limit_order(o, ru->get_active_element()->get_mode());
      fu->set_quad_order(o);
      fv->set_quad_order(o);

      Scalar *dudx, *dudy, *dvdx, *dvdy;
      fu->get_dx_dy_values(dudx, dudy);
      fv->get_dx_dy_values(dvdx, dvdy);

      double result = 0.0;
      h1_integrate_expression(sqr(dudx[i] - dvdx[i]));
      return result;
    }

    template<typename Scalar>
    inline double int_dy_error(Function<Scalar>* fu, Function<Scalar>* fv, RefMap* ru, RefMap* rv)
    {
      Quad2D* quad = fu->get_quad_2d();
      assert(quad == fv->get_quad_2d());

      int o = std::max(2*fu->get_fn_order(), 2*fv->get_fn_order()) + ru->get_inv_ref_order();
      limit_order(o, ru->get_active_element()->get_mode());
      fu->set_quad_order(o);
      fv->set_quad_order(o);

      Scalar *dudx, *dudy, *dvdx, *dvdy;
      fu->get_dx_dy_values(dudx, dudy);
      fv->get_dx_dy_values(dvdx, dvdy);

      double result = 0.0;
      h1_integrate_expression(sqr(dudy[i] - dvdy[i]));
      return result;
    }

    template<typename Scalar>
    inline double int_h1_norm(Function<Scalar>* fu, RefMap* ru)
    {
      Quad2D* quad = fu->get_quad_2d();

      int o = 2*fu->get_fn_order() + ru->get_inv_ref_order();
      limit_order(o, ru->get_active_element()->get_mode());
      fu->set_quad_order(o);

      Scalar* fnu = fu->get_fn_values();
      Scalar *dudx, *dudy;
      fu->get_dx_dy_values(dudx, dudy);

      double result = 0.0;
      h1_integrate_expression(sqr(fnu[i]) + sqr(dudx[i]) + sqr(dudy[i]));
      return result;
    }

    template<typename Scalar>
    inline double int_h1_seminorm(Function<Scalar>* fu, RefMap* ru)
    {
      Quad2D* quad = fu->get_quad_2d();

      int o = 2*fu->get_fn_order() + ru->get_inv_ref_order();
      limit_order(o, ru->get_active_element()->get_mode());
      fu->set_quad_order(o);

      Scalar* fnu = fu->get_fn_values();
      Scalar *dudx, *dudy;
      fu->get_dx_dy_values(dudx, dudy);

      double result = 0.0;
      h1_integrate_expression(sqr(dudx[i]) + sqr(dudy[i]));
      return result;
    }

    template<typename Scalar>
    inline double int_l2_norm(Function<Scalar>* fu, RefMap* ru)
    {
      Quad2D* quad = fu->get_quad_2d();

      int o = 2*fu->get_fn_order() + ru->get_inv_ref_order();
      limit_order(o, ru->get_active_element()->get_mode());
      fu->set_quad_order(o, H2D_FN_VAL);

      Scalar* fnu = fu->get_fn_values();

      double result = 0.0;
      h1_integrate_expression(sqr(fnu[i]));
      return result;
    }
  }
}
#endif