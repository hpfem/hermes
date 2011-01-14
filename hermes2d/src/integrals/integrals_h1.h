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
#include "../weakform/weakform.h"

//// the following integrals can be used in both volume and surface forms ////

template<typename Real, typename Scalar>
Scalar int_v(int n, double *wt, Func<Real> *v)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (v->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar int_u_v(int n, double *wt, Func<Real> *u, Func<Real> *v)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar int_F_v(int n, double *wt, Real (*F)(Real x, Real y), Func<Real> *v, Geom<Real> *e)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ((*F)(e->x[i], e->y[i]) * v->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar int_grad_u_grad_v(int n, double *wt, Func<Real> *u, Func<Real> *v)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar int_dudx_v(int n, double *wt, Func<Real> *u, Func<Real> *v)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->dx[i] * v->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar int_dudy_v(int n, double *wt, Func<Real> *u, Func<Real> *v)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->dy[i] * v->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar int_u_dvdx(int n, double *wt, Func<Real> *u, Func<Real> *v)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (v->dx[i] * u->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar int_u_dvdy(int n, double *wt, Func<Real> *u, Func<Real> *v)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (v->dy[i] * u->val[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar int_dudx_dvdx(int n, double *wt, Func<Real> *u, Func<Real> *v)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->dx[i] * v->dx[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar int_dudy_dvdy(int n, double *wt, Func<Real> *u, Func<Real> *v)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->dy[i] * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar int_dudx_dvdy(int n, double *wt, Func<Real> *u, Func<Real> *v)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->dx[i] * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar int_dudy_dvdx(int n, double *wt, Func<Real> *u, Func<Real> *v)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (v->dx[i] * u->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar int_w_nabla_u_v(int n, double *wt, Func<Real> *w1, Func<Real> *w2, 
                       Func<Real> *u, Func<Real> *v)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (w1->val[i] * u->dx[i] + w2->val[i] * u->dy[i]) * v->val[i];
  return result;
}

//// error calculation for adaptivity  ////

template<typename Real, typename Scalar>
Scalar h1_error_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Scalar> *u, 
               Func<Scalar> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * conj(v->val[i]) + u->dx[i] * conj(v->dx[i]) 
                       + u->dy[i] * conj(v->dy[i]));
  return result;
}

template<typename Real, typename Scalar>
Scalar h1_error_semi_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Scalar> *u, 
                    Func<Scalar> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->dx[i] * conj(v->dx[i]) + u->dy[i] * conj(v->dy[i]));
  return result;
}

//// error & norm integrals  ////////////////////////////////////////////////////////////////////////

// the inner integration loops for both constant and non-constant jacobian elements
// for expression without partial derivatives - the variables e, quad, o must be already
// defined and initialized
#define h1_integrate_expression(exp) \
  {double3* pt = quad->get_points(o); \
  int np = quad->get_num_points(o); \
  if (ru->is_jacobian_const()){ \
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

template<typename T>
inline double int_h1_error(Function<T>* fu, Function<T>* fv, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = fu->get_quad_2d();
  assert(quad == fv->get_quad_2d());

  int o = std::max(2*fu->get_fn_order(), 2*fv->get_fn_order()) + ru->get_inv_ref_order();
  limit_order(o);
  fu->set_quad_order(o);
  fv->set_quad_order(o);

  scalar* fnu = fu->get_fn_values();
  scalar* fnv = fv->get_fn_values();

  scalar *dudx, *dudy, *dvdx, *dvdy;
  fu->get_dx_dy_values(dudx, dudy);
  fv->get_dx_dy_values(dvdx, dvdy);

  double result = 0.0;
  h1_integrate_expression(sqr(fnu[i] - fnv[i]) + sqr(dudx[i] - dvdx[i]) + sqr(dudy[i] - dvdy[i]));
  return result;
}


template<typename T>
inline double int_h1_semi_error(Function<T>* fu, Function<T>* fv, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = fu->get_quad_2d();
  assert(quad == fv->get_quad_2d());

  int o = std::max(2*fu->get_fn_order(), 2*fv->get_fn_order()) + ru->get_inv_ref_order();
  limit_order(o);
  fu->set_quad_order(o);
  fv->set_quad_order(o);

  scalar* fnu = fu->get_fn_values();
  scalar* fnv = fv->get_fn_values();

  scalar *dudx, *dudy, *dvdx, *dvdy;
  fu->get_dx_dy_values(dudx, dudy);
  fv->get_dx_dy_values(dvdx, dvdy);

  double result = 0.0;
  h1_integrate_expression(sqr(dudx[i] - dvdx[i]) + sqr(dudy[i] - dvdy[i]));
  return result;
}


template<typename T>
inline double int_l2_error(Function<T>* fu, Function<T>* fv, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = fu->get_quad_2d();
  assert(quad == fv->get_quad_2d());

  int o = std::max(2*fu->get_fn_order(), 2*fv->get_fn_order()) + ru->get_inv_ref_order();
  limit_order(o);
  fu->set_quad_order(o, H2D_FN_VAL);
  fv->set_quad_order(o, H2D_FN_VAL);

  scalar* fnu = fu->get_fn_values();
  scalar* fnv = fv->get_fn_values();

  double result = 0.0;
  h1_integrate_expression(sqr(fnu[i] - fnv[i]));
  return result;
}


template<typename T>
inline double int_dx_error(Function<T>* fu, Function<T>* fv, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = fu->get_quad_2d();
  assert(quad == fv->get_quad_2d());

  int o = std::max(2*fu->get_fn_order(), 2*fv->get_fn_order()) + ru->get_inv_ref_order();
  limit_order(o);
  fu->set_quad_order(o);
  fv->set_quad_order(o);

  scalar *dudx, *dudy, *dvdx, *dvdy;
  fu->get_dx_dy_values(dudx, dudy);
  fv->get_dx_dy_values(dvdx, dvdy);

  double result = 0.0;
  h1_integrate_expression(sqr(dudx[i] - dvdx[i]));
  return result;
}


template<typename T>
inline double int_dy_error(Function<T>* fu, Function<T>* fv, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = fu->get_quad_2d();
  assert(quad == fv->get_quad_2d());

  int o = std::max(2*fu->get_fn_order(), 2*fv->get_fn_order()) + ru->get_inv_ref_order();
  limit_order(o);
  fu->set_quad_order(o);
  fv->set_quad_order(o);

  scalar *dudx, *dudy, *dvdx, *dvdy;
  fu->get_dx_dy_values(dudx, dudy);
  fv->get_dx_dy_values(dvdx, dvdy);

  double result = 0.0;
  h1_integrate_expression(sqr(dudy[i] - dvdy[i]));
  return result;
}


template<typename T>
inline double int_h1_norm(Function<T>* fu, RefMap* ru)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = 2*fu->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  fu->set_quad_order(o);

  scalar* fnu = fu->get_fn_values();
  scalar *dudx, *dudy;
  fu->get_dx_dy_values(dudx, dudy);

  double result = 0.0;
  h1_integrate_expression(sqr(fnu[i]) + sqr(dudx[i]) + sqr(dudy[i]));
  return result;
}


template<typename T>
inline double int_h1_seminorm(Function<T>* fu, RefMap* ru)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = 2*fu->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  fu->set_quad_order(o);

  scalar* fnu = fu->get_fn_values();
  scalar *dudx, *dudy;
  fu->get_dx_dy_values(dudx, dudy);

  double result = 0.0;
  h1_integrate_expression(sqr(dudx[i]) + sqr(dudy[i]));
  return result;
}


template<typename T>
inline double int_l2_norm(Function<T>* fu, RefMap* ru)
{
  Quad2D* quad = fu->get_quad_2d();

  int o = 2*fu->get_fn_order() + ru->get_inv_ref_order();
  limit_order(o);
  fu->set_quad_order(o, H2D_FN_VAL);

  scalar* fnu = fu->get_fn_values();

  double result = 0.0;
  h1_integrate_expression(sqr(fnu[i]));
  return result;
}


#endif
