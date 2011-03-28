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

#ifndef __H2D_INTEGRALS_HCURL_H
#define __H2D_INTEGRALS_HCURL_H

//// new volume integrals //////////////////////////////////////////////////////////////////////////////

template<typename real, typename scalar>
scalar int_e_f(int n, double *wt, Func<real> *u, Func<real> *v)
{
  scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
  return result;
}

template<typename real, typename scalar>
scalar int_curl_e_curl_f(int n, double *wt, Func<real> *u, Func<real> *v)
{
  scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->curl[i] * conj(v->curl[i]));
  return result;
}

template<typename real, typename scalar>
scalar int_v1(int n, double *wt, Func<real> *v)
{
  scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (v->val1[i]);
  return result;
}

template<typename real, typename scalar>
scalar int_F_e_f(int n, double *wt, double (*F)(int marker, real x, real y), Func<real> *u, Func<real> *v, Geom<real> *e)
{
  scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (*F)(e->elem_marker, e->x[i], e->y[i]) * (u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
  return result;
}

template<typename real, typename scalar>
scalar int_e_tau_f_tau(int n, double *wt, Func<real> *u, Func<real> *v, Geom<real> *e)
{
  scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (    (u->val0[i] * e->tx[i] + u->val1[i] * e->ty[i]) *
                       conj(v->val0[i] * e->tx[i] + v->val1[i] * e->ty[i]));
  return result;
}

#endif

