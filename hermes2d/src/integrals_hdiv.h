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

#ifndef __H2D_INTEGRALS_HDIV_H
#define __H2D_INTEGRALS_HDIV_H

// a
//// volume integrals //////////////////////////////////////////////////////////////////////////////

template<typename T>
inline T int_g_h(Function<T>* fg, Function<T>* fh, RefMap* rg, RefMap* rh)
{
  Quad2D* quad = rg->get_quad_2d();
  int o = fg->get_fn_order() + fh->get_fn_order() + 2 + rg->get_inv_ref_order();
  limit_order(o);
  fg->set_quad_order(o, H2D_FN_VAL);
  fh->set_quad_order(o, H2D_FN_VAL);

  T *g0 = fg->get_fn_values(0), *g1 = fg->get_fn_values(1);
  T *h0 = fh->get_fn_values(0), *h1 = fh->get_fn_values(1);

  T result = 0.0;

  double3* pt = quad->get_points(o);
  int np = quad->get_num_points(o);
  double2x2* mg, *mh;
  mg = rg->get_inv_ref_map(o);
  mh = rh->get_inv_ref_map(o);
  double* jac = rg->get_jacobian(o);
  for (int i = 0; i < np; i++, mg++, mh++)
      result += pt[i][2] * jac[i] * (( (*mg)[1][1]*g0[i] - (*mg)[1][0]*g1[i]) * ( (*mh)[1][1]*h0[i] - (*mh)[1][0]*h1[i]) +
                                     (-(*mg)[0][1]*g0[i] + (*mg)[0][0]*g1[i]) * (-(*mh)[0][1]*h0[i] + (*mh)[0][0]*h1[i]));

  return result;
}

//// error calculation for adaptivity  //////////////////////////////////////////////////////////////////////////////

template<typename Real, typename Scalar>
Scalar hdiv_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Scalar> *u, Func<Scalar> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{

  error("hdiv error form not implemented yet in integrals_hdiv.h.");

  // this is Hcurl code:
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->curl[i] * conj(v->curl[i]) +
                       u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
  return result;
}


#endif
