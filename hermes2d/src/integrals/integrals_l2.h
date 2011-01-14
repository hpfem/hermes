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

#ifndef __H2D_INTEGRALS_L2_H
#define __H2D_INTEGRALS_L2_H

#include "../quadrature/limit_order.h"
#include "../weakform/weakform.h"

//// some l2 integrals ////

template<typename Real, typename Scalar>
Scalar l2_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
               Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * conj(v->val[i]));
  return result;
}

template<typename Real, typename Scalar>
Scalar l2_residual_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                        Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* u_prev = u_ext[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (u_prev->val[i] * conj(v->val[i]));
  return result;
}

//// error calculation for adaptivity  ////

template<typename Real, typename Scalar>
Scalar l2_error_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Scalar> *u, 
               Func<Scalar> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * conj(v->val[i]));
  return result;
}




#endif
