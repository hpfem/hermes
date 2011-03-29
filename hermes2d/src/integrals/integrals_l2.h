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
class MatrixFormVolL2 : public WeakForm::MatrixFormVol
{
public:
    MatrixFormVolL2(int i, int j, SymFlag sym = HERMES_SYM) : MatrixFormVol(i, j, sym) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const
    {
      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->val[i] * conj(v->val[i]));
      return result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                 Geom<double> *e, ExtData<scalar> *ext) const
    {
        return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
            Geom<Ord> *e, ExtData<Ord> *ext) const
    {
        return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
};

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
#endif
