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

template<typename Real, typename Scalar>
Scalar int_e_f(int n, double *wt, Func<Real> *u, Func<Real> *v)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
  return result;
}

class MatrixFormVolHCurl : public WeakForm::MatrixFormVol
{
public:
    MatrixFormVolHCurl(int i, int j, SymFlag sym = HERMES_SYM) : MatrixFormVol(i, j, sym) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const
    {
      return int_e_f<Real, Scalar>(n, wt, u, v);
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
Scalar int_curl_e_curl_f(int n, double *wt, Func<Real> *u, Func<Real> *v)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->curl[i] * conj(v->curl[i]));
  return result;
}


template<typename Real, typename Scalar>
Scalar int_v0(int n, double *wt, Func<Scalar> *v)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * v->val0[i];
  return result;
}

template<typename Real, typename Scalar>
Scalar int_v1(int n, double *wt, Func<Real> *v)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * v->val1[i];
  return result;
}

template<typename Real, typename Scalar>
Scalar int_F_e_f(int n, double *wt, double (*F)(int marker, Real x, Real y), Func<Real> *u, Func<Real> *v, Geom<Real> *e)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (*F)(e->elem_marker, e->x[i], e->y[i]) * (u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
  return result;
}

template<typename Real, typename Scalar>
Scalar int_e_tau_f_tau(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (    (u->val0[i] * e->tx[i] + u->val1[i] * e->ty[i]) *
                       conj(v->val0[i] * e->tx[i] + v->val1[i] * e->ty[i]));
  return result;
}

#endif

