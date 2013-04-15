// This file is part of Hermes2D
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, see <http://www.gnu.prg/licenses/>.

#include "norm_calculator.h"
#include "config.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename FormType, typename Scalar>
    Scalar NormCalculator::NormForm<FormType, Func<Scalar>, HERMES_L2_NORM>::value(int n, double *wt, Func<Scalar> *u_ext[],
      Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e,
      Func<Scalar> **ext) const
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * u->val[i] * v->val[i];
      return result;
    }

    template<typename FormType, typename Scalar>
    Hermes::Ord NormCalculator::NormForm<FormType, Func<Scalar>, HERMES_L2_NORM>::ord(int n, double *wt, Func<Hermes::Ord> *u_ext[],
      Func<Hermes::Ord> *u, Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e,
      Func<Ord> **ext) const
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * u->val[i] * v->val[i];
      return result;
    }

    template<typename FormType, typename Scalar>
    Scalar NormCalculator::NormForm<FormType, Func<Scalar>, HERMES_H1_NORM>::value(int n, double *wt, Func<Scalar> *u_ext[],
      Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e,
      Func<Scalar> **ext) const
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->val[i] * conj(v->val[i]) + u->dx[i] * conj(v->dx[i]) + u->dy[i] * conj(v->dy[i]));
      return result;
    }

    template<typename FormType, typename Scalar>
    Hermes::Ord NormCalculator::NormForm<FormType, Func<Scalar>, HERMES_H1_NORM>::ord(int n, double *wt, Func<Hermes::Ord> *u_ext[],
      Func<Hermes::Ord> *u, Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e,
      Func<Ord> **ext) const
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->val[i] * conj(v->val[i]) + u->dx[i] * conj(v->dx[i]) + u->dy[i] * conj(v->dy[i]));
      return result;
    }

    template<typename FormType, typename Scalar>
    Scalar NormCalculator::NormForm<FormType, Func<Scalar>, HERMES_H1_SEMINORM>::value(int n, double *wt, Func<Scalar> *u_ext[],
      Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e,
      Func<Scalar> **ext) const
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->dx[i] * conj(v->dx[i]) + u->dy[i] * conj(v->dy[i]));
      return result;
    }

    template<typename FormType, typename Scalar>
    Hermes::Ord NormCalculator::NormForm<FormType, Func<Scalar>, HERMES_H1_SEMINORM>::ord(int n, double *wt, Func<Hermes::Ord> *u_ext[],
      Func<Hermes::Ord> *u, Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e,
      Func<Ord> **ext) const
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->dx[i] * conj(v->dx[i]) + u->dy[i] * conj(v->dy[i]));
      return result;
    }

    template<typename FormType, typename Scalar>
    Scalar NormCalculator::NormForm<FormType, Func<Scalar>, HERMES_HCURL_NORM>::value(int n, double *wt, Func<Scalar> *u_ext[],
      Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e,
      Func<Scalar> **ext) const
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->curl[i] * conj(v->curl[i]) + u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
      return result;
    }

    template<typename FormType, typename Scalar>
    Hermes::Ord NormCalculator::NormForm<FormType, Func<Scalar>, HERMES_HCURL_NORM>::ord(int n, double *wt, Func<Hermes::Ord> *u_ext[],
      Func<Hermes::Ord> *u, Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e,
      Func<Ord> **ext) const
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->curl[i] * conj(v->curl[i]) + u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
      return result;
    }

    template<typename FormType, typename Scalar>
    Scalar NormCalculator::NormForm<FormType, Func<Scalar>, HERMES_HDIV_NORM>::value(int n, double *wt, Func<Scalar> *u_ext[],
      Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e,
      Func<Scalar> **ext) const
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->curl[i] * conj(v->curl[i]) + u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
      return result;
    }

    template<typename FormType, typename Scalar>
    Hermes::Ord NormCalculator::NormForm<FormType, Func<Scalar>, HERMES_HDIV_NORM>::ord(int n, double *wt, Func<Hermes::Ord> *u_ext[],
      Func<Hermes::Ord> *u, Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e,
      Func<Ord> **ext) const
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->curl[i] * conj(v->curl[i]) + u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
      return result;
    }
  }
}