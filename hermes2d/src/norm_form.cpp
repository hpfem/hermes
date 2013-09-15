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

#include "norm_form.h"
#include "forms.h"
#include "config.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Real, typename Scalar>
    static Scalar l2_norm(int n, double *wt, Func<Scalar> *u, Func<Real> *v)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * u->val[i] * v->val[i];
      return result;
    }

    template<typename Real, typename Scalar>
    static Scalar h1_norm(int n, double *wt, Func<Scalar> *u, Func<Real> *v)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->val[i] * conj(v->val[i]) + u->dx[i] * conj(v->dx[i]) + u->dy[i] * conj(v->dy[i]));
      return result;
    }

    template<typename Real, typename Scalar>
    static Scalar h1_seminorm(int n, double *wt, Func<Scalar> *u, Func<Real> *v)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->val[i] * conj(v->val[i]) + u->dx[i] * conj(v->dx[i]) + u->dy[i] * conj(v->dy[i]));
      return result;
    }

    template<typename Real, typename Scalar>
    static Scalar hcurl_norm(int n, double *wt, Func<Scalar> *u, Func<Real> *v)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->curl[i] * conj(v->curl[i]) + u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
      return result;
    }

    template<typename Real, typename Scalar>
    static Scalar hdiv_norm(int n, double *wt, Func<Scalar> *u, Func<Real> *v)
    {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->div[i] * conj(v->div[i]) + u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
      return result;
    }

    NormForm::NormForm(int i, int j, FunctionsEvaluatedType functionType) : i(i), j(j), functionType(functionType)
    {
    }

    void NormForm::set_area(std::string area)
    {
      this->area = area;
    }

    template<typename Scalar>
    NormFormVol<Scalar>::NormFormVol(int i, int j) : NormForm(i, j)
    {
    }

    template<typename Scalar>
    NormFormSurf<Scalar>::NormFormSurf(int i, int j) : NormForm(i, j)
    {
    }

    template<typename Scalar>
    NormFormDG<Scalar>::NormFormDG(int i, int j) : NormForm(i, j)
    {
    }
    
    template<typename Scalar>
    DefaultNormFormVol<Scalar>::DefaultNormFormVol(int i, int j, NormType normType) : NormFormVol<Scalar>(i, j), normType(normType)
    {
    }

    template<typename Scalar>
    Scalar DefaultNormFormVol<Scalar>::value(int n, double *wt, Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e) const
    {
      switch(this->normType)
      {
      case HERMES_L2_NORM:
        return l2_norm<Scalar, Scalar>(n, wt, u, v);
      case HERMES_H1_NORM:
        return h1_norm<Scalar, Scalar>(n, wt, u, v);
      case HERMES_H1_SEMINORM:
        return h1_seminorm<Scalar, Scalar>(n, wt, u, v);
      case HERMES_HCURL_NORM:
        return hcurl_norm<Scalar, Scalar>(n, wt, u, v);
      case HERMES_HDIV_NORM:
        return hdiv_norm<Scalar, Scalar>(n, wt, u, v);
      default:
        throw Hermes::Exceptions::Exception("Unknown norm in DefaultNormFormVol<Scalar>::value.");
        return 0.0;
      }
    }
    template<typename Scalar>
    DefaultNormFormSurf<Scalar>::DefaultNormFormSurf(int i, int j, NormType normType) : NormFormSurf<Scalar>(i, j), normType(normType)
    {
    }

    template<typename Scalar>
    Scalar DefaultNormFormSurf<Scalar>::value(int n, double *wt, Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e) const
    {
      switch(this->normType)
      {
      case HERMES_L2_NORM:
        return l2_norm<Scalar, Scalar>(n, wt, u, v);
      case HERMES_H1_NORM:
        return h1_norm<Scalar, Scalar>(n, wt, u, v);
      case HERMES_H1_SEMINORM:
        return h1_seminorm<Scalar, Scalar>(n, wt, u, v);
      case HERMES_HCURL_NORM:
        return hcurl_norm<Scalar, Scalar>(n, wt, u, v);
      case HERMES_HDIV_NORM:
        return hdiv_norm<Scalar, Scalar>(n, wt, u, v);
      default:
        throw Hermes::Exceptions::Exception("Unknown norm in DefaultNormFormSurf<Scalar>::value.");
        return 0.0;
      }
    }


    template<typename Scalar>
    MatrixDefaultNormFormVol<Scalar>::MatrixDefaultNormFormVol(int i, int j, NormType normType) : MatrixFormVol<Scalar>(i, j), normType(normType)
    {
      this->setSymFlag(HERMES_SYM);
    }

    template<typename Scalar>
    Scalar MatrixDefaultNormFormVol<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
      Func<double> *v, Geom<double> *e, Func<Scalar> **ext) const
    {
      switch(this->normType)
      {
      case HERMES_L2_NORM:
        return l2_norm<double, double>(n, wt, u, v);
      case HERMES_H1_NORM:
        return h1_norm<double, double>(n, wt, u, v);
      case HERMES_H1_SEMINORM:
        return h1_seminorm<double, double>(n, wt, u, v);
      case HERMES_HCURL_NORM:
        return hcurl_norm<double, double>(n, wt, u, v);
      case HERMES_HDIV_NORM:
        return hdiv_norm<double, double>(n, wt, u, v);
      default:
        throw Hermes::Exceptions::Exception("Unknown norm in MatrixDefaultNormFormVol<Scalar>::value.");
        return 0.0;
      }
    }

    template<typename Scalar>
    Ord MatrixDefaultNormFormVol<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      Geom<Ord> *e, Func<Ord> **ext) const
    {
      switch(this->normType)
      {
      case HERMES_L2_NORM:
        return l2_norm<Ord, Ord>(n, wt, u, v);
      case HERMES_H1_NORM:
        return h1_norm<Ord, Ord>(n, wt, u, v);
      case HERMES_H1_SEMINORM:
        return h1_seminorm<Ord, Ord>(n, wt, u, v);
      case HERMES_HCURL_NORM:
        return hcurl_norm<Ord, Ord>(n, wt, u, v);
      case HERMES_HDIV_NORM:
        return hdiv_norm<Ord, Ord>(n, wt, u, v);
      default:
        throw Hermes::Exceptions::Exception("Unknown norm in MatrixDefaultNormFormVol<Scalar>::ord.");
        return Ord(0);
      }
    }

    template<typename Scalar>
    MatrixFormVol<Scalar>* MatrixDefaultNormFormVol<Scalar>::clone() const
    {
      return new MatrixDefaultNormFormVol(*this);
    }

    template<typename Scalar>
    VectorDefaultNormFormVol<Scalar>::VectorDefaultNormFormVol(int i, NormType normType) : VectorFormVol<Scalar>(i), normType(normType)
    {
    }

    template<typename Scalar>
    Scalar VectorDefaultNormFormVol<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v, Geom<double> *e, Func<Scalar> **ext) const
    {
      switch(this->normType)
      {
      case HERMES_L2_NORM:
        return l2_norm<double, Scalar>(n, wt, ext[0], v);
      case HERMES_H1_NORM:
        return h1_norm<double, Scalar>(n, wt, ext[0], v);
      case HERMES_H1_SEMINORM:
        return h1_seminorm<double, Scalar>(n, wt, ext[0], v);
      case HERMES_HCURL_NORM:
        return hcurl_norm<double, Scalar>(n, wt, ext[0], v);
      case HERMES_HDIV_NORM:
        return hdiv_norm<double, Scalar>(n, wt, ext[0], v);
      default:
        throw Hermes::Exceptions::Exception("Unknown norm in VectorDefaultNormFormVol<Scalar>::value.");
        return 0.0;
      }
    }

    template<typename Scalar>
    Ord VectorDefaultNormFormVol<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
      Geom<Ord> *e, Func<Ord> **ext) const
    {
      switch(this->normType)
      {
      case HERMES_L2_NORM:
        return l2_norm<Ord, Ord>(n, wt, ext[0], v);
      case HERMES_H1_NORM:
        return h1_norm<Ord, Ord>(n, wt, ext[0], v);
      case HERMES_H1_SEMINORM:
        return h1_seminorm<Ord, Ord>(n, wt, ext[0], v);
      case HERMES_HCURL_NORM:
        return hcurl_norm<Ord, Ord>(n, wt, ext[0], v);
      case HERMES_HDIV_NORM:
        return hdiv_norm<Ord, Ord>(n, wt, ext[0], v);
      default:
        throw Hermes::Exceptions::Exception("Unknown norm in VectorDefaultNormFormVol<Scalar>::ord.");
        return Ord(0);
      }
    }

    template<typename Scalar>
    VectorFormVol<Scalar>* VectorDefaultNormFormVol<Scalar>::clone() const
    {
      return new VectorDefaultNormFormVol(*this);
    }

    template HERMES_API class NormFormVol<double>;
    template HERMES_API class NormFormVol<std::complex<double> >;
    template HERMES_API class NormFormSurf<double>;
    template HERMES_API class NormFormSurf<std::complex<double> >;
    template HERMES_API class NormFormDG<double>;
    template HERMES_API class NormFormDG<std::complex<double> >;
    template HERMES_API class DefaultNormFormVol<double>;
    template HERMES_API class DefaultNormFormVol<std::complex<double> >;
    template HERMES_API class DefaultNormFormSurf<double>;
    template HERMES_API class DefaultNormFormSurf<std::complex<double> >;
    template HERMES_API class MatrixDefaultNormFormVol<double>;
    template HERMES_API class MatrixDefaultNormFormVol<std::complex<double> >;
    template HERMES_API class VectorDefaultNormFormVol<double>;
    template HERMES_API class VectorDefaultNormFormVol<std::complex<double> >;
  }
}